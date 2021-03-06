/*
  Copyright 1993-2014 Medical Image Processing Group
              Department of Radiology
            University of Pennsylvania

This file is part of CAVASS.

CAVASS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

CAVASS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with CAVASS.  If not, see <http://www.gnu.org/licenses/>.

*/

/*****************************************************************************
 * Title: scale_based_filtering_2D
 * Created: 1998 by Punam K Saha
 * Modified: February, 2004 by Punam K Saha
 * Function: Spherical-scale-based anisotropic fltering in two dimentions
 * Assumptions: All parameters are properly specified and the input file is readable
 * Parameters: homogeneity
 * Other effects: none
 * Output: FilteredFile ScaleFile 
 *****************************************************************************/

/***************************************************************/
#include <stdio.h>
#include <math.h>
#include <cv3dv.h>
#define DEFAULT_MAX_SCALE 20 /* Minimum value 2 i.e., without scale */
#define MIN_SCALE 1        /* Minimum value 1 */
#define MIN_THR -1
#define OBJECT_FRACTION_THRESHOLD 25.0
#define MEDIAN_LEVEL 4

#define Handle_error(message) \
{ \
  printf(message); \
  fflush(stdout); \
  exit(1); \
}

int compute_filter_image(int), compute_feature_scale(), get_max(int, int),
 get_min(int, int);

int             max_scale=DEFAULT_MAX_SCALE;
ViewnixHeader   vh;
double          sigma;
unsigned short **pimage_16;
unsigned char **pimage_8;
int             cslice, prow, pcol;
unsigned short *feature_scale, **pfeature_scale;
int            *circle_no_points, (**circle_points)[2];
unsigned short **pfeature_data_filter;
unsigned short *feature_data_filter;
unsigned short *slice_data_16;
unsigned char **pfeature_data_filt_8, *feature_data_filt_8, *slice_data_8;
double          transformation_scale[65536];
long            histogram_count;
int             size, bytes;
char           *scale_file_name;
/*******************************************************************/
/*    Modified: 2/27/04 8-bit scenes handled by Dewey Odhner. */
int main(argc, argv)
    int             argc;
    char           *argv[];
{
    int             i, j, k, tti1, tti2, error, **pptti1;
    FILE           *in, *out;
    char            group[6], elem[6];
    double          inv_sigma;
    double          tt1, circle_row_coefficient, circle_column_coefficient;


	if (argc==6 && sscanf(argv[5], "%d", &max_scale)==1 && max_scale>=2)
		argc = 5;
    if (argc != 5 || sscanf(argv[4], "%lf", &sigma) != 1) {
        printf("Usage: %s image_file output_file scale_file homogeneity [max_scale]\n", argv[0]);
        exit(-1);
    }

    scale_file_name = argv[3];
    in = fopen(argv[1], "rb");
    if (in == NULL)
        Handle_error("Couldn't open input file\n");

    error = VReadHeader(in, &vh, group, elem);
    if (error > 0 && error <= 104)
        Handle_error("Fatal error in reading header\n");

    if (vh.gen.data_type != IMAGE0)
        Handle_error("This is not an IMAGE0 file\n");

    bytes = vh.scn.num_of_bits / 8;
    if (vh.scn.num_of_bits!=16 && vh.scn.num_of_bits!=8)
        Handle_error("This program only handles 8- and 16-bit images\n");

    if (!vh.scn.num_of_subscenes_valid)
        Handle_error("Number of subscenes is not valid!\n");

    if (vh.scn.dimension != 3)
        Handle_error("This program handles only 3D scenes!\n");

    if (!vh.scn.xypixsz_valid)
        Handle_error("Resolution of this scenes is not valid!\n");

    pcol = vh.scn.xysize[0];
    prow = vh.scn.xysize[1];
    size = vh.scn.xysize[0] * vh.scn.xysize[1];

    /******* Computation of transformation functions ************************/

    inv_sigma = -0.5 / pow(sigma, 2.0);

    for (i = 0; i <= 65535; i++)
        transformation_scale[i] = exp(inv_sigma *
                          pow((double) i, 2.0));


    /***********************************************/
	if (bytes == 2)
	{
	    slice_data_16 = (unsigned short *) malloc(size * bytes);
	    if (slice_data_16 == NULL)
	        Handle_error("Couldn't allocate memory (execution terminated)\n");

	    pimage_16 = (unsigned short **) malloc(vh.scn.xysize[1] *
                           sizeof(unsigned short *));
	    if (pimage_16 == NULL)
	        Handle_error("Couldn't allocate memory (execution terminated)\n");

    	for (j = 0; j < vh.scn.xysize[1]; j++)
	        pimage_16[j] = slice_data_16 + j * vh.scn.xysize[0];
	}
	else
	{
	    slice_data_8 = (unsigned char *) malloc(size * bytes);
	    if (slice_data_8 == NULL)
	        Handle_error("Couldn't allocate memory (execution terminated)\n");

	    pimage_8 = (unsigned char **) malloc(vh.scn.xysize[1] *
                           sizeof(unsigned char *));
	    if (pimage_8 == NULL)
	        Handle_error("Couldn't allocate memory (execution terminated)\n");

    	for (j = 0; j < vh.scn.xysize[1]; j++)
	        pimage_8[j] = slice_data_8 + j * vh.scn.xysize[0];
	}

    pptti1 = (int **) malloc(2 * (max_scale + 5) * sizeof(int *));
    if (pptti1 == NULL)
        Handle_error("COULD ALLOCATE MEMORY\n");
    pptti1[0] = (int *) malloc(2 * (max_scale + 5) * 2 * (max_scale + 5) * sizeof(int));
    if (pptti1[0] == NULL)
        Handle_error("COULD ALLOCATE MEMORY\n");
    for (i = 0; i < 2 * (max_scale + 5); i++)
        pptti1[i] = pptti1[0] + i * 2 * (max_scale + 5);
    tti1 = max_scale + 5;
    for (i = 0; i < 2 * (max_scale + 5); i++)
        for (j = 0; j < 2 * (max_scale + 5); j++)
            pptti1[i][j] = 0;

    circle_no_points = (int *) malloc(max_scale * sizeof(int));
    if (circle_no_points == NULL)
        Handle_error("Couldn't allocate memory (execution terminated)\n");

    circle_points = (void *) malloc(max_scale * sizeof(void *));
    if (circle_points == NULL)
        Handle_error("Couldn't allocate memory (execution terminated)\n");

    if (vh.scn.xypixsz[0] == vh.scn.xypixsz[1]) {
        circle_column_coefficient = 1.0;
        circle_row_coefficient = 1.0;
    } else if (vh.scn.xypixsz[0] > vh.scn.xypixsz[1]) {
        circle_row_coefficient = 1.0;
        circle_column_coefficient = vh.scn.xypixsz[0] / vh.scn.xypixsz[1];
    } else {
        circle_column_coefficient = 1.0;
        circle_row_coefficient = vh.scn.xypixsz[1] / vh.scn.xypixsz[0];
    }

	printf("Anisotropy: row = %f, column = %f\n",
	        circle_row_coefficient,circle_column_coefficient);

    for (k = 0; k < max_scale; k++) {

        circle_no_points[k] = 0;
        for (i = -k - 2; i <= k + 2; i++)
            for (j = -k - 2; j <= k + 2; j++) 
			  if (pptti1[tti1 + i][tti1 + j]==0)  {
                tt1 = sqrt(pow(((double) i) * circle_row_coefficient, 2.0)
                       + pow(((double) j) * circle_column_coefficient, 2.0));
                if (tt1 <= ((double) k) + 0.5) {
                    circle_no_points[k] = circle_no_points[k] + 1;
                    pptti1[tti1 + i][tti1 + j] = 2;
                }
            }

        circle_points[k] = (void *) malloc(2 * circle_no_points[k] * sizeof(int));
        if (circle_points[k] == NULL)
            Handle_error("Couldn't allocate memory (execution terminated)\n");

        tti2 = 0;
        for (i = -k - 2; i <= k + 2; i++)
            for (j = -k - 2; j <= k + 2; j++)
                if (pptti1[tti1 + i][tti1 + j] == 2) {
                    pptti1[tti1 + i][tti1 + j] = 1;
                    circle_points[k][tti2][0] = i;
                    circle_points[k][tti2][1] = j;
                    tti2 = tti2 + 1;
                }
    }

    /***********************************************/
    feature_data_filter = NULL;
	feature_data_filt_8 = NULL;
    feature_scale = NULL;
    error = VSeekData(in, 0);
	if (error && error <= 104)
	    Handle_error("Seek error\n");
    out = fopen(argv[2], "w+b");
    if (out == NULL)
        Handle_error("Couldn't open output file\n");

    error = VWriteHeader(out, &vh, group, elem);
    if (error && error <= 104)
        Handle_error("Fatal error in writing header\n");
    for (cslice=0; cslice<vh.scn.num_of_subscenes[0]; cslice++)
	{
		if (VReadData(bytes==2? (void *)slice_data_16:(void *)slice_data_8,
			    bytes, size, in, &j)!=0 || j!=size)   {
		    if(j!=0) {
	           printf("Could read %d voxels out of %d\n", j, size);
			   fflush(stdout);
			   exit(-1);
			   }
			else   
	           Handle_error("Could not read data\n");
			}    

	    compute_filter_image(strcmp(argv[2], "/dev/null") == 0);
	
        if (strcmp(argv[2], "/dev/null") == 0)
            continue;
	    /************ WRITE FILTER IMAGE ************/

	    if (VWriteData(bytes==2?(void *)feature_data_filter:(void *)feature_data_filt_8, bytes, size, out, &j)!=0 || j!=size)
	        Handle_error("Could not read data\n");

    }
	fclose(out);

    /***********************************************************/

    printf("\n");

    exit(0);
}
/*****************************************************************************
 * FUNCTION: compute_filter_image
 * DESCRIPTION: Computes the scale-based filtered image and
 *    stores them at feature_data_filter[feature].
 * PARAMETERS:
 *    scale_only: If non-zero, returns after computing scale image.
 * SIDE EFFECTS: Will print working row number on screen while filtering.
 * ENTRY CONDITIONS: (1) Scale image (2D array) is computed,
 *                   (2) Image (2D array) is assigned.
 *                       The variables vh,slice_data, bytes, cslice
 *                       must be properly set.
 *                       If feature_data_filter or feature_data_down is null.
 * RETURN VALUE: None
 * EXIT CONDITIONS: None
 * HISTORY:
 *    Created: 4/19/99  by Punam K. Saha
 *    Modified: 2/27/04 8-bit scenes handled by Dewey Odhner.
 *    Modified: 2/27/04 scale_only flag added by Dewey Odhner.
 *
 *****************************************************************************/
int compute_filter_image(int scale_only)
{
    int             row, col, i, k, tti1, iscale, xx, yy, x, y;
    unsigned short **outptr, **pin_16, **pscale;
	unsigned char **outptr_8, **pin_8;
    double          tt1, count, sum, temp_sum, inv_k, inv_half_scale;


	if (feature_data_filter == NULL && feature_data_filt_8 == NULL)
	{
	  if (bytes == 2)
	  {
	    pfeature_data_filter = (unsigned short **)
	        malloc(prow * sizeof(unsigned short *));
	    if (pfeature_data_filter == NULL)
	        Handle_error("Couldn't allocate memory (execution terminated)\n");

	    feature_data_filter = (unsigned short *) malloc(
	               vh.scn.xysize[0] * vh.scn.xysize[1] * sizeof(short));
	    if (feature_data_filter == NULL)
	        Handle_error("Couldn't allocate memory\n");

	    for (i = 0; i < prow; i++)
	        pfeature_data_filter[i] =
	            feature_data_filter + i * pcol;
	  }
	  else
	  {
	    pfeature_data_filt_8 = (unsigned char **)
	        malloc(prow * sizeof(unsigned char *));
	    if (pfeature_data_filt_8 == NULL)
	        Handle_error("Couldn't allocate memory (execution terminated)\n");

	    feature_data_filt_8 = (unsigned char *) malloc(
	               vh.scn.xysize[0] * vh.scn.xysize[1]);
	    if (feature_data_filt_8 == NULL)
	        Handle_error("Couldn't allocate memory\n");

	    for (i = 0; i < prow; i++)
	        pfeature_data_filt_8[i] =
	            feature_data_filt_8 + i * pcol;
	  }
	}

    printf("\rSCALE (slice): %5d", cslice);
    fflush(stdout);
    compute_feature_scale();
	if (scale_only)
		return (1);
    pscale = pfeature_scale;
    pin_16 = pimage_16;
	pin_8 = pimage_8;

    outptr = pfeature_data_filter;
	outptr_8 = pfeature_data_filt_8;
    printf("\rFILTERING (slice): %5d", cslice);
    fflush(stdout);
    for (row = 0; row < vh.scn.xysize[1]; row++) {
        for (col = 0; col < vh.scn.xysize[0]; col++)
            if ((bytes==2? pin_16[row][col]:pin_8[row][col]) > MIN_THR) {
                sum = 0.0;
                count = 0.00001;
                iscale = pscale[row][col];
                for (k = 0; k < iscale; k++) {
                    temp_sum = 0.0;

                    tt1 = (double) iscale;
                    tt1 = 0.5 * tt1;
                    inv_half_scale = -0.5 / pow(tt1, 2.0);
					inv_half_scale = -0.5 / pow(10.0*tt1,2.0);
                    inv_k = exp(inv_half_scale * pow((double) k, 2.0));

                    for (i = 0; i < circle_no_points[k]; i++) {
                        xx = circle_points[k][i][0];
                        yy = circle_points[k][i][1];
                        x = row + xx;
                        y = col + yy;
                        if (x >= 0 && x < prow && y >= 0 && y < pcol &&
                            (bytes==2? pin_16[x][y]:pin_8[x][y]) > MIN_THR) {
                            tti1 = bytes==2? pin_16[x][y]:pin_8[x][y];

                            temp_sum = temp_sum + tti1;
                            count = count + inv_k;
                        }
                    }
                    sum = sum + temp_sum * inv_k;
                }
				if (bytes == 2)
                	outptr[row][col] = (unsigned short) (sum / count + 0.5);
				else
					outptr_8[row][col] = (unsigned char)(sum / count + 0.5);
            } else
                if (bytes == 2)
					outptr[row][col] = pin_16[row][col];
				else
					outptr_8[row][col] = pin_8[row][col];
    }

    return (1);
}
/*****************************************************************************
 * FUNCTION: compute_feature_scale
 * DESCRIPTION: computes scale representation of the image.
 * PARAMETERS:
 *    feature: whether it is an intensity or an edge or a fractal image.
 * SIDE EFFECTS:
 *    feature_scale_data_valid is set.
 * ENTRY CONDITIONS:
 *    sigma, bytes must be set.
 *    image data is valid.
 * RETURN VALUE: None
 * EXIT CONDITIONS: Undefined if entry condition is not met.
 * HISTORY:
 *    Created: 10/3/97 by Punam K. Saha
 *    Modified: 2/27/04 8-bit scenes handled by Dewey Odhner.
 *
 *****************************************************************************/
int compute_feature_scale()
{
    int             tti2, tti5, error;
    char            group[6], elem[6];
    ViewnixHeader   vh_scale;
    static FILE           *fp_scale;
    short           bit_fields[] = {0, 15};
    float           smallest_density_value_scale[] = {0.0},
                    largest_density_value_scale[] = {(float) max_scale};
    int             row, col, i, j, k, x, y, xx, yy, mean_g, tti1, tli1[9];
    static unsigned short **pin_16, *temp_16, **ptemp_16;
	static unsigned char **pin_8, *temp_8, **ptemp_8;
    int             flag;
    double          count_obj, count_nonobj;

    if (feature_scale == NULL)
	{
	  pfeature_scale = (unsigned short **) malloc(
                       prow * sizeof(unsigned short *));
      if (pfeature_scale == NULL)
        Handle_error("Couldn't allocate memory (execution terminated)\n");

      feature_scale = (unsigned short *) malloc(
                      prow * pcol * sizeof(unsigned short));
      if (feature_scale == NULL)
        Handle_error("Couldn't allocate memory\n");

      for (i = 0; i < prow; i++)
        pfeature_scale[i] = feature_scale + i * pcol;

	  if (bytes == 2)
	  {
	    temp_16= (unsigned short *)malloc(prow* pcol * sizeof(unsigned short));
	    if (temp_16 == NULL)
	        Handle_error("Couldn't allocate memory\n");

	    ptemp_16 = (unsigned short **) malloc(prow * sizeof(unsigned short *));
	    if (ptemp_16 == NULL)
	        Handle_error("Couldn't allocate memory\n");

	    for (i = 0; i < prow; i++)
	        ptemp_16[i] = temp_16 + i * pcol;

	  }
	  else
	  {
	    temp_8 = (unsigned char *) malloc(prow * pcol);
	    if (temp_8 == NULL)
	        Handle_error("Couldn't allocate memory\n");

	    ptemp_8 = (unsigned char **) malloc(prow * sizeof(unsigned char *));
	    if (ptemp_8 == NULL)
	        Handle_error("Couldn't allocate memory\n");

	    for (i = 0; i < prow; i++)
	        ptemp_8[i] = temp_8 + i * pcol;

	  }
	}

	if (bytes == 2)
	    for (row = 0; row < prow; row++)
	        for (col = 0; col < pcol; col++)
	            ptemp_16[row][col] = pimage_16[row][col];
	else
	    for (row = 0; row < prow; row++)
	        for (col = 0; col < pcol; col++)
	            ptemp_8[row][col] = pimage_8[row][col];
    for (row = 0; row < prow; row++)
        for (col = 0; col < pcol; col++) {
            tti1 = 0;
            for (xx = -1; xx <= 1; xx++)
                for (yy = -1; yy <= 1; yy++) {
                    if (xx < 0) {
                        x = get_max(row + xx, 0);
                    } else {
                        x = get_min(row + xx, prow - 1);
                    }
                    if (yy < 0) {
                        y = get_max(col + yy, 0);
                    } else {
                        y = get_min(col + yy, pcol - 1);
                    }
                    j = 0;
					if (bytes == 2)
                    	while (j < tti1 && ptemp_16[x][y] < tli1[j])
                        	j = j + 1;
					else
						while (j < tti1 && ptemp_8[x][y] < tli1[j])
							j = j + 1;
                    tti2 = j;
                    for (j = tti1 - 1; j >= tti2; j--)
                        tli1[j + 1] = tli1[j];
                    tli1[tti2] = bytes==2? ptemp_16[x][y]: ptemp_8[x][y];
                    tti1 = tti1 + 1;
                }
			if (bytes == 2)
			{
	            for (j = 0; j < tti1; j++)
	                if (tli1[j] == pimage_16[row][col]) {
	                    tti2 = j;
	                    tti1 = 9;
	                }
	            if (tti2 < MEDIAN_LEVEL)
	                pimage_16[row][col] = tli1[MEDIAN_LEVEL];
	            else if (tti2 > 8-MEDIAN_LEVEL)
	                pimage_16[row][col] = tli1[8-MEDIAN_LEVEL];
	            else
	                pimage_16[row][col] = ptemp_16[row][col];
			}
			else
			{
	            for (j = 0; j < tti1; j++)
	                if (tli1[j] == pimage_8[row][col]) {
	                    tti2 = j;
	                    tti1 = 9;
	                }
	            if (tti2 < MEDIAN_LEVEL)
	                pimage_8[row][col] = tli1[MEDIAN_LEVEL];
	            else if (tti2 > 8-MEDIAN_LEVEL)
	                pimage_8[row][col] = tli1[8-MEDIAN_LEVEL];
	            else
	                pimage_8[row][col] = ptemp_8[row][col];
			}
        }

    pin_16 = pimage_16;
	pin_8 = pimage_8;
    for (row = 0; row < prow; row++) {
        for (col = 0; col < pcol; col++)
            if ((bytes==2? pimage_16[row][col]:pimage_8[row][col]) <= MIN_THR)
                pfeature_scale[row][col] = MIN_SCALE;
            else {
                flag = 0;
                mean_g = bytes==2? pimage_16[row][col]:pimage_8[row][col];
                count_obj = 0;
                count_nonobj = 0;
                for (k = MIN_SCALE; (k < max_scale) && !flag; k++) {
                    for (i = 0; i < circle_no_points[k]; i++) {
                        xx = row + circle_points[k][i][0];
                        yy = col + circle_points[k][i][1];
                        if (xx >= 0 && xx < prow && yy >= 0 && yy < pcol &&
                            (bytes==2? pin_16[xx][yy]:pin_8[xx][yy]) > MIN_THR)
						{
                            tti5 = (int)(bytes==2?
							    pin_16[xx][yy]: pin_8[xx][yy])- mean_g;
                            if (tti5 < 0)
                                tti5 = -tti5;
                            count_obj = count_obj + transformation_scale[tti5];
                            count_nonobj = count_nonobj + 1.0 - transformation_scale[tti5];
                        }
                    }
                    if (100.0 * count_nonobj >= OBJECT_FRACTION_THRESHOLD
                        * (count_nonobj + count_obj)) {
                        pfeature_scale[row][col] = k;
                        flag = 1;
                    }
                }
                if (!flag)
                    pfeature_scale[row][col] = max_scale - 1;
            }
    }

    /************ WRITE SCALE IMAGE ************/
	if (strcmp(scale_file_name, "/dev/null"))
	{
		if (fp_scale == NULL)
		{
		    vh_scale = vh;
		    vh_scale.scn.num_of_bits = 16;
		    vh_scale.scn.bit_fields = bit_fields;
		    vh_scale.scn.smallest_density_value = smallest_density_value_scale;
		    vh_scale.scn.largest_density_value = largest_density_value_scale;
		    fp_scale = fopen(scale_file_name, "w+b");
		    error = VWriteHeader(fp_scale, &vh_scale, group, elem);
		    if (error <= 104)
		        Handle_error("Fatal error in writing header\n");
		}

	    if (VWriteData((char *)feature_scale, 2, size, fp_scale, &j)!=0 ||
		        j!=size)
	        Handle_error("Could not write data\n");

	}

    /***********************************************************/

    return (1);
}
/*****************************************************************************
   * FUNCTION: get_max
   * DESCRIPTION: returns the maximum between two integers
   * PARAMETERS:
   *    integer1 and integer2
   * SIDE EFFECTS:
   *    Nil.
   * ENTRY CONDITIONS:
   *    Nil.
   * RETURN VALUE: see DESCRIPTION
   * EXIT CONDITIONS:
   * HISTORY:
   *    Created: 9/5/97 by Punam K. Saha
   *
   *****************************************************************************/
int get_max(int a, int b)
{

    if (a > b)
        return (a);
    else
        return (b);
}
/*****************************************************************************
 * FUNCTION: get_min
 * DESCRIPTION: returns the minimum between two integers
 * PARAMETERS:
 *    integer1 and integer2
 * SIDE EFFECTS:
 *    Nil.
 * ENTRY CONDITIONS:
 *    Nil.
 * RETURN VALUE: see DESCRIPTION
 * EXIT CONDITIONS:
 * HISTORY:
 *   Created: 9/5/97 by Punam K. Saha
 *
 *****************************************************************************/
int get_min(int a, int b)
{

    if (a < b)
        return (a);
    else
        return (b);
}
