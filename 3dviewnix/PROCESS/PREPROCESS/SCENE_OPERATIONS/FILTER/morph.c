/*
  Copyright 1993-2015 Medical Image Processing Group
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

/*******************
*  In this file:   *
*  TAB = 4 spaces  *
* (Set you editor) *
*******************/


/*
	This program performs filtering on an nD Scene (3DViewnix format).

	The type of filter applied is:



	Each slice is is loaded and filtered individually.



Author: Roberto J. Goncalves
Date  : 11/02/92


Modified for morphological operations by Dewey Odhner
Date: 11/10/00

*/


#include <math.h>
 
#include <cv3dv.h>
 
 

#include "slices.c"
#include "fff.c"



/*-------------------------------------------------------------------------*/
/* Read the file header and its length */
int get_file_info(file, vh, len)
ViewnixHeader *vh;
char *file; /* filename */
int *len;   /* length of header */
{
    char group[5], element[5];
    FILE *fp;
    int error;
 
    /* open file */
    if ((fp=fopen(file,"rb"))==NULL)
    {
        fprintf(stderr,"ERROR: Can't open file [%s] !\n", file);
        VDeleteBackgroundProcessInformation();
        exit(0);
    }
 
 
    error = VReadHeader(fp, vh, group, element);
    if( error>0 && error < 106)
    {
        fprintf(stderr,"ERROR #%d: Inconsistent 3DVIEWNIX header on file [%s] !\n", error, file);
        fprintf(stderr,"group=%s,  element=%s\n", group, element);
        VDeleteBackgroundProcessInformation();
        exit(0);
    }
 
    error = VGetHeaderLength(fp, len);
    if( error>0 )
    {
        fprintf(stderr,"ERROR #%d: Can't get header length on file [%s] !\n", error, file);
        VDeleteBackgroundProcessInformation();
        exit(0);
    }
 
    fclose(fp);
 
    return(0);
 
}



/* Modified: 11/14/00 for morphological operations by Dewey Odhner */
int main(argc, argv)
int argc;
char *argv[];
{
	FILE *fpin, *fpout;	/* inpput/output files */
	int execution_mode;	/* execution mode */
	int operation;		/* execution mode */
	ViewnixHeader vh;	/* 3DViewnix header */
	SLICES	sl;			/* Structure containing information about the slices of the scene */
	char group[5],		/* Used in VWriteHeader */
		element[5];
	int length;			/* length of a slice */
	int hlength;		/* length of the input header */
	int width, height;	/* dimensions of a slice */
	int nbits;			/* Number of bits of input data */
	int i,j,k;			/* general use */
	int error;			/* error code */
	char *comments;     /* used to modify the header (description field) */
	float space,pixel;
	int ll, m;

	unsigned char *in_buffer1;
	unsigned char *in_buffer8a, *in_buffer8b, *in_buffer8c, *out_buffer8, *temp8;
	unsigned short *in_buffer16a, *in_buffer16b, *in_buffer16c, *out_buffer16, *temp16;


	if(argc != 5)
	{
		printf("Usage:\n");
		printf("%s input output operation mode\n", argv[0]);
		printf("where:\n");
		printf("input    : name of input file;\n");
		printf("output   : name of output file;\n");
		printf("operation: [+|-][5|7|9|19|27];\n");
		printf("mode     : mode of operation (0=foreground, 1=background);\n");
		exit(1);
	}
	


	 
 
    /* Open INPUT and OUTPUT Files */
    if( (fpin = fopen(argv[1], "rb")) == NULL)
    {
        printf("ERROR: Can't open INPUT file !\n");
        exit(1);
    }
    if( (fpout = fopen(argv[2], "wb")) == NULL)
    {
        printf("ERROR: Can't open OUTPUT file !\n");
        exit(1);
    }

    /* Get operation */
    sscanf(argv[3], "%d", &operation);

    /* Get EXECUTION MODE */
    sscanf(argv[4], "%d", &execution_mode);


	/* If in background mode, then place an entry in the BG_PROCESS.COM file */
    if(execution_mode == 1)
        VAddBackgroundProcessInformation("ndinterpolate");


 
    /*-----------------------*/
    /* Read 3DViewnix header */
    /*-----------------------*/
    get_file_info(argv[1], &vh, &hlength);
 

	/* Comoute information about the slices of the scene (number, locations, etc...) */
	compute_slices(&vh, &sl);


	/* Calculate length of a slice */
	width =  vh.scn.xysize[0];
	height =  vh.scn.xysize[1];
	nbits = vh.scn.num_of_bits;
	length = width * height * nbits / 8;
	if(nbits == 1)
		length = (width * height + 7) / 8;
	pixel = vh.scn.xypixsz[0];
	space = sl.Min_spacing3;


	/* Allocate memory */
	if(nbits == 1)
	{
    	/* create buffer for one binary image */
    	if( (in_buffer1 = (unsigned char *) calloc(1, length) ) == NULL)
    	{
       		printf("ERROR: Can't allocate input image buffer.\n");
       		exit(1);
    	}

		length = length*8;

    	/* create buffer for one grey image */
    	if( (in_buffer8b = (unsigned char *) calloc(1, length) ) == NULL  ||
    		(in_buffer8a = (unsigned char *) calloc(1, length) ) == NULL  ||
    		(in_buffer8c = (unsigned char *) calloc(1, length) ) == NULL  ||
    		(out_buffer8 = (unsigned char *) calloc(1, length) ) == NULL)
    	{
       		printf("ERROR: Can't allocate input image buffer.\n");
       		exit(1);
    	}
 
	}
	else
	if(nbits == 8)
	{
    	/* create buffer for one grey image */
    	if( (in_buffer8b = (unsigned char *) calloc(1, length) ) == NULL  ||
    		(in_buffer8a = (unsigned char *) calloc(1, length) ) == NULL  ||
    		(in_buffer8c = (unsigned char *) calloc(1, length) ) == NULL  ||
    		(out_buffer8 = (unsigned char *) calloc(1, length) ) == NULL)
    	{
       		printf("ERROR: Can't allocate input image buffer.\n");
       		exit(1);
    	}
 
	}
	else
	{
    	/* create buffer for one grey image */
    	if( (in_buffer16b = (unsigned short *) calloc(1, length) ) == NULL  ||
    		(in_buffer16a = (unsigned short *) calloc(1, length) ) == NULL  ||
    		(in_buffer16c = (unsigned short *) calloc(1, length) ) == NULL  ||
    		(out_buffer16 = (unsigned short *) calloc(1, length) ) == NULL)
    	{
       		printf("ERROR: Can't allocate input image buffer.\n");
       		exit(1);
    	}
 
	}

	/*-------------------------*/
	/* Modify 3DViewnix Header */
    strncpy(vh.gen.filename1, argv[1], sizeof(vh.gen.filename1)-1);
	vh.gen.filename1[sizeof(vh.gen.filename1)-1] = 0;
    strncpy(vh.gen.filename, argv[2], sizeof(vh.gen.filename)-1);
	vh.gen.filename[sizeof(vh.gen.filename)-1] = 0;
    /* Build "description" header entry */
	comments = (char *)calloc(1,1);
    for(i=0; i<argc; i++)
    {
        comments=(char *)realloc(comments, strlen(comments)+strlen(argv[i])+2);
		strcat(comments,argv[i]);
        strcat(comments," ");
    }
    vh.scn.description = (char *) malloc( strlen(comments) + 1);
    strcpy(vh.scn.description, comments);
    vh.scn.description_valid = 0x1;

	/* Write output 3DViewnix Header */
	error = VWriteHeader(fpout, &vh, group, element);
	if(error < 100)
	{
		printf("ERROR: Can't write output Header (ERROR #%d, %s-%s)!\n", error,group, element);
		exit(error);
	}


	/************************/	
	/* Traverse ALL VOLUMES */
	k=0;	/* k=index of the current slice (within all slices) */
	for(j=0; j<sl.volumes; j++)
	{
		/* Seek the appropriate location */
		if(nbits == 1)
			fseek(fpin, (k*length/8)+hlength, 0);
		else
			fseek(fpin, k*length+hlength, 0);


	   /* Read slices for filtering first slice on the scene (a=b=1st, c=2nd slice) */
	   if(nbits == 1)
	   {
		   /* read first slice */
		   if(fread(in_buffer1, 1, (length/8), fpin) != length/8)
       	   {
       		   printf("ERROR: Couldn't read slice #%d of volume #%d.\n", i+1, j+1);
       		   exit(2);
       	   }
		   bin_to_grey(in_buffer1, (length/8), in_buffer8b, 0, 255);

		   for(i=length-1; i>=0; i--)
			   in_buffer8a[i] = 0;
			
	   }
	   /* 8 BITS/PIXEL */
	   else
	   if (nbits <= 8)
	   {
		   if(fread(in_buffer8b, 1, length, fpin) != length)
       	   {
       		   printf("ERROR: Couldn't read slice #%d of volume #%d.\n", i+1, j+1);
       		   exit(2);
       	   }

		   for(i=length-1; i>=0; i--)
			   in_buffer8a[i] = 0;
	   }
	   /* 16 BITS/PIXEL */
	   else
	   {
			VReadData((char *)in_buffer16b, 2, length/2, fpin, &ll);
			if( ll != length/2)
       	   	{
       		   printf("ERROR: Couldn't read slice #%d of volume #%d.\n", i+1, j+1);
       		   exit(2);
       	   	}
		   	for(i=(length/2)-1; i>=0; i--)
			   in_buffer16a[i] = 0;
	   }
	
	  /*--------------------------------------*/
	  /* For each Volume, traverse ALL SLICES */
	  for(i=0; i<sl.slices[j]; i++)
	  {


		if(execution_mode == 0)
		{
			if(sl.volumes > 1)
			printf("Filtering volume #%d/%d, slice #%d/%d ...\n", j+1,sl.volumes,i+1, sl.slices[j]);
			else
			printf("Filtering slice #%d/%d ...\n", i+1, sl.slices[j]);

			fflush(stdout);
		}


		/*--------------------*/
		/* LOAD NEXT SLICE IN */

		if(i>0)
		{
			if(nbits <= 8)
			{
				temp8 = in_buffer8a;
				in_buffer8a = in_buffer8b;
				in_buffer8b = in_buffer8c;
				in_buffer8c = temp8;
			}
			else
			{
				temp16 = in_buffer16a;
				in_buffer16a = in_buffer16b;
				in_buffer16b = in_buffer16c;
				in_buffer16c = temp16;
			}
		}


		/* 1 or 8 Bits/Pixel */
		if(nbits <= 8)
		{
			/* Load input slice */
			/* BINARY */
			if(nbits == 1)
			{
				/* 3rd slice */
				if (i < sl.slices[j]-1)
				{
				  if (fread(in_buffer1, 1, (length/8), fpin) != length/8)
       			  {
       				printf("ERROR: Couldn't read slice #%d of volume #%d.\n",
					  i+1, j+1);
       				exit(2);
       			  }
				}
				else
				  memset(in_buffer1, 0, length/8);
				bin_to_grey(in_buffer1, (length/8), in_buffer8c, 0, 255);
			}
			/* 8 BITS/PIXEL */
			else
			{
				/* 3rd slice */
				if (i < sl.slices[j]-1)
				{
				  if (fread(in_buffer8c, 1, length, fpin) != length)
       			  {
       				printf("ERROR: Couldn't read slice #%d of volume #%d.\n", i+1, j+1);
       				exit(2);
       			  }
				}
				else
				  memset(in_buffer8c, 0, length);
			}

			/* Filter */
			morph_8(in_buffer8a, in_buffer8b, in_buffer8c, width,height,
								operation, out_buffer8);

			/* Save output slice */
        	if (nbits == 1)
			{
				for (ll=0; ll<width*height/8; ll++)
					out_buffer8[ll] =
						(out_buffer8[ll*8]? 128:0) |
						(out_buffer8[ll*8+1]? 64:0) |
						(out_buffer8[ll*8+2]? 32:0) |
						(out_buffer8[ll*8+3]? 16:0) |
						(out_buffer8[ll*8+4]? 8:0) |
						(out_buffer8[ll*8+5]? 4:0) |
						(out_buffer8[ll*8+6]? 2:0) |
						(out_buffer8[ll*8+7]? 1:0);
				if (width*height%8)
				{
					m = 0;
					for (ll=0; ll<width*height%8; ll++)
						if (out_buffer8[width*height/8*8+ll])
							m |= 128>>ll;
					out_buffer8[width*height/8] = m;
				}
				if (fwrite(out_buffer8, 1, (width*height+7)/8, fpout) !=
						(width*height+7)/8)
        		{
           			printf("ERROR: Couldn't write slice #%d of volume #%d.\n",
						i+1, j+1);
           			exit(3);
        		}
			}
			else if(fwrite(out_buffer8, 1,width*height, fpout) != width*height)
        	{
           		printf("ERROR: Couldn't write slice #%d of volume #%d.\n",
					i+1, j+1);
           		exit(3);
        	}
		}
		/* 16 Bits/Pixel */
		else
		{
			/* 3rd slice */
			if( i < sl.slices[j] - 1)
			{
				VReadData((char *)in_buffer16c, 2, length/2, fpin, &ll);
				if( ll != length/2)
       			{
       				printf("ERROR: Couldn't read slice #%d of volume #%d.\n",
						i+1, j+1);
       				exit(2);
       			}
			}
			else
				memset(in_buffer16c, 0, length);

			/* Filter */
			morph_16(in_buffer16a, in_buffer16b, in_buffer16c,
								width, height, operation, out_buffer16);

			/* Save output slice */
			VWriteData((char *)out_buffer16, 2, width*height, fpout, &ll);
			if(ll != length/2)
        	{
           		printf("ERROR: Couldn't write slice #%d of volume #%d.\n", i+1, j+1);
           		exit(3);
        	}
		}

		k++;

	  } /* end for-loop for all slices[i] */
	} /* end for-loop for volumes[j] */

	VCloseData(fpout);

	if(execution_mode == 0)
	{
		printf("Done.\n");
		fflush(stdout);
	}

 
	/* If in Background mode remove the entry from the BG_PROCESS.COM */
    if(execution_mode == 1)
        VDeleteBackgroundProcessInformation();


	exit(0);
}
