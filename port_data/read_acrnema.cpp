/*
  Copyright 1993-2014, 2016 Medical Image Processing Group
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


#include  <Viewnix.h>
#include  "cv3dv.h"
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include  "port_data/from_dicom.h"



bool sequence_tag(unsigned tag);

/*-------------------------------------------------------------------------------------*/
/*
 *    Modified: 4/24/95 VReadData called by Dewey Odhner
 *    Modified: 8/25/95 byte swapping cases corrected by Dewey Odhner
 */
int read_BI(FILE *fp, unsigned short *bi, int rec_arch_type)
{
	int items_read;
	union {unsigned short s; char c[2];} t;
	char c;

	if( VReadData((char *)bi, sizeof(unsigned short), 1, fp, &items_read) != 0)
		return(-1);
	switch (rec_arch_type)
	{
		case TYPE2:
		case TYPE3:
		case TYPE6:
		case TYPE7:
			t.s = *bi;
			c = t.c[0];
			t.c[0] = t.c[1];
			t.c[1] = c;
			*bi = t.s;
	}
	return(1);
}

/*-------------------------------------------------------------------------------------*/
/*
 *    Modified: 4/24/95 VReadData called by Dewey Odhner
 */
int read_BD(FILE *fp, unsigned int *bd, int rec_arch_type)
{
	int items_read;
	union {int i; char c[4];} t, u;

	if( VReadData((char *)bd, sizeof(int), 1, fp, &items_read) != 0)
		return(-1);
	switch (rec_arch_type)
	{
		case TYPE1:
		case TYPE5:
			break;
		case TYPE2:
		case TYPE6:
			t.i = *bd;
			u.c[0] = t.c[1];
			u.c[1] = t.c[0];
			u.c[2] = t.c[3];
			u.c[3] = t.c[2];
			*bd = u.i;
			break;
		case TYPE3:
		case TYPE7:
			t.i = *bd;
			u.c[0] = t.c[3];
			u.c[1] = t.c[2];
			u.c[2] = t.c[1];
			u.c[3] = t.c[0];
			*bd = u.i;
			break;
		case TYPE4:
		case TYPE8:
			t.i = *bd;
			u.c[0] = t.c[2];
			u.c[1] = t.c[3];
			u.c[2] = t.c[0];
			u.c[3] = t.c[1];
			*bd = u.i;
			break;
	}
	return(1);
}


/*-------------------------------------------------------------------------------------*/
int read_AN(FILE *fp, char *an, int n, int rec_arch_type)
{
	char temp;
	int j;

	if( fread(an, n, 1, fp) == 0)
	{	an[0] = 0;
		return(-1);
	}
	if (rec_arch_type >= TYPE5)
		for (j=0; j<n; j+=2)
		{	temp = an[j];
			an[j] = an[j+1];
			an[j+1] = temp;
		}
	an[n] = 0;
	return(1);
}



/*-------------------------------------------------------------------------------------*/
int extract_floats(char *text, int n, float floats[])
{
	int i, j;
	int begin, end;
	char cfloat[20];

	strcpy(cfloat, "");

	begin = end = 0;
	for(i=0; i<n; i++)
	{
		while(text[end] != '\\' &&  text[end] != '\0')
			end++;

		j = end - begin;
		strncpy(cfloat, &text[begin], j);
		cfloat[j] = 0;
		sscanf(cfloat, "%f", &(floats[i]));
		begin = end+1;
		end = begin;
	}
	return 0;
}

/*-------------------------------------------------------------------------------------*/
/*
 * RETURN VALUE:
 *    -1: fp is NULL.
 *    0: success
 *    2: read error
 *    5: seek error
 *    100: group or element not found.
 *    235: element length greater than maxlen.
 * HISTORY:
 *    Modified: 4/26/95 element length limited by Dewey Odhner
 *    Modified: 9/5/96 to work without correct group length by George Grevera
 *    Modified: 12/9/98 to work without group length by Dewey Odhner
 *    Modified: 12/11/98 to work with explicit VR encoding by Dewey Odhner
 *    Modified: 12/15/98 to use correct transfer syntax by Dewey Odhner
 *    Modified: 2/13/03 undefined element length/SQ items skipped
 *       by Dewey Odhner.
 *    Modified: 6/16/03 corrected in case transfer syntax UID not given
 *       (no file meta information group) by Dewey Odhner.
 *    Modified: 11/29/05 *items_read initialized to 0 by Dewey Odhner.
 *    Modified: 1/26/10 corrected in case transfer syntax 1.2.840.10008.1.2.2
 *       by Dewey Odhner.
 */
int get_element(FILE *fp, unsigned short group, unsigned short element,
	int type /* type of element being read (BI=16bits, BD=32bits,
	AN=ASCIInumeric, AT=ASCIItext)*/, void *result,
	unsigned int maxlen /* bytes at result */, int *items_read)
{
	int return_val=get_element(fp, group, element, type, result, maxlen,
		items_read, true, 0);
	if (return_val && return_val!=235)
	{
		int return_val2 = get_element(fp, group, element, type, result, maxlen,
			items_read, false, 0);
		if (return_val2==0 || return_val2>return_val)
			return_val = return_val2;
	}
	return return_val;
}


int get_contour_struct_set(FILE *fp, contour_struct_set *ss)
{
	int items_read;
	unsigned short g, e, slen;
	unsigned int len, tag, cbuf_len=0;
	int explicit_vr=1, rec_arch_type=2;
	char *cp, vr[3]={0,0,0}, uid[65], *cbuf;

	ss->num_rois = 0;
	if (fseek(fp, 128, 0)==0 && fread(vr, 1, 4, fp)==4 &&
			strncmp(vr, "DICM", 4)==0)
	{
		if (get_element(fp, 2, 0x10, AT, uid, 64, &items_read, false, 2))
		/* transfer syntax UID not given */
		{
			explicit_vr = 0;
			fseek(fp, 132, 0);
		}
		else
		{
			if (strcmp(uid, "1.2.840.10008.1.2") == 0)
				explicit_vr = 0;
			else if (strcmp(uid, "1.2.840.10008.1.2.2") == 0)
				rec_arch_type = 0;
		}
	}

	// find Structure Set ROI Sequence
	int return_val=get_element(fp, 0x3006, 0x20, BI, 0, 0, &items_read,true,0);
	if (return_val != 235)
		return return_val;
	for (;;)
	{
		/* Check Group */
		if (read_BI(fp, &g, rec_arch_type) < 0)
			return (2);
		/* Check Element */
		if (read_BI(fp, &e, rec_arch_type) < 0)
			return (2);
		if (g==0xfffe && e==0xe0dd) // end of sequence
			break;
		if (g!=0xfffe || e!=0xe000) // beginning of item
			return (100);
		if ((explicit_vr? read_BI(fp, &slen, rec_arch_type):
				read_BD(fp, &len, rec_arch_type)) < 0)
			return (2);
		while (g!=0xfffe || e!=0xe00d) // end of item
		{
			/* Check Group */
			if (read_BI(fp, &g, rec_arch_type) < 0)
				return (2);
			/* Check Element */
			if (read_BI(fp, &e, rec_arch_type) < 0)
				return (2);
			if (explicit_vr)
			{
				if (fread(vr, 1, 2, fp) != 2)
					return (2);
				if (read_BI(fp, &slen, rec_arch_type) < 0)
					return (2);
				len = slen;
			}
			else
				if (read_BD(fp, &len, rec_arch_type) < 0)
					return (2);
			if (g==0x3006 && e==0x22) // ROI number
			{
				if (ss->num_rois == 0)
					ss->roi = (contour_roi *)malloc(sizeof(contour_roi));
				else
					ss->roi = (contour_roi *)realloc(ss->roi,
						(ss->num_rois+1)*sizeof(contour_roi));
				ss->roi[ss->num_rois].roi_name = NULL;
				ss->roi[ss->num_rois].num_contours = 0;
				if (cbuf_len <= len)
				{
					if (cbuf_len)
						free(cbuf);
					cbuf = (char *)malloc(len+1);
					cbuf_len = len+1;
				}
				if (read_AN(fp, cbuf, len, rec_arch_type) < 0)
					return (2);
				ss->roi[ss->num_rois].roi_number = atoi(cbuf);
				ss->num_rois++;
			}
			else if (g==0x3006 && e==0x26) // ROI name
			{
				if (ss->num_rois <= 0)
					return (100);
				if (cbuf_len < len+1)
				{
					if (cbuf_len)
						free(cbuf);
					cbuf = (char *)malloc(len+1);
					cbuf_len = len+1;
				}
				if (read_AN(fp, cbuf, len, rec_arch_type) < 0)
					return (2);
				if (cbuf[len-1] == ' ')
					cbuf[len-1] = 0;
				for (cp=cbuf+1; *cp; cp++)
					if (*cp==' ' || *cp=='/' || *cp=='\\')
						*cp = '_';
				cp = cbuf;
				if (*cp == ' ')
					cp++;
				ss->roi[ss->num_rois-1].roi_name= (char *)malloc(strlen(cp)+1);
				strcpy(ss->roi[ss->num_rois-1].roi_name, cp);
			}
			else
				if (fseek(fp, len, 1))
					return (5);
		}
	}

	// find ROI Contour Sequence
	return_val = get_element(fp, 0x3006, 0x39, BI, 0, 0, &items_read, true, 0);
	if (return_val != 235)
		return return_val;
	int sequence_depth=1;
	int num_contours=0;
	int roi_number=0;
	unsigned int *contour_num_points;
	float (**contour_coord)[3];
	for (;;)
	{
		/* Check Group */
		if (read_BI(fp, &g, rec_arch_type) < 0)
			return (2);
		/* Check Element */
		if (read_BI(fp, &e, rec_arch_type) < 0)
			return (2);
		if (explicit_vr && g!=0xfffe)
		{
			if (fread(vr, 1, 2, fp) != 2)
				return (2);
			if (read_BI(fp, &slen, rec_arch_type) < 0)
				return (2);
			len = slen;
		}
		else
			if (read_BD(fp, &len, rec_arch_type) < 0)
				return (2);
		if (g==0xfffe && e==0xe0dd) // end of sequence
		{
			if (sequence_depth == 1)
				break;
			sequence_depth--;
			continue;
		}
		tag = g;
		tag = (tag<<16)|e;
		if (sequence_tag(tag))
		{
			sequence_depth++;
			continue;
		}
		if (g==0xfffe && e==0xe00d) // end of item
		{
			if (sequence_depth == 1)
			{
				if (num_contours && roi_number)
					for (int j=0; j<ss->num_rois; j++)
						if (roi_number == ss->roi[j].roi_number)
						{
							ss->roi[j].num_contours = num_contours;
							ss->roi[j].contour_num_points = contour_num_points;
							ss->roi[j].contour_coord = contour_coord;
							break;
						}
				roi_number = 0;
				num_contours = 0;
			}
			continue;
		}
		if (g==0xfffe && e==0xe000) // beginning of item
			continue;
		if (g==0x3006 && e==0x46) // Number of Contour Points
		{
			if (num_contours == 0)
			{
				contour_num_points = (unsigned int *)malloc(sizeof(int));
				contour_coord = (float (**)[3])malloc(sizeof(float (*)[3]));
			}
			else
			{
				contour_num_points= (unsigned int *)realloc(contour_num_points,
					(num_contours+1)*sizeof(int));
				contour_coord = (float (**)[3])realloc(contour_coord,
					(num_contours+1)*sizeof(float (*)[3]));
			}
			if (read_AN(fp, cbuf, len, rec_arch_type) < 0)
				return (2);
			contour_num_points[num_contours] = atoi(cbuf);
			contour_coord[num_contours] = (float (*)[3])
				malloc(contour_num_points[num_contours]*sizeof(float [3]));
			if (contour_coord[num_contours] == NULL)
				return (1);
			num_contours++;
			continue;
		}
		if (g==0x3006 && e==0x50) // Contour Data
		{
			if (cbuf_len < 27*contour_num_points[num_contours-1])
			{
				free(cbuf);
				cbuf = (char *)malloc(27*contour_num_points[num_contours-1]);
				cbuf_len = 27*contour_num_points[num_contours-1];
				if (cbuf == NULL)
					return (1);
			}
			if (read_AN(fp, cbuf, len, rec_arch_type) < 0)
				return (2);
			extract_floats(cbuf, 3*contour_num_points[num_contours-1],
				contour_coord[num_contours-1][0]);
			continue;
		}
		if (g==0x3006 && e==0x84) // Referenced ROI Number
		{
			if (read_AN(fp, cbuf, len, rec_arch_type) < 0)
				return (2);
			roi_number = atoi(cbuf);
			continue;
		}
		/* Skip to next Element */
		if (len!=0xffffffff && fseek(fp, len, 1))
			return (5);
	}
	if (cbuf_len)
		free(cbuf);
	return 0;
}

int get_element(FILE *fp, unsigned short group, unsigned short element,
	int type /* type of element being read (BI=16bits, BD=32bits,
	AN=ASCIInumeric, AT=ASCIItext)*/, void *result,
	unsigned int maxlen /* bytes at result */, int *items_read,
	bool skip_sequences, int rec_arch_type)
{
	unsigned short g, e, slen;
	unsigned int len, tag;
	int items, explicit_vr=0, orig_arch;
	char *cp, vr[3]={0,0,0}, uid[65];
	int sequence_depth=0;

	*items_read = 0;
	if(fp == NULL) return (-1);
	orig_arch = rec_arch_type;
	if (fseek(fp, 128, 0)==0 && fread(uid, 1, 4, fp)==4 &&
			strncmp(uid, "DICM", 4)==0)
	{
		explicit_vr = 1;
		rec_arch_type = 2;
		if (group>2 &&
				get_element(fp, 2, 0x10, AT, uid, 64, &items,skip_sequences,2))
		/* transfer syntax UID not given */
		{
			fseek(fp, 132, 0);
			/* Check Group */
			if (read_BI(fp, &g, rec_arch_type) < 0)
				return (2);
			/* Check Element */
			if (read_BI(fp, &e, rec_arch_type) < 0)
				return (2);
			if (g==2 && e==0)
			{
				if (fread(vr, 1, 2, fp) != 2)
					return (2);
				if (strncmp(vr, "UL", 2))
					if (fseek(fp, -2, 1))
						return (5);
				if (read_BD(fp, &len, rec_arch_type) < 0)
					return (2);
				if (fseek(fp, len, 1))
					return (5);
			}
			else
				if (fseek(fp, -4, 1))
					return (5);
			explicit_vr = 0;
			rec_arch_type = orig_arch;
		}
	}
	else
		fseek(fp, 0, 0);

	/* Find Group and Element */
	int ii;
	for (ii=0; ii>=0; ii++)
	{
		/* Check Group */
		if (read_BI(fp, &g, rec_arch_type) < 0)
			return (2);
		/* Check transfer syntax */
		if (explicit_vr && g>2)
		{
			if (strcmp(uid, "1.2.840.10008.1.2") == 0)
				explicit_vr = 0;
			else if (strcmp(uid, "1.2.840.10008.1.2.2") == 0)
			{
				rec_arch_type = orig_arch = 0;
				if (fseek(fp, -2, 1))
					return (5);
				if (read_BI(fp, &g, rec_arch_type) < 0)
					return (2);
			}
		}
		/* Check Element */
		if (read_BI(fp, &e, rec_arch_type) < 0)
			return (2);
		/* Check value representation */
		if (fread(vr, 1, 2, fp) != 2)
			return (2);
		tag = g;
		tag = (tag<<16)|e;
		if (g==0xfffe ||
				(!explicit_vr && (!sequence_tag(tag)||strncmp(vr, "SQ", 2))))
		// Explicit VR not found. (Sometimes sequence VR is given explicitly
		// even when transfer syntax with implicit VR is specified.)
		{
			if (fseek(fp, -2, 1))
				return (5);
			vr[0] = 0;
		}
		if (strncmp(vr, "OB", 2)==0 || strncmp(vr, "OD", 2)==0 ||
				strncmp(vr, "OF", 2)==0 || strncmp(vr, "OL", 2)==0 ||
				strncmp(vr, "OW", 2)==0 || strncmp(vr, "SQ", 2)==0 ||
				strncmp(vr, "UN", 2)==0 || (explicit_vr&&e==0))
		{
			if (fseek(fp, 2, 1))
				return (5);
		}
		/* Check length */
		if (e==0 && vr[0])
			len = 4;
		else if (strncmp(vr, "OB", 2)==0 || strncmp(vr, "OW", 2)==0 ||
				strncmp(vr, "SQ", 2)==0 || strncmp(vr, "UN", 2)==0 || vr[0]==0)
		{
			if (read_BD(fp, &len, rec_arch_type) < 0)
				return (2);
		}
		else
		{
			if (read_BI(fp, &slen, rec_arch_type) < 0)
				return (2);
			len = slen;
		}
		if (sequence_depth==0 && skip_sequences && g>group)
			return (100);
		if (g==group && e==element)
			break;
		if (g==0xfffe && e==0xe0dd) // end of sequence
			sequence_depth--;
		else if (sequence_tag(tag))
			sequence_depth++;
		if (len && len!=0xffffffff && fseek(fp, len, 1))
			return (5);
	}

	/* Read Element */
	switch (type)
	{
		case BI:
			items = (len>maxlen? maxlen:len)/sizeof(short);
			for (*items_read=0; *items_read<items; (*items_read)++)
				if (read_BI(fp, (unsigned short *)result+ *items_read,
						rec_arch_type) < 0)
					return (2);
			if (len > maxlen)
				return (235);
			break;
		case BD:
			items = (len>maxlen? maxlen:len)/sizeof(int);
			for (*items_read=0; *items_read<items; (*items_read)++)
				if (read_BD(fp, (unsigned int *)result+ *items_read,
						rec_arch_type) < 0)
					return (2);
			if (len > maxlen)
				return (235);
			break;
		case AN:
		case AT:
			*items_read = 0;
			if (read_AN(fp, (char *)result, len<maxlen? len: maxlen-1,
					rec_arch_type) < 0)
				return (2);
			for (cp=(char *)result, items=FALSE; ; cp++)
				if (*cp == '\\' || *cp == 0)
				{	if (items)
						(*items_read)++;
					if (*cp == 0)
						break;
					items = FALSE;
				}
				else if (*cp != ' ')
					items = TRUE;
			if (len >= maxlen)
				return (235);
			break;
	}

	return (0);
}


/*****************************************************************************
 * FUNCTION: check_acrnema_file
 * DESCRIPTION: Checks an ACR-NEMA file for image size.
 * PARAMETERS:
 *    path: The directory where the file is
 *    file: The file name
 * SIDE EFFECTS: None
 * ENTRY CONDITIONS: The variable rec_arch_type must be properly set.
 * RETURN VALUE:
 *    0: success
 *    4: Cannot open file
 *    104: Cannot get image size
 * EXIT CONDITIONS: None
 * HISTORY:
 *    Created: 2/16/99 by Dewey Odhner
 *    Modified: 10/18/07 null path allowed by Dewey Odhner.
 *
 *****************************************************************************/
int check_acrnema_file(const char *path, const char *file)
{
    FILE *fp;
    char filename[500];
	int items_read, err;
	unsigned short bi;
 
    if (path)
	{
		strcpy(filename, path);
	    strcat(filename, "/");
	    strcat(filename, file);
	}
	else
		strcpy(filename, file);

    if( (fp = fopen(filename, "rb")) == NULL)
        return (4);
 
	/* IMAGE DIMENSIONS */
	err = get_element(fp, 0x0028, 0x0011, BI,
		&bi, sizeof(bi), &items_read) || bi <= 0 ||
		 get_element(fp, 0x0028, 0x0010, BI,
			&bi, sizeof(bi), &items_read) || bi <= 0;
	fclose(fp);
	return (err? 104: 0);
}

/*****************************************************************************
 * FUNCTION: names_are_similar
 * DESCRIPTION: Checks whether two names differ only by three
 *    consecutive digits or less.
 * PARAMETERS:
 *    name1, name2: the names to compare
 * SIDE EFFECTS: None
 * ENTRY CONDITIONS: None
 * RETURN VALUE: TRUE or FALSE
 * EXIT CONDITIONS: None
 * HISTORY:
 *    Created: 12/1/98 by Dewey Odhner
 *    Modified: 1/18/99 difference not restricted to last three digits
 *       by Dewey Odhner
 *
 *****************************************************************************/
int names_are_similar(const char *name1, const char *name2)
{
	char *e1, *e2;
	char *namep1, *namep2;

	namep1 = (char *)name1;
	namep2 = (char *)name2;
	while (*namep1 && *namep2==*namep1)
		namep1++, namep2++;
	for (e1=namep1+strlen(namep1), e2=namep2+strlen(namep2);
			e1>namep1 && e2>namep2 && e1[-1]==e2[-1]; e1--, e2--)
		;
	if (e1-namep1>3 || e2-namep2>3)
		return (FALSE);
	for (; namep1<e1; namep1++)
		if (!isdigit(*namep1))
			return (FALSE);
	for (; namep2<e2; namep2++)
		if (!isdigit(*namep2))
			return (FALSE);
	return (TRUE);
}

/*****************************************************************************
 * FUNCTION: get_series_uid
 * DESCRIPTION: Returns the series instance UID from a DICOM file.
 * PARAMETERS:
 *    dir: the directory in which the file is
 *    filename: the file name
 * SIDE EFFECTS: An error message may be printed.
 * ENTRY CONDITIONS: rec_arch_type should be set.
 * RETURN VALUE: A string containing the series instance UID, if found, or
 *    NULL.  The memory should be freed by the caller after use.
 * EXIT CONDITIONS: Returns NULL on error without notice.
 * HISTORY:
 *    Created: 2/25/03 by Dewey Odhner
 *    Modified: 4/10/03 dir parameter added by Dewey Odhner.
 *    Modified: 5/19/03 buf_size updated by Dewey Odhner.
 *    Modified: 10/16/07 buf_size updated sooner by Dewey Odhner.
 *    Modified: 10/18/07 null dir allowed by Dewey Odhner.
 *
 *****************************************************************************/
char *get_series_uid(const char dir[], const char filename[])
{
	int items_read, valid, dir_len;
	char *suid;
	static int buf_size;
	static char *buf;
	FILE *fp;

	dir_len = dir? strlen(dir): 0;
	if (buf_size == 0)
	{
		buf = (char *)malloc(dir_len+strlen(filename)+2);
		if (buf == NULL)
			return NULL;
		buf_size = dir_len+strlen(filename)+2;
		buf[buf_size-1] = 0;
	}
	if ((unsigned)buf_size < dir_len+strlen(filename)+2)
	{
		suid = (char *)realloc(buf, dir_len+strlen(filename)+2);
		if (suid == NULL)
			return NULL;
		buf = suid;
		buf_size = dir_len+strlen(filename)+2;
		buf[buf_size-1] = 0;
	}
	if (dir)
		sprintf(buf, "%s/%s", dir, filename);
	else
		strcpy(buf, filename);
	fp = fopen(buf, "rb");
	if (fp == NULL)
	{
		fprintf(stderr, "Can't open file %s\n", filename);
		return NULL;
	}
	valid = get_element(fp, 0x0020, 0x000e, AT, buf, buf_size-1, &items_read);
	while (valid && strlen(buf)+1>=(unsigned)buf_size-1)
	{
		suid = (char *)realloc(buf, buf_size+1000);
		if (suid == NULL)
		{
			fclose(fp);
			return NULL;
		}
		buf = suid;
		buf_size += 1000;
		buf[buf_size-1] = 0;
		valid =
			get_element(fp, 0x0020, 0x000e, AT, buf, buf_size-1, &items_read);
	}
	suid = (char *)malloc(strlen(buf)+1);
	if (suid)
		strcpy(suid, buf);
	fclose(fp);
	return (suid);
}

/*****************************************************************************
 * FUNCTION: get_time
 * DESCRIPTION: Returns the time stamp from a DICOM file.
 * PARAMETERS:
 *    filename: the file name
 * SIDE EFFECTS: An error message may be printed.
 * ENTRY CONDITIONS: rec_arch_type should be set.
 * RETURN VALUE: the time stamp
 * EXIT CONDITIONS:
 * HISTORY:
 *    Created: 9/3/10 by Dewey Odhner
 *
 *****************************************************************************/
float get_time(const char filename[])
{
	int items_read;
	float ft;
	FILE *fp;
	char an[500];

	fp = fopen(filename, "rb");
	if (fp == NULL)
	{
		fprintf(stderr, "Can't open file %s\n", filename);
		return 0.0;
	}
	if (get_element(fp, 0x0008, 0x0033, AN, an, sizeof(an), &items_read) &&
		get_element(fp, 0x0008, 0x0032, AN, an, sizeof(an), &items_read) &&
		get_element(fp, 0x0008, 0x0031, AN, an, sizeof(an), &items_read))
		ft = 0;
	else
		extract_floats(an, 1, &ft);
	fclose(fp);
	return (ft);
}

/*****************************************************************************
 * FUNCTION: get_num_frames
 * DESCRIPTION: Returns the number of frames from a DICOM file.
 * PARAMETERS:
 *    filename: the file name
 * SIDE EFFECTS: An error message may be printed.
 * ENTRY CONDITIONS: rec_arch_type should be set.
 * RETURN VALUE: the number of frames
 * EXIT CONDITIONS:
 * HISTORY:
 *    Created: 4/4/11 by Dewey Odhner
 *
 *****************************************************************************/
int get_num_frames(const char filename[])
{
	int items_read;
	float ft;
	FILE *fp;
	char an[500];

	fp = fopen(filename, "rb");
	if (fp == NULL)
	{
		fprintf(stderr, "Can't open file %s\n", filename);
		return 0;
	}
	if (get_element(fp, 0x0028, 0x0008, AN, an, sizeof(an), &items_read))
		ft = 1;
	else
		extract_floats(an, 1, &ft);
	fclose(fp);
	return (int)ft;
}

bool sequence_tag(unsigned tag)
{
    switch (tag)
    {
        case 0x00041220:
        case 0x00080082:
        case 0x00081032:
        case 0x00081084:
        case 0x00081100:
        case 0x00081110:
        case 0x00081111:
        case 0x00081115:
        case 0x00081120:
        case 0x00081125:
        case 0x00081130:
        case 0x00081140:
        case 0x00081145:
        case 0x0008114A:
        case 0x00081198:
        case 0x00081199:
        case 0x00082112:
        case 0x00082218:
        case 0x00082220:
        case 0x00082228:
        case 0x00082229:
        case 0x00082230:
        case 0x00082240:
        case 0x00082242:
        case 0x00082244:
        case 0x00082246:
        case 0x00100050:
        case 0x00180012:
        case 0x00180014:
        case 0x00180026:
        case 0x00180029:
        case 0x0018002A:
        case 0x00180036:
        case 0x00185104:
        case 0x00186011:
        case 0x00283000:
        case 0x00283010:
        case 0x00283110:
        case 0x00285000:
        case 0x00286100:
        case 0x00321064:
        case 0x00380004:
        case 0x00380044:
        case 0x003A0200:
        case 0x003A0208:
        case 0x003A0209:
        case 0x003A020A:
        case 0x003A0211:
        case 0x00400008:
        case 0x00400100:
        case 0x00400220:
        case 0x00400260:
        case 0x00400270:
        case 0x00400275:
        case 0x00400293:
        case 0x00400295:
        case 0x00400296:
        case 0x00400320:
        case 0x00400321:
        case 0x00400324:
        case 0x00400330:
        case 0x00400340:
        case 0x00400550:
        case 0x00400555:
        case 0x0040059A:
        case 0x0040071A:
        case 0x004008D8:
        case 0x004008DA:
        case 0x004008EA:
        case 0x0040A043:
        case 0x0040A073:
        case 0x0040A088:
        case 0x0040A168:
        case 0x0040A195:
        case 0x0040A300:
        case 0x0040A360:
        case 0x0040A370:
        case 0x0040A372:
        case 0x0040A375:
        case 0x0040A385:
        case 0x0040A504:
        case 0x0040A525:
        case 0x0040A730:
        case 0x0040B020:
        case 0x00500010:
        case 0x00500030:
        case 0x00540012:
        case 0x00540013:
        case 0x00540016:
        case 0x00540022:
        case 0x00540032:
        case 0x00540052:
        case 0x00540062:
        case 0x00540063:
        case 0x00540072:
        case 0x00540220:
        case 0x00540222:
        case 0x00540300:
        case 0x00540302:
        case 0x00540304:
        case 0x00540306:
        case 0x00540410:
        case 0x00540412:
        case 0x00540414:
        case 0x00603000:
        case 0x00700001:
        case 0x00700008:
        case 0x00700009:
        case 0x0070005A:
        case 0x00700060:
        case 0x00880200:
        case 0x2000001E:
        case 0x200000A2:
        case 0x200000A4:
        case 0x200000A8:
        case 0x20000500:
        case 0x20000510:
        case 0x20100500:
        case 0x20100510:
        case 0x20100520:
        case 0x20200110:
        case 0x20200111:
        case 0x20200130:
        case 0x20200140:
        case 0x20400010:
        case 0x20400020:
        case 0x20400500:
        case 0x20500010:
        case 0x20500500:
        case 0x21000500:
        case 0x21200050:
        case 0x21200070:
        case 0x21300010:
        case 0x21300015:
        case 0x21300030:
        case 0x21300040:
        case 0x21300050:
        case 0x21300060:
        case 0x21300080:
        case 0x213000A0:
        case 0x213000C0:
        case 0x30020030:
        case 0x30040010:
        case 0x30040050:
        case 0x30040060:
        case 0x30060010:
        case 0x30060012:
        case 0x30060014:
        case 0x30060016:
        case 0x30060020:
        case 0x30060030:
        case 0x30060039:
        case 0x30060040:
        case 0x30060080:
        case 0x30060086:
        case 0x300600A0:
        case 0x300600B0:
        case 0x300600C0:
        case 0x30080010:
        case 0x30080020:
        case 0x30080030:
        case 0x30080040:
        case 0x30080050:
        case 0x30080060:
        case 0x30080070:
        case 0x30080080:
        case 0x30080090:
        case 0x300800A0:
        case 0x300800B0:
        case 0x300800C0:
        case 0x300800D0:
        case 0x300800E0:
        case 0x30080100:
        case 0x30080110:
        case 0x30080120:
        case 0x30080130:
        case 0x30080140:
        case 0x30080150:
        case 0x30080160:
        case 0x30080220:
        case 0x300A0010:
        case 0x300A0040:
        case 0x300A0048:
        case 0x300A0070:
        case 0x300A00B0:
        case 0x300A00B6:
        case 0x300A00CA:
        case 0x300A00D1:
        case 0x300A00E3:
        case 0x300A00F4:
        case 0x300A0107:
        case 0x300A0111:
        case 0x300A0116:
        case 0x300A011A:
        case 0x300A0180:
        case 0x300A0190:
        case 0x300A01A0:
        case 0x300A01B4:
        case 0x300A0206:
        case 0x300A0210:
        case 0x300A0230:
        case 0x300A0260:
        case 0x300A0280:
        case 0x300A02B0:
        case 0x300A02D0:
        case 0x300C0002:
        case 0x300C0004:
        case 0x300C000A:
        case 0x300C0020:
        case 0x300C0040:
        case 0x300C0042:
        case 0x300C0050:
        case 0x300C0055:
        case 0x300C0060:
        case 0x300C0080:
        case 0x300C00B0:
        case 0x40080050:
        case 0x40080111:
        case 0x40080117:
        case 0x40080118:
        case 0x54000100:
            return true;
        default:
        //case 0x50XX2600:
			if (tag>=0x50002600 && tag<=0x50ff2600 && (tag&0xffff)==0x2600)
				return true;
            break;
    }
	return false;
}
