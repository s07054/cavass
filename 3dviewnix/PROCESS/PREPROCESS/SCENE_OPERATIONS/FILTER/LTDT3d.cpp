/*
  Copyright 1993-2012, 2015 Medical Image Processing Group
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
	This program performs LTDT (linear time distance transform) on an nD Scene (3DViewnix format).

	
Author: Xinjian Chen
Date  : 03/05/2009

*/


#include <math.h>

#include <Viewnix.h>
 
#include  "cv3dv.h"

#include "slices.c"
#include "fff.c"

extern "C" 
{
 int VGetHeaderLength ( FILE* fp, int* hdrlen );
 int VAddBackgroundProcessInformation ( char* command );
}

#define INTFINITY  300000000
#define MAXDIM    10
#define MAXDIMLENGTH  10000
#define MAXDIST  INTFINITY

#define LTSDT_DIST

int nDim;

int  GetCoordinates( int index, int nCoord[], int *pnDimEach );
int GetCoordIndex( int nCoord[], int *pnDimEach );
int Check( int u, int v, int w, int d, int nCoord[], int *pnDimEach );

double GetCenter_d( int nCoordu[], int nCoordv[], int nCoord[], int d)
{
	double dCenterd;
	double dSum = 0;
	int i;

	for( i=0; i<nDim; i++ )
	{
		if( i== d )
			continue;

		dSum += (nCoordv[i] - nCoord[i] )*(nCoordv[i] - nCoord[i] ) - (nCoordu[i] - nCoord[i] )*(nCoordu[i] - nCoord[i] ); // 
	}

	if( nCoordu[d] == nCoordv[d] )
		dCenterd = nCoordv[d];
	else
		dCenterd = ( dSum + nCoordv[d]*nCoordv[d] - nCoordu[d]*nCoordu[d] ) / (2*( nCoordv[d] - nCoordu[d] ) );

	return dCenterd;
}

int Check( int u, int v, int w, int d, int nCoord[], int *pnDimEach )
{
	int nCoordu[MAXDIM]; 	
	int nCoordv[MAXDIM]; 
	int nCoordw[MAXDIM]; 
	/*double x_uvd;
	double x_vwd;*/

	GetCoordinates( u, nCoordu, pnDimEach );
	GetCoordinates( v, nCoordv, pnDimEach );	
	GetCoordinates( w, nCoordw, pnDimEach );

	int a = nCoordv[d] - nCoordu[d];
	int b = nCoordw[d] - nCoordv[d];
	int c = a+b;

	int vD2 = 0;
	for( int i=0; i<nDim; i++ )
	{
		if( i== d )
			continue;

		vD2 += (nCoordv[i] - nCoord[i] )*(nCoordv[i] - nCoord[i] ); // 
	}

	int uD2 = 0;
	for( int i=0; i<nDim; i++ )
	{
		if( i== d )
			continue;

		uD2 += (nCoordu[i] - nCoord[i] )*(nCoordu[i] - nCoord[i] ); // 
	}

	int wD2 = 0;
	for( int i=0; i<nDim; i++ )
	{
		if( i== d )
			continue;

		wD2 += (nCoordw[i] - nCoord[i] )*(nCoordw[i] - nCoord[i] ); // 
	}

	return ( c*vD2 - b*uD2 - a*wD2 - a*b*c > 0 );

//	x_uvd = (int)( GetCenter_d(nCoordu, nCoordv, nCoord, d) +0.999999999 );  // LTDT
//	x_vwd = (int)( GetCenter_d(nCoordv, nCoordw, nCoord, d) );               // LTDT

	/*x_uvd = GetCenter_d(nCoordu, nCoordv, nCoord, d);
	x_vwd = GetCenter_d(nCoordv, nCoordw, nCoord, d);*/


	/*if( x_uvd > x_vwd )
		return 1;
	else
		return 0;*/

}

int GetCoordIndex( int nCoord[], int *pnDimEach )
{
	int index = 0;
	int nMul = 1;
	int i;
	for( i=0; i<nDim; i++ )
	{
		index += nCoord[i] * nMul;
		nMul *= pnDimEach[i];
	}

	return index;
}


int Trim(  int c, int d, int F[], int q[], int *pnDimEach, int nVoxelsNum  )
{
	int m = -1;
	int xi;
	int i;

	int nCoord[MAXDIM]; // suppose max dimension 
	GetCoordinates( c, nCoord, pnDimEach );

	for( i=0; i< pnDimEach[d]; i++ )
	{
		nCoord[d] = i;
		xi = GetCoordIndex( nCoord, pnDimEach );

		if( F[xi] != -1 )  // -1 is empty (NULL) set
		{
			if( m < 1 )
			{
				m++;
			    q[m] = F[xi];
			}
			else
			{
				while ( m >= 1 && Check( q[m-1], q[m], F[xi], d, nCoord, pnDimEach ) )   // start from 0, so +1
				{					
					m--;
				}

				m++;
				q[m] = F[xi];

			}
			//m = m+1;
			//q[m] = F[xi];

			//if ( m+1 > 1)  // start from 0, so +1
			//{
			//	u = q[m-1];
			//	v = q[m];
			//	
			//	GetCoordinates( u, nCoordu, pnDimEach );				
			//	GetCoordinates( v, nCoordv, pnDimEach );

			//	x_uvd = GetCenter_d(nCoordu, nCoordv, nCoord, d);
			//	if( x_uvd >= pnDimEach[d]-1 )
			//		m = m-1;
			//	else
			//	{
			//		while ( m+1>2 && Check( q[m-2], q[m-1], q[m], d, nCoord, pnDimEach ) )   // start from 0, so +1
			//		{
			//			q[m-1] = q[m];
			//			m = m-1;
			//		}

			//		if ( m+1 == 2 )
			//		{
			//			u = q[0];
			//			v = q[1];
			//			
			//			GetCoordinates( u, nCoordu, pnDimEach );						
			//			GetCoordinates( v, nCoordv, pnDimEach );

			//			x_uvd = GetCenter_d(nCoordu, nCoordv, nCoord, d);

			//			if( x_uvd < 0 )
			//			{
			//				q[0] = q[1];
			//				m = 0;
			//			}
			//		}
			//	}
			//} //if ( m+1 > 1)
		}
	}

	return m+1;  //m+1;
}

int  GetCoordinates( int index, int nCoord[], int *pnDimEach )
{
	int dim;
	int posi = index;	
	for( dim = 0; dim < nDim; dim++ )
	{
		nCoord[dim] = posi % pnDimEach[dim];
		posi /= pnDimEach[dim];
	}

	return 1;
}


double dist( int i, int j, int *pnDimEach )
{
	// get coordinates of i, j
	int nCoordi[MAXDIM]; // suppose max dimension 
	int nCoordj[MAXDIM];
	double dist2 = 0; 
	int dim;

	GetCoordinates( i, nCoordi, pnDimEach );
	GetCoordinates( j, nCoordj, pnDimEach );

	dist2 = 0;        // here, suppose p = 2
	for( dim = 0; dim < nDim; dim++ )
	{
		dist2 += ( nCoordj[dim] - nCoordi[dim] ) * ( nCoordj[dim] - nCoordi[dim] );
	}

//	dist2 = sqrt( dist2 );

	return dist2;
}

int DimUp( int c, int d, int F[], int *pnDimEach, int nVoxelsNum )
{
	int q[MAXDIMLENGTH];
	int m;
	int l = 0;
	int nCoord[MAXDIM]; // suppose max dimension 
	int i;
	int xi;
	
	m = -1;	
	
	GetCoordinates( c, nCoord, pnDimEach );

	for( i=0; i< pnDimEach[d]; i++ )
	{
		nCoord[d] = i;
		xi = GetCoordIndex( nCoord, pnDimEach );

		if( F[xi] != -1 )  // -1 is empty (NULL) set
		{
			if( m < 1 )
			{
				m++;
			    q[m] = F[xi];
			}
			else
			{
				while ( m >= 1 && Check( q[m-1], q[m], F[xi], d, nCoord, pnDimEach ) )   // start from 0, so +1  // take 2 more secs
				{					
					m--;
				}

				m++;
				q[m] = F[xi];

			}
		
		}
	}
//	m = Trim( c, d, F, q, pnDimEach, nVoxelsNum  );

	if( m >= 0 )
	{		
		GetCoordinates( c, nCoord, pnDimEach );

		//k=0;
		l = 0;

		for( i=0; i<pnDimEach[d]; i++ )
		{
			nCoord[d] = i;
			xi = GetCoordIndex( nCoord, pnDimEach );

			while( l < m && dist( xi, q[l], pnDimEach ) > dist( xi, q[l+1], pnDimEach ) )   // take 4 more secs
				l++;

			F[xi] = q[l];
		}

		//l = 0;
		//while( k> 0 )  // Queue is not empty
		//{
		//	x = Queue[0];

		//	for( i=0; i< k-1; i++ )
		//		Queue[i] = Queue[i+1];

		//	k--;

		//	while( l < m-1 && dist( x, q[l], pnDimEach ) >= dist( x, q[l+1], pnDimEach ) )
		//		l++;

		//	F[x] = q[l];
		//}

	}

	return 1;
}

int LTDT(unsigned char pInImg[],  unsigned char pBinaryImg[], int F[], float DT[], int *pnDimEach, int nVoxelsNum )
{
	int i, d;
	int x;
	int nCoordi[MAXDIM]; // suppose max dimension 


	//for( i= 0; i< nVoxelsNum; i++ )  //// Bg --> Fg distance
	//{
	//	if( pBinaryImg[i] > 0 )   
	//		F[i] = i;
	//	else
	//		F[i] = -1; // -1 represent empty set
	//}

	for( i= 0; i< nVoxelsNum; i++ )  //// Fg --> Bg distance
	{
		if( pBinaryImg[i] > 0 )   
			F[i] = -1;
		else
			F[i] = i; // -1 represent empty set
	}

	///////// Design for 3D Feature Transform //////////
	int DimTable1[3];
	int DimTable2[3];
	DimTable1[0] = pnDimEach[1];	DimTable2[0] = pnDimEach[2];  // x == 0
	DimTable1[1] = pnDimEach[0];	DimTable2[1] = pnDimEach[2];  // y == 0
	DimTable1[2] = pnDimEach[0];	DimTable2[2] = pnDimEach[1];  // z == 0

	int index1[3] = {1, 0, 0};
	int index2[3] = {2, 2, 1};
	
	for( d=0; d<nDim; d++ )
	{		
		//for( i= 0; i< nVoxelsNum; i++ )
		//{			
		//	GetCoordinates( i, nCoordi, pnDimEach );

		//	if( nCoordi[d] == 0 )  // passing (0,0,0) && x(d) = 0
		//	{				
		//		x = i;

		//		DimUp( x, d, F, pnDimEach,	nVoxelsNum);
		//	}
		//}		
		nCoordi[d] = 0;
		for( int i = 0; i< DimTable1[d]; i++)
		{
			nCoordi[ index1[d] ] = i;
			for( int j = 0; j< DimTable2[d]; j++)
			{
				nCoordi[ index2[d] ] = j;

				x = GetCoordIndex( nCoordi, pnDimEach );

				DimUp( x, d, F, pnDimEach,	nVoxelsNum);
			}
		}
		
	}	
	
	return 1;

}



bool ITKCheck( float d1, float d2, float d3, int ud, int vd, int wd )
{	
	int a = vd - ud; 
	int b = wd - vd;
	int c = wd - ud;
	
	return ( c*d2 - b*d1 - a*d3 - a*b*c > 0 );
}


int ITKDimUp( int c, int d, float DT[], int *pnDimEach, int nVoxelsNum )
{
	float q[MAXDIMLENGTH];
	int   indexD[MAXDIMLENGTH];

	int m = -1;
	int l = 0;
	int nCoord[MAXDIM]; // suppose max dimension 
	int i;
	int xi;	
	
	GetCoordinates( c, nCoord, pnDimEach );

	for( i=0; i< pnDimEach[d]; i++ )
	{
		nCoord[d] = i;
		xi = GetCoordIndex( nCoord, pnDimEach );

		if( DT[xi] != INTFINITY )  // -1 is empty (NULL) set
		{
			if( m < 1 )
			{
				m++;
			    q[m] = DT[xi];	indexD[m] = i;
			}
			else
			{
				while ( m >= 1 && ITKCheck( q[m-1], q[m], DT[xi], indexD[m-1], indexD[m], i ) )   // start from 0, so +1
				{					
					m--;
				}

				m++;
				q[m] = DT[xi];  indexD[m] = i;

			}
		
		}
	}	
	
	if( m >= 0 )
	{	
		l = 0;

		for( i=0; i<pnDimEach[d]; i++ )
		{
			nCoord[d] = i;
			xi = GetCoordIndex( nCoord, pnDimEach );

			while( l < m && ( q[l] + (indexD[l] -i) * (indexD[l] -i) > q[l+1] + (indexD[l+1] -i) * (indexD[l+1] -i) ) )
				l++;
			
			DT[xi] = q[l] + (indexD[l] -i) * (indexD[l] -i);
			
		}
	}

	return 1;
}

int ITKLTDT(unsigned char pInImg[],  unsigned char pBinaryImg[], float DT[], int *pnDimEach, int nVoxelsNum )
{
	int i, d;
	int x;
	int nCoordi[MAXDIM]; // suppose max dimension 	

	for( i= 0; i< nVoxelsNum; i++ )  //// Fg --> Bg distance
	{
		if( pBinaryImg[i] > 0 )   
			DT[i] = INTFINITY;
		else
			DT[i] = 0; // -1 represent empty set
	}

	///////// Design for 3D Feature Transform //////////
	int DimTable1[3];
	int DimTable2[3];
	DimTable1[0] = pnDimEach[1];	DimTable2[0] = pnDimEach[2];  // x == 0
	DimTable1[1] = pnDimEach[0];	DimTable2[1] = pnDimEach[2];  // y == 0
	DimTable1[2] = pnDimEach[0];	DimTable2[2] = pnDimEach[1];  // z == 0

	int index1[3] = {1, 0, 0};
	int index2[3] = {2, 2, 1};
	
	for( d=0; d<nDim; d++ )
	{	
		nCoordi[d] = 0;
		for( int i = 0; i< DimTable1[d]; i++)
		{
			nCoordi[ index1[d] ] = i;
			for( int j = 0; j< DimTable2[d]; j++)
			{
				nCoordi[ index2[d] ] = j;

				x = GetCoordIndex( nCoordi, pnDimEach );

				ITKDimUp( x, d, DT, pnDimEach,	nVoxelsNum);
			}
		}
		
	}	

	for( i= 0; i< nVoxelsNum; i++ )  //// Fg --> Bg distance
	{
		DT[i] = (float)sqrt( double(DT[i]) );
	}
	
	return 1;

}



int FTDTDimUp( int c, int d, int F[], float DT[], int *pnDimEach, int nVoxelsNum )
{
	float q[MAXDIMLENGTH];
	int   indexD[MAXDIMLENGTH];
	int   g[MAXDIMLENGTH];

	int m = -1;
	int l = 0;
	int nCoord[MAXDIM]; // suppose max dimension 
	int i;
	int xi;	
	
	GetCoordinates( c, nCoord, pnDimEach );

	for( i=0; i< pnDimEach[d]; i++ )
	{
		nCoord[d] = i;
		xi = GetCoordIndex( nCoord, pnDimEach );

		if( DT[xi] != INTFINITY )  // -1 is empty (NULL) set
		{
			if( m < 1 )
			{
				m++;
			    q[m] = DT[xi];	indexD[m] = i;	g[m] = F[xi];  //g[m]  saved the feature transform
			}
			else
			{
				while ( m >= 1 && ITKCheck( q[m-1], q[m], DT[xi], indexD[m-1], indexD[m], i ) )   // start from 0, so +1
				{					
					m--;
				}

				m++;
				q[m] = DT[xi];  indexD[m] = i;	g[m] = F[xi];

			}
		
		}
	}	
	
	if( m >= 0 )
	{	
		l = 0;

		for( i=0; i<pnDimEach[d]; i++ )
		{
			nCoord[d] = i;
			xi = GetCoordIndex( nCoord, pnDimEach );

			while( l < m && ( q[l] + (indexD[l] -i) * (indexD[l] -i) > q[l+1] + (indexD[l+1] -i) * (indexD[l+1] -i) ) )
				l++;
		
			DT[xi] = q[l] + (indexD[l] -i) * (indexD[l] -i);
			F[xi]  = g[l];

		}
	}

	return 1;
}

int LTFTDT(unsigned char pInImg[],  unsigned char pBinaryImg[], int F[], float DT[], int *pnDimEach, int nVoxelsNum )
{
	int i, d;
	int x;
	int nCoordi[MAXDIM]; // suppose max dimension 	

	for( i= 0; i< nVoxelsNum; i++ )  //// Fg --> Bg distance
	{
		if( pBinaryImg[i] > 0 )   
		{
			DT[i] = INTFINITY;
			F[i]  = -1;
		}
		else
		{
			DT[i] = 0; // -1 represent empty set
			F[i]  = i;
		}
	}

	///////// Design for 3D Feature Transform //////////
	int DimTable1[3];
	int DimTable2[3];
	DimTable1[0] = pnDimEach[1];	DimTable2[0] = pnDimEach[2];  // x == 0
	DimTable1[1] = pnDimEach[0];	DimTable2[1] = pnDimEach[2];  // y == 0
	DimTable1[2] = pnDimEach[0];	DimTable2[2] = pnDimEach[1];  // z == 0

	int index1[3] = {1, 0, 0};
	int index2[3] = {2, 2, 1};
	
	for( d=0; d<nDim; d++ )
	{	
		nCoordi[d] = 0;
		for( int i = 0; i< DimTable1[d]; i++)
		{
			nCoordi[ index1[d] ] = i;
			for( int j = 0; j< DimTable2[d]; j++)
			{
				nCoordi[ index2[d] ] = j;

				x = GetCoordIndex( nCoordi, pnDimEach );

				FTDTDimUp( x, d, F, DT, pnDimEach,	nVoxelsNum);
			}
		}
		
	}	

	for( i= 0; i< nVoxelsNum; i++ )  //// Fg --> Bg distance
	{
		DT[i] = (float)sqrt( double(DT[i]) );
	}

	return 1;

}


int  GetDT( int F[], float DT[], int nLength, int *pnDimEach )
{
	int i;
	for( i= 0; i< nLength; i++ )
	{
		if( F[i] == -1 ) // -1 represent empty set
			DT[i] = INTFINITY;
		else
			DT[i] = (float)dist( i, F[i], pnDimEach ); 
	}

	return 1;
}


int LTSDT(  unsigned char pBinaryImg[], int pF[], float pDT[], int *pnDimEach, int nVoxelsNum, int distType )
{
	int i, j, k;
	int x, y;
	int nCoordx[MAXDIM]; // suppose max dimension 
	int nCoordy[MAXDIM]; // suppose max dimension 	
//	int *pF;
	unsigned char * pBoundaryImg;
	

	time_t ltime;
    time( &ltime );	
	
	pBoundaryImg = (unsigned char *)malloc( nVoxelsNum * sizeof( unsigned char) );
	if( pBoundaryImg == NULL )
		return 0;

	memset(pBoundaryImg, 1, nVoxelsNum  * sizeof( unsigned char) );

	int boundaryValue = 0;
	if( distType == 0 )  // background to goreground	
		boundaryValue = 255;
	else
		boundaryValue = 0;

	for( x= 0; x< nVoxelsNum; x++ )
	{
		if( pBinaryImg[x] == boundaryValue )
		{
			GetCoordinates( x, nCoordx, pnDimEach );			

			for( i= nCoordx[0] -1; i<= nCoordx[0] +1; i++ )   // here, only support 3D, 18 neighbour adjacency
			{
				if( i< 0 || i> pnDimEach[0] - 1 )
					continue;

				nCoordy[0] = i;

				for( j= nCoordx[1] -1; j<= nCoordx[1] +1; j++ )
				{
					if( j< 0 || j> pnDimEach[1] - 1 )
						continue;

					nCoordy[1] = j;

					for( k= nCoordx[2] -1; k<= nCoordx[2] +1; k++ )
					{
						if( k< 0 || k> pnDimEach[2] - 1 )
							continue;

						nCoordy[2] = k;

						if( nCoordy[0] != nCoordx[0] && nCoordy[1] != nCoordx[1] && nCoordy[2] != nCoordx[2] )  // get rid of 8 corners
							continue;

						y = GetCoordIndex( nCoordy, pnDimEach );

						if( pBinaryImg[y] == 255-boundaryValue )
						{
							pBoundaryImg[x] = 0;
							break;
						}
					}

					if( pBoundaryImg[x] == 0 )						
						break;						
				}

				if( pBoundaryImg[x] == 0 )						
						break;
			}
		}

	}

	//time_t ltime2;
 //   time( &ltime2 );

	//int secs = ltime2-ltime;

	///// output result
	//FILE *fp = NULL;
	//fp = fopen("C:/xinjian/Program/LinearDistanceTransform/LTSDT_WithFTDT.txt", "a+");
	//if( fp == NULL )
	//return 0;

	//fprintf(fp, "initial time = %d \n", secs );		

	//pF = (int *)malloc( nVoxelsNum * sizeof( int) );
	//if( pF == NULL )
	//	return 0;

//	LTDT(  NULL,  pBoundaryImg, pF, pDT, pnDimEach,	nVoxelsNum );  // 1.1 old version, slow

	//time_t ltime3;
 //   time( &ltime3 );

	//secs = ltime3-ltime2;
	//fprintf(fp, "FT time = %d \n", secs );		

	
	if( pDT == NULL )
	{
		pDT = (float *)malloc( (nVoxelsNum) * sizeof( float) );  // nor enough memory for super dataset
		if( pDT == NULL )
			return 0;
	}
	
//	GetDT( pF, pDT,nVoxelsNum, pnDimEach );   // 1.2  Get DT

//	ITKLTDT(  NULL,  pBoundaryImg, pDT, pnDimEach,	nVoxelsNum );  // 2. new version, according to PAMI Maurer paper

	LTFTDT(  NULL,  pBoundaryImg, pF, pDT, pnDimEach,	nVoxelsNum );  // 3. Based on 2, including Feature transform	
	
	for( x= 0; x< nVoxelsNum; x++ )
	{
		if( distType == 0 ) //Distance Type, 0: background --> foreground; 1: foreground --> background; 2: both;\n");
		{
			if( pBinaryImg[x] != 0 )				
				pDT[x] = 0; 
		}
		else if( distType == 1 ) 
		{
			if( pBinaryImg[x] == 0 )
				pDT[x] = 0;			
		}
		else
		{
			if( pBinaryImg[x] == 0 )  // background to foreground, positive
				pDT[x] = pDT[x];
			else
				pDT[x] = -pDT[x];     // foreground to background, negative
		}

	}
	
	//if( pF != NULL )
	//	free(pF);
	
	//time_t ltime4;
 //   time( &ltime4 );

	//secs = ltime4-ltime3;

	///// output result	

	//fprintf(fp, "DT time = %d, Total time = %d \n", secs, ltime4 - ltime );		

	//fclose(fp);

	return 1;

}

// GoldenStandard test

int GoldenLTDT(unsigned char pInImg[],  unsigned char pBinaryImg[], int F[], float DT[], int *pnDimEach, int nVoxelsNum )
{

	int i, j;
	double distij;
	double minDist;
	int minj;
	

	for( i= 0; i< nVoxelsNum; i++ )
	{
		if( pBinaryImg[i] == 0 )
			F[i] = i;
		else
			F[i] = -1; // -1 represent empty set
	}
	
	for( i= 0; i< nVoxelsNum; i++ )
	{
		if( F[i] == i )  // background voxel
			continue;
		
		minDist = MAXDIST;
		minj = -1;
		
		for( j= 0; j< nVoxelsNum; j++ )
		{			
			if( F[j] == j )  // background voxel
			{
				distij = dist( i, j, pnDimEach  );
				if( distij < minDist )
				{
					minDist = distij;
					minj = j;
				}
			}
		}	

		F[i] = minj;

	}

	return 1;

}

int GBDT(  unsigned char pBinaryImg[], float pDT[], int *pnDimEach, int nVoxelsNum )
{
	int i, d;
	int x, y, c;
	int nVoxelsNumX2 = nVoxelsNum; 
	int nCoordx[MAXDIM]; // suppose max dimension 
	int nCoordy[MAXDIM]; // suppose max dimension 
	int nCoordc[MAXDIM]; // suppose max dimension 
	int pnDimEachX2[MAXDIM];

	time_t ltime;
    time( &ltime );	

	unsigned char* pBinaryImgX2 = NULL;
	for( int i=0; i<nDim; i++ )
	{
		nVoxelsNumX2   *= 2;
		pnDimEachX2[i] = 2*pnDimEach[i];
	}

	pBinaryImgX2 = (unsigned char *)malloc( nVoxelsNumX2 * sizeof( unsigned char) );
	if( pBinaryImgX2 == NULL )
		return 0;

	memset(pBinaryImgX2, 1, nVoxelsNumX2  * sizeof( unsigned char) );

	for( x= 0; x< nVoxelsNum; x++ )
	{
		for( d=0; d<nDim; d++ )
		{
			GetCoordinates( x, nCoordx, pnDimEach );

			for( i = 0; i< nDim; i++ )
			{				
				if( i != d )
					nCoordy[i] = nCoordx[i];
				else
				{
					if( nCoordx[i] < pnDimEach[i] - 1 )  // boundary point, need check ???
						nCoordy[i] = nCoordx[i] +1;
					else
						nCoordy[i] = nCoordx[i];

				}
			}

			y = GetCoordIndex( nCoordy, pnDimEach );

			if( pBinaryImg[x] != pBinaryImg[y] )
			{
				for( i = 0; i< nDim; i++ )
				{				
					if( i != d )
						nCoordc[i] = 2*nCoordx[i];
					else
						nCoordc[i] = 2*nCoordx[i] +1;
				}

				c = GetCoordIndex( nCoordc, pnDimEachX2 );

				pBinaryImgX2[c] = 0;

			}
		}		
	}

	time_t ltime2;
    time( &ltime2 );

	int secs = (int)(ltime2-ltime);

	/// output result
	//FILE *fp = NULL;
	//fp = fopen("C:/xinjian/Program/LinearDistanceTransform/DistTransformResult05.08.txt", "a+");
	//if( fp == NULL )
	//return 0;

	//fprintf(fp, "initial time = %d \n", secs );		
	

	int *pFX2  = NULL;  // Feature Transform Image,   Here is the closest(nearest)  point coordinates
	float *pDTX2  = NULL;  // Feature Transform Image

	pFX2 = (int *)malloc( nVoxelsNumX2 * sizeof( int) );
	if( pFX2 == NULL )
		return 0;

	/////////LTDT(  NULL,  pBinaryImgX2, pFX2, pDTX2, pnDimEachX2, nVoxelsNumX2 ); // old version	
	

	pDTX2 = (float *)malloc( nVoxelsNumX2 * sizeof( float) );
	if( pDTX2 == NULL )
		return 0;

    LTFTDT(  NULL, pBinaryImgX2, pFX2, pDTX2, pnDimEachX2, nVoxelsNumX2 );
	//ITKLTDT(  NULL,  pBinaryImgX2, pDTX2, pnDimEachX2, nVoxelsNumX2 );
//	GoldenLTDT(  NULL,  pBinaryImgX2, pFX2, pDTX2, pnDimEachX2, nVoxelsNumX2 );

	//GetDT( pFX2, pDTX2,nVoxelsNumX2, pnDimEachX2 );

	for( x= 0; x< nVoxelsNum; x++ )
	{
		GetCoordinates( x, nCoordx, pnDimEach );
		for( i=0; i<nDim; i++)
			nCoordx[i] *= 2;

		y = GetCoordIndex( nCoordx, pnDimEachX2 );


		if( pBinaryImg[x] == 0 )
			pDT[x] = pDTX2[y]/2;
		else
			pDT[x] = -pDTX2[y]/2;

	}


	if( pBinaryImgX2 != NULL )
		free(pBinaryImgX2);
	if( pFX2 != NULL )
		free(pFX2);
	if( pDTX2 != NULL )
		free(pDTX2);

	time_t ltime3;
    time( &ltime3 );

	secs = (int)(ltime3-ltime2);

	/// output result	

	//fprintf(fp, "GBDT time = %d \n", secs );		

	//fclose(fp);

	return 1;

}

/*-------------------------------------------------------------------------*/
/* Read the file header and its length */
int get_file_info(char *file, ViewnixHeader *vh, int *len)
//ViewnixHeader *vh;
//char *file; /* filename */
//int *len;   /* length of header */
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
 


int main(int argc,char* argv[])
{
	FILE *fpin, *fpout;	/* inpput/output files */
	int execution_mode;	/* execution mode */
	ViewnixHeader vh;	/* 3DViewnix header */
	SLICES	sl;			/* Structure containing information about the slices of the scene */
	char group[5],		/* Used in VWriteHeader */
		element[5];
	int length;			/* length of a slice */
	int hlength;		/* length of the input header */
	int width, height;	/* dimensions of a slice */
	int distType;		/* Distance Type, 0: background --> foreground; 1: foreground --> background; 2: both */
	int i,j,k;			/* general use */
	int error;			/* error code */
	char *comments;     /* used to modify the header (description field) */
	int binary_flag=0;	/* indicates if scene is binary or not */

	unsigned char *in_buffer1b;
	unsigned char *in_buffer8b,	*out_buffer8;
	int *pF  = NULL;  // Feature Transform Image,   Here is the closest(nearest)  point coordinates
	float *pDT  = NULL;  // Distance Transform Image

		
	int pnDimEach[MAXDIM];
	int nVoxelsNum;
	unsigned char *pBinaryImg = NULL;
	int ii, jj;
	int tmpLen;

	float min = 65535;
	float max = -65535;
	float smallest, largest;
	float scaleLarge;
	float var;
	bool bFTSave = false;


	if (argc < 6)
	{
		printf("Usage:\n");
		printf("%% LTDT3D input output mode distType bFTSave\n");
		printf("where:\n");
		printf("input    : name of input file;\n");
		printf("output   : name of output file;\n");
		printf("mode     : mode of operation (0=foreground, 1=background);\n");
		printf("distType : Distance Type, 0: background --> foreground; 1: foreground --> background; 2: both; 3: Digital Boundary\n");
		printf("bFTSave  : 1: Save FT Image 0; Don't save FT Image. \n");
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
 
    /* Get EXECUTION MODE */
    sscanf(argv[3], "%d", &execution_mode);

	/* Get distType */
	sscanf(argv[4], "%d", &distType);

	/* Get bFTSave */
	sscanf(argv[5], "%d", (int *)&bFTSave );
 

	/* If in background mode, then place an entry in the BG_PROCESS.COM file */
    if(execution_mode == 1)
        VAddBackgroundProcessInformation((char *)"ndLTDT");


 
    /*-----------------------*/
    /* Read 3DViewnix header */
    /*-----------------------*/
    get_file_info(argv[1], &vh, &hlength);
 
	if( vh.scn.num_of_bits != 1 )
	{
		printf("ERROR: The input image need to be a binary image.\n");
		exit(1);
	}

	/* Comoute information about the slices of the scene (number, locations, etc...) */
	compute_slices(&vh, &sl);


	/* Calculate length of a slice */
	width =  vh.scn.xysize[0];
	height =  vh.scn.xysize[1];

	length = (width * height + 7)  / 8;


	/* Allocate memory */

	/* create buffer for one binary image */
	if( (in_buffer1b = (unsigned char *) calloc(1, length) ) == NULL)
	{
		printf("ERROR: Can't allocate input image buffer.\n");
		exit(1);
	}

	length = length*8;
	binary_flag = 1;

	/* create buffer for one grey image */
	if( (in_buffer8b = (unsigned char *) calloc(1, length) ) == NULL)
	{
		printf("ERROR: Can't allocate input image buffer.\n");
		exit(1);
	}

	/* create buffer for one 16 bits grey-level image */
	if( (out_buffer8 = (unsigned char *) calloc(1, length ) ) == NULL)
	{
		printf("ERROR: Can't allocate output image buffer.\n");
		exit(1);
	}

	

	nDim = vh.scn.dimension;  // input dimension num	

	pnDimEach[0] = vh.scn.xysize[0]; //256;  x
	pnDimEach[1] = vh.scn.xysize[1]; //256;  y
	pnDimEach[2] = vh.scn.num_of_subscenes[0];

	nVoxelsNum = 1;
	for( i=0; i< nDim; i++ )  		nVoxelsNum *= pnDimEach[i];

	pBinaryImg = (unsigned char *)malloc( nVoxelsNum * sizeof( unsigned char) );
	if( pBinaryImg == NULL )
		return 0;

	/* Traverse ALL VOLUMES */	
	k = 0;
	for(j=0; j<sl.volumes; j++)
	{
		fseek(fpin, (k*length/8)+hlength, 0);
		/* For each Volume, traverse ALL SLICES */
		for(i=0; i<sl.slices[j]; i++)
		{
			tmpLen = fread(in_buffer1b, 1, (length/8), fpin);
			if( tmpLen != length/8)
			{
				printf("ERROR: Couldn't read slice #%d of volume #%d.\n", i+1, j+1);
				exit(2);
			}
			
			/* Convert to 8 Bits */
			bin_to_grey(in_buffer1b, (length/8), in_buffer8b, 0, 255);

			for(  ii=0; ii< pnDimEach[0]; ii++ )
				for( jj=0; jj< pnDimEach[1]; jj++ )
					pBinaryImg[ i * pnDimEach[1]*pnDimEach[0] + jj*pnDimEach[0] + ii ] = in_buffer8b[jj*pnDimEach[0] + ii];

		}
	}



#ifdef LTSDT_DIST
	pF = (int *)malloc( nVoxelsNum * sizeof( int) );
	if( pF == NULL )
		return 0;

	pDT = (float *)malloc( (nVoxelsNum) * sizeof( float) );  // nor enough memory for super dataset
	if( pDT == NULL )
		return 0;

	//LTSDT(  pBinaryImg, pF, pDT, pnDimEach,	nVoxelsNum, distType );
	////GBDT(  pBinaryImg, pDT, pnDimEach,	nVoxelsNum );
	if( distType == 3 )
		GBDT(  pBinaryImg, pDT, pnDimEach,	nVoxelsNum );
	else
		LTSDT(  pBinaryImg, pF, pDT, pnDimEach,	nVoxelsNum, distType );	

#else

	pF = (int *)malloc( nVoxelsNum * sizeof( int) );
	if( pF == NULL )
		return 0;


	LTDT(  NULL,  pBinaryImg, pF, pDT, pnDimEach,	nVoxelsNum );

	if( pBinaryImg != NULL )
		free(pBinaryImg);

	pDT = (float *)malloc( (nVoxelsNum) * sizeof( float) );  // nor enough memory for super dataset
	if( pDT == NULL )
		return 0;

	GetDT( pF, pDT,nVoxelsNum, pnDimEach );
#endif

	float fResolutionF = (float)0.4;
	for(j=0; j<sl.volumes; j++)
	{
		for(i=0; i<sl.slices[j]; i++)
		{
			for(  ii=0; ii< pnDimEach[0]; ii++ )
				for( jj=0; jj< pnDimEach[1]; jj++ )
				{
					float value = pDT[ i * pnDimEach[1]*pnDimEach[0] + jj*pnDimEach[0] + ii ]/fResolutionF;
					if(  value < min )
						min = value;
					
					if (  value > max )
						max = value;
				}

		}
	}

	/*-------------------------*/
	/* Modify 3DViewnix Header */
	/* In case of binary, change header to 8 bits */
	smallest = 0;
    largest = max;
    if (-min > largest) largest = -min;
	scaleLarge = largest;
    if (largest > 255)  largest = 255;
    vh.scn.num_of_density_values  = 1;
    vh.scn.smallest_density_value = &smallest;
    vh.scn.largest_density_value  = &largest;
    vh.scn.num_of_bits = 8;
	vh.scn.bit_fields[1] = 7;

	/* Get the filenames right (own and parent) */
    strcpy(vh.gen.filename1, argv[1]);
    strcpy(vh.gen.filename, argv[2]);
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

	for(j=0; j<sl.volumes; j++)
	{
		for(i=0; i<sl.slices[j]; i++)
		{
			for(  ii=0; ii< pnDimEach[0]; ii++ )
				for( jj=0; jj< pnDimEach[1]; jj++ )
				{
					var = pDT[ i * pnDimEach[1]*pnDimEach[0] + jj*pnDimEach[0] + ii ]/fResolutionF;
					if ( var < 0 )
						var = -var;

					if ( scaleLarge > 255 )					
						((unsigned char*) out_buffer8)[jj*pnDimEach[0] + ii] = (unsigned char)((var/scaleLarge)*255.0);					
					else
						((unsigned char*) out_buffer8)[jj*pnDimEach[0] + ii] = (unsigned char)var;
				}

			if (VWriteData( (char *)out_buffer8, 1, pnDimEach[1]*pnDimEach[0],fpout,&k))
				printf("Could not write data\n");

		}
	}

	VCloseData(fpout);

	// The following is: Output the feature transform image as IM0, the value is the coordinate index in the image
	if( bFTSave )
	{
		char FT_FileName[200];
		strcpy( FT_FileName, argv[2] );
		strcat( FT_FileName, "_FT.IM0");
		if( (fpout = fopen( FT_FileName, "wb")) == NULL)
		{
			printf("ERROR: Can't open OUTPUT file !\n");
			exit(1);
		}

		smallest = 0;
		largest = (float)nVoxelsNum;
		vh.scn.num_of_density_values  = 1;
		vh.scn.smallest_density_value = &smallest;
		vh.scn.largest_density_value  = &largest;
		vh.scn.num_of_bits = 32;
		vh.scn.bit_fields[1] = 31;

		/* Write output 3DViewnix Header */
		error = VWriteHeader(fpout, &vh, group, element);
		if(error < 100)
		{
			printf("ERROR: Can't write output Header (ERROR #%d, %s-%s)!\n", error,group, element);
			exit(error);
		}	

		if( VWriteData( (char *)pF, 4, nVoxelsNum,fpout,&k) )
			printf("Could not write data\n");

		VCloseData(fpout);

	}

	if( pF != NULL )
	{
		free(pF);
		pF = NULL;
	}
	if( pDT != NULL )
	{
		free(pDT);
		pDT = NULL;
	}

	if( out_buffer8 != NULL )
	{
		free(out_buffer8);
		out_buffer8 = NULL;
	}
	if( in_buffer8b != NULL )
	{
		free(in_buffer8b);
		in_buffer8b = NULL;
	}
	if( in_buffer1b != NULL )
	{
		free(in_buffer1b);
		in_buffer1b = NULL;
	}

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
