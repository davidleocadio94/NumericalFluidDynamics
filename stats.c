#include <math.h>
/*
 * ========================= stats =====================
 * Stats computes and prints out the max S values
 * Arguments:
 *
 *	s2	real array	Latest data. Check i1..i2;
 *				  [i1-1],[i2+1] = ghost zones
 *	i1,i2	integers	indices bounding array data
 *	nx	integer		size of data array without ghost zones
 *	n	integer		time step counter
 *	smax	real		holds max absolute value of s2
 */

/*void stats(float s2[],int i1,int i2,int nx,int n,float *smax)*/
void stats(int nx,int ny, int i1, int i2, int j1, int j2,float s[][ny],float  *max,float * min)
{
	int i,j;

	*max = s[i1][j1];
	*min = s[i1][j1];
	for (i=i1; i<=i2; i++) {
		for(j=j1;j<=j2;j++){
	  if (s[i][j] > * max){ * max = s[i][j];}
	  if (s[i][j] < * min){ * min = s[i][j];}
	  }
	}


	return;
}

