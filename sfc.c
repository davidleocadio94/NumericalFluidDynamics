/*
 * ===================================  sfc ===============================
 *   sfc shows 2-d field as a surface and shows the max/min, time & label.
 *
 *   Note: this routine plots the entire field, including any ghost zones.
 *   All array contents must be set (including ghost zones, if any)
 *
 *   Note change in argument list compared to prior versions.
 *   RE-updated to transpose X,Y for correct view when plotted
 *
 *   ATMS 502 / CSE 566, Spring 2017
 *
 *   Arguments:
 * 
 * 	nx,ny	input	integers	dimensions of 's', including ghost zones
 * 	nyuse	input	integer		max value of 2nd dimension actually used
 * 					  ... so nyuse is <= ny
 * 	s	input	real array	field to be displayed.
 * 	simtime	input	real		time of integration
 * 	angh	input	real		horizontal viewing angle counter-
 * 					 clockwise from x-axis
 * 	angv	input	real		vertical viewing angle; >0 = above
 * 					 plane of field average.
 * 	label	input	character*	character label
 *	name	input	character*	your name, to place on plot
 */

#include <ncarg/ncargC.h>
#include <ncarg/gks.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define IWTYPE 1
#define WKID   1
/* sfc(NXDIM,NYDIM,NYDIM,s1[NXDIM][NYDIM],0.0,90.0,30.0,"ic",name);      */
void sfc(int nx,int ny,int nyuse,float s[nx][ny],float simtime,
	 float angh,float angv,char *label,char *name)
{

	int i,j;
	float smin,smax,*workarray;
	float *sTMP,*ptr;
        Gcolr_rep rgb1,rgb2;
	char mmlabel[80],tlabel[40];

	/*
	 * Find min,max
	 */

	smin = s[0][0];
	smax = s[0][0];
	for (j=0; j<nyuse; j++) {
	  for (i=0; i<nx; i++) {
	    if (s[i][j] < smin) smin=s[i][j];
	    if (s[i][j] > smax) smax=s[i][j];
	  }
  	}
	/*printf("Sfc: min %.2f, max %.2f\n",smin,smax);*/

	/*
	 * Create min/max and time labels
	 */

	sprintf(mmlabel,"MIN=%.3f  MAX=%.3f",smin,smax);
	sprintf( tlabel,"TIME = %.4f",simtime);

	/*
	 * Plot labels
	 */

	c_set(0.,1.,0.,1.,0.,1.,0.,1.,1);
	c_pcmequ(0.50,0.97,label  ,0.020,0.0, 0.0);
	c_pcmequ(0.95,0.02,mmlabel,0.015,0.0, 1.0);
	c_pcmequ(0.05,0.02, tlabel,0.015,0.0,-1.0);

	/*
	 * Additional labels
	 */

        c_pcmequ(0.02,0.99,name,0.01,90.,1.);
        c_pcmequ(0.98,0.06,"ATMS 502/CSE 566",0.01,0.,1.);
        c_pcmequ(0.02,0.06,"Spring 2017",0.01,0.,-1.);

	/*
	 * Set colors to default (commented out for now)
	 */

	/*
        rgb1.rgb.red = 1.; rgb1.rgb.green = 1.; rgb1.rgb.blue = 1.;
        rgb2.rgb.red = 0.; rgb2.rgb.green = 0.; rgb2.rgb.blue = 0.;
        gset_colr_rep(WKID,0,&rgb1);
        gset_colr_rep(WKID,1,&rgb2);
	*/

        /*
         * Allocate and fill temporary array of correct dimensions
         * This is needed if 2nd dimension not used in full,
         *   true in program 1 with time counter as 2nd dimension.
         * This is done so routine displays correctly in C
         */

        sTMP = (float *) malloc(nx*nyuse*sizeof(float));
	ptr = sTMP;

        for (j=0; j<nyuse; j++) {
          for (i=0; i<nx; i++) {
            *ptr++ = s[i][j];
          }
        }

	/*
	 * Allocate scratch work array for ezsrfc
	 */

	workarray = (float *) malloc((2*nx*nyuse+nx+nyuse)*sizeof(float));

	/*
	 * Plot surface
	 */

	c_ezsrfc(sTMP,nx,nyuse,angh,angv,workarray);

	free(workarray);
	free(sTMP);
	return;
}
