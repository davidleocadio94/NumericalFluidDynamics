/*
 *  .............................. contr ............................
 * 
 * CONTR contours the given field and shows the max/min & title, with box option for nesting.
 * ATMS 502/CSE 566  Spring, 2017
 *
 * >>> NOTE: if you use boundary points (even though passing a NXxNY array to here),
 * then your nestX1,nestX2,nestY1,nestY2 are "off" by BC_WIDTH as far as plotting
 * coordinates for the nest are concerned. If so, call this program with:
 * nestX1-BC_WIDTH,nestX2-BC_WIDTH,nestY1-BC_WIDTH, and nestY2-BC_WIDTH
 * (subtract the BC_WIDTH!) as the last arguments if you want to plot the nest domain.
 * This is because your code really has nestX1 as grid points offset from I1.
 * 
 *   Arguments:
 * 
 * 	nx,ny	input	integers	dimensions of 's' without ghost points
 * 	s	input	float array	field to be contoured.
 * 	cint	input	float		contour interval
 *	simtime	input	float		integration time - real value
 * 	title	input	char *		character string with plot title
 *      colors	input	integer		if =0, positive values colored red, <0 blue;
 *					if >0, colors reversed;
 *					if <0, all colors black (but negative values dashed)
 *	pltzero	input	integer		if zero, the zero contour is plotted.
 *					if nonzero, zero contour is omitted.
 *	nestX1,	input	integers	for nesting : bounding box to highlight on plot;
 *	  nestX2,			  set nestX1 < 0 to skip (if nestX1<0, no box plotted)
 *	  nestY1,			  change box_thickness option below if desired.
 *	  nestY2			If no nesting: set all 4 values to zero.
 *	name	input	char*		character string with your name
 */ 

#define LIWK 3000
#define LRWK 5000

#include <ncarg/ncargC.h>
#include <ncarg/gks.h>
#include <stdlib.h>
#include <math.h>
#define IWTYPE 1
#define WKID   1

void contr(int nx,int ny,float s[][ny],float cint,float simtime,
	   char *title,int colors,int pltzero,int nestX1,int nestX2,
	   int nestY1,int nestY2,char *name)
{

	int i,j,ncon,iwrk[LIWK],s_imin,s_jmin,s_imax,s_jmax;
	float rwrk[LRWK],cval,smin,smax,cmin,cmax,csize;
	float *workarray,*ptr;
        Gcolr_rep rgb1,rgb2,rgb3,rgb4,rgb5;
	char mmlabel[80],tlabel[40];

/*
 * .. Internal controls
 *      If nonzero, "debug" produces text information
 *      label_interval : 0=none, 1=all, 2=every other contour labeled, etc.
 *      high_low_labels: 0=none; 1=H/L only; 2=H/L with value; 3=value only
 *      other variables: high_low_size set the size of the high/low labels;
 *                       cntr_smoothing sets contour spline tension (0: none)
 */

	int debug 	     = 0;
	int label_interval   = 0;
	int high_low_labels  = 2;
	float box_thickness  = 3.0;
        float high_low_size  = 0.025;
	float cntr_smoothing = 0.000;

	/*
	 * Find min,max
	 */

        s_imin = 0; s_jmin = 0;
        s_imax = 0; s_jmax = 0;
	smin = s[0][0];
	smax = s[0][0];
	for (j=0; j<ny; j++) {
	  for (i=0; i<nx; i++) {
	    if (s[i][j] < smin) {
	      smin=s[i][j];
	      s_imin = i;
	      s_jmin = j;
	    }
	    if (s[i][j] > smax) {
	      smax=s[i][j];
	      s_imax = i;
	      s_jmax = j;
	    }
	  }
	}
	if (debug) printf("Contr: %dx%d, min %.2f, max %.2f : %s\n",nx,ny,smin,smax,title);

	/*
	 * Create min/max and time labels
	 */

	if (smin == smax) {
	  csize = 0.014;
	  if (fabs(smin)<999.0) {
	    sprintf(mmlabel,"CONSTANT FIELD = %10.5f",smin);
	  } else {
	    sprintf(mmlabel,"CONSTANT FIELD = %.3f",smin);
	  }

	} else if (smin > -999.0 && smax < 999.0) {
	  csize = 0.014;
	  sprintf(mmlabel,"MIN =%8.3f (%4d,%3d), MAX =%8.3f (%4d,%3d)",
	    smin,s_imin,s_jmin,smax,s_imax,s_jmax);

	} else {
	  csize = 0.013;
	  sprintf(mmlabel,"MIN =%.3f (%4d,%3d), MAX =%.3f (%4d,%3d)",
	    smin,s_imin,s_jmin,smax,s_imax,s_jmax);
	}
	sprintf(tlabel,"TIME=%8.3f",simtime);
	if (debug) {
	  printf("mmlabel '%s',tlabel '%s'\n",mmlabel,tlabel);
	}

	/*
	 * Allocate and fill temporary array with 
	 * dimensions reversed for NCAR routine
	 */

	workarray = (float *) malloc(nx*ny*sizeof(float));
	ptr = workarray;

	for (j=0; j<ny; j++) {
	  for (i=0; i<nx; i++) {
	    *ptr++ = s[i][j];
	  }
	}

	/*
	 * Plot labels
	 */

	c_set(0.,1.,0.,1.,0.,1.,0.,1.,1);
	c_pcmequ(0.50,0.97,title  ,0.020,0.0, 0.0);
	c_pcmequ(0.95,0.03,mmlabel,csize,0.0, 1.0);
	c_pcmequ(0.05,0.03, tlabel,0.014,0.0,-1.0);

        /*
         * Additional labels
         */

        c_pcmequ(0.02,0.94,name,0.01,90.,1.);
        c_pcmequ(0.98,0.06,"ATMS 502/CSE 566",0.01,90.,-1.);
        c_pcmequ(0.98,0.94,"Spring 2017",0.01,90.,1.);

	/*
	 * Set colors for contours
	 */

        rgb1.rgb.red = 1.; rgb1.rgb.green = 1.; rgb1.rgb.blue = 1.;
        rgb5.rgb.red = 0.; rgb5.rgb.green = 0.; rgb5.rgb.blue = 0.;
	if (colors == 0) {
          rgb2.rgb.red = 1.; rgb2.rgb.green = 0.; rgb2.rgb.blue = 0.;
          rgb3.rgb.red = 0.; rgb3.rgb.green = 1.; rgb3.rgb.blue = 0.;
          rgb4.rgb.red = 0.; rgb4.rgb.green = 0.; rgb4.rgb.blue = 1.;
	  if (debug) printf("Positive values contoured red, negative blue.\n");
	} else if (colors > 0) {
          rgb2.rgb.red = 0.; rgb2.rgb.green = 0.; rgb2.rgb.blue = 1.;
          rgb3.rgb.red = 0.; rgb3.rgb.green = 1.; rgb3.rgb.blue = 0.;
          rgb4.rgb.red = 1.; rgb4.rgb.green = 0.; rgb4.rgb.blue = 0.;
	  if (debug) printf("Positive values contoured blue, negative red.\n");
	} else {
          rgb2.rgb.red = 0.; rgb2.rgb.green = 0.; rgb2.rgb.blue = 0.;
          rgb3.rgb.red = 0.; rgb3.rgb.green = 0.; rgb3.rgb.blue = 0.;
          rgb4.rgb.red = 0.; rgb4.rgb.green = 0.; rgb4.rgb.blue = 0.;
	  if (debug) printf("All contours plotted black.\n");
	}
        /* gset_colr_rep(WKID,0,&rgb1); No!! */
        gset_colr_rep(WKID,1,&rgb2);
        gset_colr_rep(WKID,2,&rgb3);
        gset_colr_rep(WKID,3,&rgb4);
        gset_colr_rep(WKID,4,&rgb5);

	/*
	 * Prepare to plot contours.
	 */

        c_cpsetr("CIS - CONTOUR INTERVAL SPECIFIER",cint);
        c_cpseti("LIS - LABEL INTERVAL SPECIFIER"  ,label_interval);
        c_cpseti("HIC - High Label Color Index"    ,4);
        c_cpseti("LOC - Low Label Color Index"     ,4);
        if (high_low_labels == 0) {
          c_cpsetc("HLT","");
        } else if (high_low_labels == 1) {
          c_cpsetc("HLT","H'L");
        } else if (high_low_labels == 2) {
          c_cpsetc("HLT","H:B:$ZDV$:E:'L:B:$ZDV$:E:");
        } else {
          c_cpsetc("HLT","$ZDV$'$ZDV$");
	}
        c_cpsetr("HLS - HIGH/LOW LABEL SIZE",high_low_size);
        c_cpsetr("T2D - Tension 2-dim spline",cntr_smoothing);
        c_cpseti("HLO - High/Low Label Overlap Flag",0);
        /* c_cpsetc("CFT","") */
        c_cpsetc("ILT","");
        if (pltzero != 0 && debug)
          printf("Zero contour omitted; plotting positive values.\n");

	/*
	 * Plot positive (and possibly zero) contours
	 */

        if ( (pltzero == 0 && smax >= 0.0) ||
             (pltzero != 0 && smax >= cint) ||
             (smin == smax && smin == 0.0) ) {
          cmin = cint * (float)( (int)(smin/cint) );
          if (cmin < 0.0) cmin=0.0;
          if (cmin == 0.0 && pltzero != 0) cmin=cint;
          c_cpsetr("CMN - CONTOUR MIN",cmin);
          cmax = cint * (float)( (int)(smax/cint) );
          if (cmax < cmin) cmax=cmin;
          c_cpsetr("CMX - CONTOUR MAX",cmax);
	  c_cprect(workarray,nx,nx,ny,rwrk,LRWK,iwrk,LIWK);
	  /*c_cpback(workarray,rwrk,iwrk);*/
	  c_cppkcl(workarray,rwrk,iwrk);
	  c_cpgeti("NCL",&ncon);
	  if (debug) printf("%d positive contours to be plotted.\n",ncon);

	  /* Set color of each contour depending on <0,=0,>0 */

	  for (i=1; i<=ncon; i++) {
	    c_cpseti("PAI",i);
	    c_cpgetr("CLV",&cval);
	    if (label_interval == 0) c_cpseti("CLU",1);
	    if (cval > 0.0) {
	      c_cpseti("CLC",1);
	    } else if (cval == 0.0) {
	      c_cpseti("CLC",2);
	    } else {
	      c_cpseti("CLC",3);
	      if (colors<0) c_cpseti("CLD",21845);
	    }
	  }

	  c_cplbdr(workarray,rwrk,iwrk);
	  c_cpcldr(workarray,rwrk,iwrk);
	  c_gacolr(4,4,4,4);
	  c_gridal(nx-1,0,ny-1,0,0,0,5,0.,0.);
	}

	/*
	 * Now plot all contours for negative values
	 */

	if (smin < 0) {
	  c_cpsetr("CIS - CONTOUR INTERVAL SPECIFIER",cint);
	  c_cpseti("LIS - LABEL INTERVAL SPECIFIER"  ,label_interval);
          cmin = cint * (float)( (int)(smin/cint) );
          if (cmin > (-cint)) cmin=(-cint);
          if (smax >= 0.0) {
            cmax = -cint;
          } else {
            cmax = cint * (float)( (int)(smax/cint) );
          }
          c_cpsetr("CMN - CONTOUR MIN",cmin);
          c_cpsetr("CMX - CONTOUR MAX",cmax);

	  c_cprect(workarray,nx,nx,ny,rwrk,LRWK,iwrk,LIWK);
	  /*c_cpback(workarray,rwrk,iwrk);*/
	  c_cppkcl(workarray,rwrk,iwrk);
	  c_cpgeti("NCL",&ncon);
	  if (debug) printf("%d negative contours to be plotted.\n",ncon);

	  for (i=1; i<=ncon; i++) {
	    c_cpseti("PAI",i);
	    c_cpgetr("CLV",&cval);
            if (label_interval == 0) c_cpseti("CLU",1);
	    if (cval > 0.0) {
	      c_cpseti("CLC",1);
	    } else if (cval == 0.0) {
	      c_cpseti("CLC",2);
	    } else {
	      c_cpseti("CLC",3);
	      if (colors<0) c_cpseti("CLD",21845);
	    }
	  }

	  c_cplbdr(workarray,rwrk,iwrk);
	  c_cpcldr(workarray,rwrk,iwrk);
	  c_gacolr(4,4,4,4);
	  c_gridal(nx-1,0,ny-1,0,0,0,5,0.,0.);
	}

	/*
	 * Plot bounding box if desired (if nestX1 >= 0)
	 */

	if (nestX1 >= 0) {
	  c_plotif(0.,0.,2);
	  gset_linewidth(box_thickness);
	  c_frstpt((float) nestX1+1,(float) nestY1+1);
	  c_vector((float) nestX2+1,(float) nestY1+1);
	  c_vector((float) nestX2+1,(float) nestY2+1);
	  c_vector((float) nestX1+1,(float) nestY2+1);
	  c_vector((float) nestX1+1,(float) nestY1+1);
	  c_plotif(0.,0.,2);
	  gset_linewidth(1.0);
	}

	/*
	 * Close plot frame.
	 */

	c_frame();

	/*
	 * Restore colors
	 */

        rgb1.rgb.red = 1.; rgb1.rgb.green = 1.; rgb1.rgb.blue = 1.;
        rgb2.rgb.red = 0.; rgb2.rgb.green = 0.; rgb2.rgb.blue = 0.;
        /* gset_colr_rep(WKID,0,&rgb1); No!! */
        gset_colr_rep(WKID,1,&rgb2);

	/*
	 * Free work array and return.
	 */

	free(workarray);
	return;
}
