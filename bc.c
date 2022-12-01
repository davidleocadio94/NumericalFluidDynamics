/*
 * ============================ bc =====================
 * BC sets the boundary conditions
 * Arguments:
 *
 *	s1	real array	values at current time step
 *	i1,i2	integers	indices bounding array data
 *	nx	integer		main array size, not including
 *				extra 'ghost' zones/points
 */

/* BFJ add args */
/*void bc(int NXDIM, int NYDIM, int i1,int i2,int j1,int j2, float s[NXDIM][NYDIM])*/

void bc(int NXDIM, int NYDIM, int NZDIM, float p1[NXDIM][NYDIM][NZDIM], float p2[NXDIM][NYDIM][NZDIM], float p3[NXDIM][NYDIM][NZDIM],
float t1[NXDIM][NYDIM][NZDIM],float t2[NXDIM][NYDIM][NZDIM],
float u1[NXDIM+1][NYDIM][NZDIM],float u2[NXDIM+1][NYDIM][NZDIM],float u3[NXDIM+1][NYDIM][NZDIM],
float v1[NXDIM][NYDIM+1][NZDIM],float v2[NXDIM][NYDIM+1][NZDIM],float v3[NXDIM][NYDIM+1][NZDIM],
float w1[NXDIM][NYDIM][NZDIM+1],float w2[NXDIM][NYDIM][NZDIM+1],float w3[NXDIM][NYDIM][NZDIM+1],
int I1,int I2,int J1,int J2,int K1,int K2)
{
    	int i,j,k;

    /*** U:  X (Symmetry ... but asymmetry for U) Boundaries ***/
	/*#pragma omp parallel for shared(u1,u2,u3) private(i,j,k)*/


#pragma omp parallel for private(i,j,k)
	
    	for (j=J1; j<=J2; j++) { 
      		for (k=K1; k<=K2; k++) { 
        		u1[I1][j][k]   = -u1[I1+1][j][k];
        		u1[I2+1][j][k] = -u1[I2][j][k];
        		u2[I1][j][k]   = -u2[I1+1][j][k];
        		u2[I2+1][j][k] = -u2[I2][j][k];
        		u3[I1][j][k]   = -u3[I1+1][j][k];
        		u3[I2+1][j][k] = -u3[I2][j][k];
		}
	}

    /*** U:  Z (0-gradient) Boundaries ***/
	/*#pragma omp parallel for shared(u1,u2,u3) private(i,j,k)*/


#pragma omp parallel for private(i,j,k)
	
    	for (i=I1; i<=I2+1; i++) {
      		for (j=J1; j<=J2; j++) {
        		u1[i][j][K1-1] = u1[i][j][K1];
          		u1[i][j][K2+1] = u1[i][j][K2];
        		u2[i][j][K1-1] = u2[i][j][K1];
        		u2[i][j][K2+1] = u2[i][j][K2];
        		u3[i][j][K1-1] = u3[i][j][K1];
        		u3[i][j][K2+1] = u3[i][j][K2];
      		} 
    	}
	
	/*** U:  Y (Periodic) Boundaries ***/
	/*#pragma omp parallel for shared(u1,u2,u3) private(i,j,k)*/


#pragma omp parallel for private(i,j,k)
	
    	for (i=I1; i<=I2+1; i++) {
      		for (k=K1; k<=K2; k++) {
        		u1[i][J1-1][k] = u1[i][J2][k];
        		u1[i][J2+1][k] = u1[i][J1][k];
        		u2[i][J1-1][k] = u2[i][J2][k];
        		u2[i][J2+1][k] = u2[i][J1][k];
        		u3[i][J1-1][k] = u3[i][J2][k];
        		u3[i][J2+1][k] = u3[i][J1][k];
      		}  
    	}

    /*** W:  Z (Rigid upper/lower lid) Boundaries ***/

/*	#pragma omp parallel for shared(w1,w2,w3) private(i,j,k)*/
	

#pragma omp parallel for private(i,j,k)
	
	for (i=I1; i<=I2; i++) {
      		for (j=J1; j<=J2; j++) {
        		w1[i][j][K1]   = 0.0;
        		w1[i][j][K2+1] = 0.0;
        		w2[i][j][K1]   = 0.0;
        		w2[i][j][K2+1] = 0.0;
        		w3[i][j][K1]   = 0.0;
        		w3[i][j][K2+1] = 0.0;
      		}
    	}

    /*** W:  X (Symmetry) Boundaries ***/
/*	#pragma omp parallel for shared(w1,w2,w3) private(i,j,k)*/
 

#pragma omp parallel for private(i,j,k)
	   
	for (j=J1; j<=J2; j++) {
      		for (k=K1; k<=K2+1; k++) {
        		w1[I1-1][j][k] = w1[I1+1][j][k];
        		w1[I2+1][j][k] = w1[I2-1][j][k];
        		w2[I1-1][j][k] = w2[I1+1][j][k];
        		w2[I2+1][j][k] = w2[I2-1][j][k];
        		w3[I1-1][j][k] = w3[I1+1][j][k];
        		w3[I2+1][j][k] = w3[I2-1][j][k];					
      		}
    	}

    /*** W:  Y (Periodic) Boundaries ***/
/*	#pragma omp parallel for shared(w1,w2,w3) private(i,j,k)*/


#pragma omp parallel for private(i,j,k)
	
	for (i=I1; i<=I2; i++) {
      		for (k=K1; k<=K2+1; k++) {
        		w1[i][J1-1][k] = w1[i][J2][k];
        		w1[i][J2+1][k] = w1[i][J1][k];
        		w2[i][J1-1][k] = w2[i][J2][k];
        		w2[i][J2+1][k] = w2[i][J1][k];
        		w3[i][J1-1][k] = w3[i][J2][k];
        		w3[i][J2+1][k] = w3[i][J1][k];
      		}
    	}

    /*** P:  Z (0-gradient) Boundaries ***/
/*	#pragma omp parallel for shared(p1,p2,p3) private(i,j,k)*/


#pragma omp parallel for private(i,j,k)
	
    	for (i=I1; i<=I2; i++) {
      		for (j=J1; j<=J2; j++) {
        		p1[i][j][K1-1] = p1[i][j][K1];
        		p1[i][j][K2+1] = p1[i][j][K2];
        		p2[i][j][K1-1] = p2[i][j][K1];
        		p2[i][j][K2+1] = p2[i][j][K2];
        		p3[i][j][K1-1] = p3[i][j][K1];
        		p3[i][j][K2+1] = p3[i][j][K2];
        		
        		t1[i][j][K1-1] = t1[i][j][K1];
			t1[i][j][K1-2] = t1[i][j][K1];
        		t1[i][j][K2+1] = t1[i][j][K2];
			t1[i][j][K2+2] = t1[i][j][K2];
        		t2[i][j][K1-1] = t2[i][j][K1];
			t2[i][j][K1-2] = t2[i][j][K1];
        		t2[i][j][K2+1] = t2[i][j][K2];
			t2[i][j][K2+2] = t2[i][j][K2];
        		
        		
      		}
    	}

     /*** P: X (Symmetry) Boundaries ***/
/*	#pragma omp parallel for shared(p1,p2,p3) private(i,j,k)*/


#pragma omp parallel for private(i,j,k)
	
    	for (j=J1; j<=J2; j++) {
      		for (k=K1; k<=K2; k++) {
        		p1[I1-1][j][k] = p1[I1+1][j][k];
        		p1[I2+1][j][k] = p1[I2-1][j][k];
        		p2[I1-1][j][k] = p2[I1+1][j][k];
        		p2[I2+1][j][k] = p2[I2-1][j][k];
        		p3[I1-1][j][k] = p3[I1+1][j][k];
        		p3[I2+1][j][k] = p3[I2-1][j][k];
        		
        		
        		t1[I1-1][j][k] = t1[I1+1][j][k];
			t1[I1-2][j][k] = t1[I1+2][j][k];
        		t1[I2+1][j][k] = t1[I2-1][j][k];
			t1[I2+2][j][k] = t1[I2-2][j][k];
        		t2[I1-1][j][k] = t2[I1+1][j][k];
			t2[I1-2][j][k] = t2[I1+2][j][k];
        		t2[I2+1][j][k] = t2[I2-1][j][k];
			t2[I2+2][j][k] = t2[I2-2][j][k];
        		
        		
        		
         	}
    	}

    /*** P: Y (Periodic) Boundaries ***/
/*	#pragma omp parallel for shared(p1,p2,p3) private(i,j,k)*/


#pragma omp parallel for private(i,j,k)
	
    	for (i=I1; i<=I2; i++) {
      		for (k=K1; k<=K2; k++) {
        		p1[i][J1-1][k] = p1[i][J2][k];
        		p1[i][J2+1][k] = p1[i][J1][k];
        		p2[i][J1-1][k] = p2[i][J2][k];
        		p2[i][J2+1][k] = p2[i][J1][k];
        		p3[i][J1-1][k] = p3[i][J2][k];
        		p3[i][J2+1][k] = p3[i][J1][k];
        		
        		
        		
        		t1[i][J1-1][k] = t1[i][J2][k];
			t1[i][J1-2][k] = t1[i][J2-1][k];
        		t1[i][J2+1][k] = t1[i][J1][k];
			t1[i][J2+2][k] = t1[i][J1+1][k];
        		t2[i][J1-1][k] = t2[i][J2][k];
			t2[i][J1-2][k] = t2[i][J2-1][k];
        		t2[i][J2+1][k] = t2[i][J1][k];
			t2[i][J2+2][k] = t2[i][J1+1][k];
      		}
    	}

     /*** V: Z (0-Gradient) Boundaries ***/
/*	#pragma omp parallel for shared(v1,v2,v3) private(i,j,k)*/


#pragma omp parallel for private(i,j,k)
	
    	for (i=I1; i<=I2; i++) {
     	 	for (j=J1; j<=J2+1; j++) {
        		v1[i][j][K1-1] = v1[i][j][K1];
        		v1[i][j][K2+1] = v1[i][j][K2];
        		v2[i][j][K1-1] = v2[i][j][K1];
        		v2[i][j][K2+1] = v2[i][j][K2];
        		v3[i][j][K1-1] = v3[i][j][K1];
        		v3[i][j][K2+1] = v3[i][j][K2];
      		}
    	}

    /*** V: X (Symmetry) Boundaries ***/
/*	#pragma omp parallel for shared(v1,v2,v3) private(i,j,k)*/


#pragma omp parallel for private(i,j,k)
	
    	for (j=J1; j<=J2+1; j++) {
      		for (k=K1; k<=K2; k++) {
        		v1[I1-1][j][k] = v1[I1+1][j][k];
        		v1[I2+1][j][k] = v1[I2-1][j][k];
        		v2[I1-1][j][k] = v2[I1+1][j][k];
        		v2[I2+1][j][k] = v2[I2-1][j][k];
        		v3[I1-1][j][k] = v3[I1+1][j][k];
        		v3[I2+1][j][k] = v3[I2-1][j][k];
      		}
    	}

    /*** V: Y (Periodic) Boundaries ***/
/*	#pragma omp parallel for shared(v1,v2,v3) private(i,j,k)*/


#pragma omp parallel for private(i,j,k)
	
    	for (i=I1; i<=I2; i++) {
      		for (k=K1; k<=K2; k++) {
        		v1[i][J1-1][k] = v1[i][J2][k];
        		v1[i][J2+1][k] = v1[i][J1][k];
        		v2[i][J1-1][k] = v2[i][J2][k];
        		v2[i][J2+1][k] = v2[i][J1][k];
        		v3[i][J1-1][k] = v3[i][J2][k];
        		v3[i][J2+1][k] = v3[i][J1][k];
      		}
    	}

    	return;
}

