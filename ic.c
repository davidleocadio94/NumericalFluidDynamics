/*
 * ============================ ic =====================
 * IC sets the initial condition
 * Arguments:
 *
 *	s1	real array	IC data. Set i1..i2 here;
 *				  [i1-1],[i2+1] = ghost zones
 *				  if 1 ghost point on each side
 *	dx	real		grid spacing
 *	i1,i2	integers	indices bounding array data
 *	nx	integer		number of physical grid points
 *
 */
#include <math.h>
#include <stdlib.h>

/*
void ic(int nxdim, int nzdim,int nx, int nz, int i1,int i2,int j1,int j2,
		float t1[nxdim][nzdim],float t2[nxdim][nzdim],float p1[nxdim][nzdim],float p2[nxdim][nzdim],float p3[nxdim][nzdim],
		float u1[nx+1][nzdim], float w1[nxdim][nz+1],float u2[nx+1][nzdim], float w2[nxdim][nz+1],
		float u3[nx+1][nzdim], float w3[nxdim][nz+1],float rhobar[nzdim],
		float dx,float dz,float dt)*/
		
void ic(int nxdim, int nydim, int nzdim, float p1[nxdim][nydim][nzdim], float p2[nxdim][nydim][nzdim], float p3[nxdim][nydim][nzdim],
float t1[nxdim][nydim][nzdim], float t2[nxdim][nydim][nzdim],
float u1[nxdim+1][nydim][nzdim],float u2[nxdim+1][nydim][nzdim],float u3[nxdim+1][nydim][nzdim],
float v1[nxdim][nydim+1][nzdim],float v2[nxdim][nydim+1][nzdim],float v3[nxdim][nydim+1][nzdim],
float w1[nxdim][nydim][nzdim+1],float w2[nxdim][nydim][nzdim+1],float w3[nxdim][nydim][nzdim+1], float rhobar[nzdim],
int i1,int i2,int j1,int j2,int k1,int k2,float dx,float dy,float dz)
{
	int i,j,k,m;
	float pi = 4.0*atan(1.0);
	float tbar = 300.0;
	float pert[2]={-10.,-10.};
	float pertv[2]={-20.,20.};
	float x,y,z;
	float x0[2]={75.,14975.};
	float y0[2] ={7525.,7525.};
	float z0[2]={1225.,1225.};
	float rm;
	float xradius=2000.;
	float yradius=999999.;
	float zradius=2000.;
	
	
/*
#pragma omp parallel for private(i,j,k)
*/	
	for(i=0;i<nxdim;i++){
		for(j=0;j<nydim;j++){
			for(k=0;k<nzdim;k++){
				t1[i][j][k]=tbar;
				p1[i][j][k]=0.0;
				p2[i][j][k]=0.0;
				p3[i][j][k]=0.0;
			}
		}
	}
	
	
/*
#pragma omp parallel for private(i,j,k)
*/	
	for(i=0;i<nxdim+1;i++){
		for(j=0;j<nydim;j++){
			for(k=0;k<nzdim;k++){
				u1[i][j][k]=0.0;
				u2[i][j][k]=0.0;
				u3[i][j][k]=0.0;

			}
		}
	}
	
	
/*
#pragma omp parallel for private(i,j,k)
*/	
	for(i=0;i<nxdim;i++){
		for(j=0;j<nydim+1;j++){
			for(k=0;k<nzdim;k++){
				v1[i][j][k]=0.0;
				v2[i][j][k]=0.0;
				v3[i][j][k]=0.0;

			}
		}
	}
	
	
/*
#pragma omp parallel for private(i,j,k)
*/	
	for(i=0;i<nxdim;i++){
		for(j=0;j<nydim;j++){
			for(k=0;k<nzdim+1;k++){
				w1[i][j][k]=0.0;
				w2[i][j][k]=0.0;
				w3[i][j][k]=0.0;

			}
		}
	}
	
	
	
	
	
	
	
	
/*
#pragma omp parallel for private(i,j,k)
	*/


	for(i=i1;i<=i2;i++){
		for(j=j1;j<=j2;j++){
			for(k=k1;k<=k2;k++){
		x=dx/2.0+dx*(float)(i-i1);
		y=dy/2.0+dy*(float)(j-j1);
		z=dz/2.0+dz*(float)(k-k1);
			for(m=0;m<2;m++){
			rm = sqrt( pow((x-x0[m])/xradius,2.0)+pow((y-y0[m])/yradius,2.0)+pow((z-z0[m])/zradius,2.0));
				if(rm<=1.0){
				t1[i][j][k] = t1[i][j][k] + 0.5*pert[m]*(cos(rm*pi)+1);
				v1[i][j][k] = v1[i][j][k] + 0.5*pertv[m]*(cos(rm*pi)+1);
					}
				}
			}
		}
	}
	
	
	for(k=k1;k<=k2;k++){
	z = dz/2+dz*(float)(k-k1);
	rhobar[k]=pow(10,5)/(287*pow(300,1004./287))*pow(300-9.81/1004*z,1004./287-1.);
	}
	
	rhobar[k1-1]=rhobar[k1];
	rhobar[k2+1]=rhobar[k2];
	
	for(k=k1-1;k<k2+1;k++){
		printf("rhobar %.2f \n",rhobar[k]);
}	
	
/*
#pragma omp parallel for private(i,j,k)
*/	
	
for(i=i1;i<=i2;i++){
	for(j=j1;j<=j2;j++){
		for(k=k1;k<=k2;k++){
		t2[i][j][k]=t1[i][j][k];	
		}
	}
}


/*
#pragma omp parallel for private(i,j,k)
*/	
for(i=0;i<nxdim;i++){
	for(j=0;j<nydim+1;j++){
		for(k=0;k<nzdim+1;k++){
v2[i][j][k]=v1[i][j][k];
		}
	}
}



	
float upertur=2.;
srand(0.0);

/*
#pragma omp parallel for private(i,j,k)
*/	
for(i=i1+1;i<=i2;i++){
	for(j=j1;j<=j2;j++){
		for(k=k1;k<=k2;k++){
u1[i][j][k] = upertur * ( (float)rand()/(RAND_MAX+1.0) )-upertur*0.5;
		}
	}
}		

/*
#pragma omp parallel for private(i,j,k)
*/	

for(i=0;i<nxdim+1;i++){
	for(j=0;j<nydim;j++){
		for(k=0;k<nzdim;k++){
u2[i][j][k] = u1[i][j][k];
		}
	}
}					
							

	return;
}

