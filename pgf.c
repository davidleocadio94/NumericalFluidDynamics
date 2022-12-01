#include <math.h>


/*
void pgf(int nxdim, int nzdim,int nx, int nz, int i1,int i2,int j1,int j2,
		float t1[nxdim][nzdim],float t2[nxdim][nzdim],float p1[nxdim][nzdim],float p2[nxdim][nzdim], float p3[nxdim][nzdim],
		float u1[nx+1][nzdim], float w1[nxdim][nz+1],float u2[nx+1][nzdim], float w2[nxdim][nz+1],
		float u3[nx+1][nzdim], float w3[nxdim][nz+1],
		float dx,float dz,float dt, float rhobar[nzdim]){
*/
void pgf(int nxdim, int nydim, int nzdim, float p1[nxdim][nydim][nzdim], float p2[nxdim][nydim][nzdim], float p3[nxdim][nydim][nzdim],
float t1[nxdim][nydim][nzdim], float t2[nxdim][nydim][nzdim],
float u1[nxdim+1][nydim][nzdim],float u2[nxdim+1][nydim][nzdim],float u3[nxdim+1][nydim][nzdim],
float v1[nxdim][nydim+1][nzdim],float v2[nxdim][nydim+1][nzdim],float v3[nxdim][nydim+1][nzdim],
float w1[nxdim][nydim][nzdim+1],float w2[nxdim][nydim][nzdim+1],float w3[nxdim][nydim][nzdim+1], float rhobar[nzdim],
int i1,int i2,int j1,int j2,int k1,int k2,float dx,float dy,float dz, float dt,int n)
{

	
	
	float tbar = 300.0;
	float g = 9.81;
	float cs = 60.0;
	float max = 0.0;

int i,j,k;

if(n==1){
dt = dt/2.;
}	
/*u pgf*/	
	
#pragma omp parallel for private(i,j,k)
	
for(i=i1+1;i<=i2;i++){
	for(j=j1;j<=j2;j++){
		for(k=k1;k<=k2;k++){
u3[i][j][k] = u3[i][j][k]-2*dt*(1./rhobar[k]*(p1[i][j][k]-p1[i-1][j][k])/dx);
		}
	}
}


#pragma omp parallel for private(i,j,k)
	

for(i=i1;i<=i2;i++){
	for(j=j1;j<=j2;j++){
		for(k=k1;k<=k2;k++){
v3[i][j][k]=v3[i][j][k] - 2*dt*(1./rhobar[k]*(p1[i][j][k]-p1[i][j-1][k])/dy);	
		}
	}
}


float term1;
float term2;


#pragma omp parallel for private(i,j,k,term1,term2)
	

for(i=i1;i<=i2;i++){
	for(j=j1;j<=j2;j++){
		for(k=k1+1;k<=k2;k++){
term1 = 0.5*(rhobar[k]+rhobar[k-1]);
term1 = 1/term1;
term2 = 0.5*(t1[i][j][k]-tbar)+0.5*(t1[i][j][k-1]-tbar);
term2 = term2/tbar;


w3[i][j][k] = w3[i][j][k]-2*dt*term1/dz*(p1[i][j][k]-p1[i][j][k-1])+2*dt*g*term2;
/*
w3[i][j][k]=w1[i][j][k]-2*dt*(2/(rhobar[k]+rhobar[k-1])*(p1[i][j][k]-p1[i][j][k-1])/dz)\
+2*dt*g*((t1[i][j][k]+t1[i][j][k-1]-2.*tbar)/(2.*tbar));
*/

		}
	}
}




	
bc(nxdim,nydim, nzdim,p1,p2,p3,t1,t2,u1,u2,u3,v1,v2,v3,w1,w2,w3,i1,i2,j1,j2,k1,k2);


#pragma omp parallel for private(i,j,k,term1,term2)	

for(i=i1;i<=i2;i++){
	for(j=j1;j<=j2;j++){
		for(k=k1;k<=k2;k++){
term1 = 0.5*(rhobar[k+1]+rhobar[k])*w3[i][j][k+1];
term2 = 0.5*(rhobar[k]+rhobar[k-1])*w3[i][j][k];
p3[i][j][k]=p1[i][j][k]-2*dt*pow(cs,2)*(rhobar[k]*(u3[i+1][j][k]-u3[i][j][k])/dx+rhobar[k]*(v3[i][j+1][k]-v3[i][j][k])/dy+(term1-term2)/dz);


/*
p3[i][j][k]=p1[i][j][k]-2*dt*pow(cs,2)*(rhobar[k]*(u3[i+1][j][k]-u3[i][j][k])/dx+rhobar[k]*(v3[i][j+1][k]-v3[i][j][k])/dy+\
( 0.5*(rhobar[k+1]+rhobar[k])*w3[i][j][k+1]-0.5*(rhobar[k]+rhobar[k-1])*w3[i][j][k] )/dz);
*/

		}
	}
}




	return;
}
