#include <math.h>


void diffusion(int nxdim, int nydim, int nzdim, float p1[nxdim][nydim][nzdim], float p2[nxdim][nydim][nzdim], float p3[nxdim][nydim][nzdim],
float t1[nxdim][nydim][nzdim], float t2[nxdim][nydim][nzdim],
float u1[nxdim+1][nydim][nzdim],float u2[nxdim+1][nydim][nzdim],float u3[nxdim+1][nydim][nzdim],
float v1[nxdim][nydim+1][nzdim],float v2[nxdim][nydim+1][nzdim],float v3[nxdim][nydim+1][nzdim],
float w1[nxdim][nydim][nzdim+1],float w2[nxdim][nydim][nzdim+1],float w3[nxdim][nydim][nzdim+1], float rhobar[nzdim],
int i1,int i2,int j1,int j2,int k1,int k2,float dx,float dy,float dz, float dt,int Km, int Kt, int n)
{	
int i,j,k;
	

if(n==1){
dt = dt/2;
}

float diffx,diffy,diffz;

#pragma omp parallel for private(i,j,k,diffx,diffy,diffz)
	
for(i=i1+1;i<=i2;i++){
	for(j=j1;j<=j2;j++){
		for(k=k1;k<=k2;k++){

diffx = (u1[i+1][j][k]-2*u1[i][j][k]+u1[i-1][j][k])/pow(dx,2);
diffy = (u1[i][j+1][k]-2*u1[i][j][k]+u1[i][j-1][k])/pow(dy,2);
diffz = (u1[i][j][k+1]-2*u1[i][j][k]+u1[i][j][k-1])/pow(dz,2);

u3[i][j][k] = u3[i][j][k] + 2*dt*Km*(diffx+diffy+diffz);
		}
	}
}

#pragma omp parallel for private(i,j,k,diffx,diffy,diffz)
	
for(i=i1;i<=i2;i++){
	for(j=j1;j<=j2;j++){
		for(k=k1;k<=k2;k++){

	
diffx = (v1[i+1][j][k]-2*v1[i][j][k]+v1[i-1][j][k])/pow(dx,2);
diffy = (v1[i][j+1][k]-2*v1[i][j][k]+v1[i][j-1][k])/pow(dy,2);
diffz = (v1[i][j][k+1]-2*v1[i][j][k]+v1[i][j][k-1])/pow(dz,2);

v3[i][j][k] = v3[i][j][k] + 2*dt*Km*(diffx+diffy+diffz);
		}
	}
}
 

#pragma omp parallel for private(i,j,k,diffx,diffy,diffz)
	
for(i=i1;i<=i2;i++){
	for(j=j1;j<=j2;j++){
		for(k=k1+1;k<=k2;k++){

diffx = (w1[i+1][j][k]-2*w1[i][j][k]+w1[i-1][j][k])/pow(dx,2);
diffy = (w1[i][j+1][k]-2*w1[i][j][k]+w1[i][j-1][k])/pow(dy,2);
diffz = (w1[i][j][k+1]-2*w1[i][j][k]+w1[i][j][k-1])/pow(dz,2);

w3[i][j][k] = w3[i][j][k] + 2*dt*Km*(diffx+diffy+diffz);
		}
	}
}

if(n==1){
dt = dt*2;
}

#pragma omp parallel for private(i,j,k,diffx,diffy,diffz)
	
for(i=i1;i<=i2;i++){
	for(j=j1;j<=j2;j++){
		for(k=k1;k<=k2;k++){

	
diffx = (t1[i+1][j][k]-2*t1[i][j][k]+t1[i-1][j][k])/pow(dx,2);
diffy = (t1[i][j+1][k]-2*t1[i][j][k]+t1[i][j-1][k])/pow(dy,2);
diffz = (t1[i][j][k+1]-2*t1[i][j][k]+t1[i][j][k-1])/pow(dz,2);

t2[i][j][k] = t2[i][j][k] + dt*Kt*(diffx+diffy+diffz);
		}
	}
}



 
			
		return;
}
