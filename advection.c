#include <math.h>
/*
void advection(int nxdim, int nzdim,int nx, int nz, int i1,int i2,int j1,int j2,
		float t1[nxdim][nzdim],float t2[nxdim][nzdim],float p1[nxdim][nzdim],float p2[nxdim][nzdim],
		float u1[nx+1][nzdim], float w1[nxdim][nz+1],float u2[nx+1][nzdim], float w2[nxdim][nz+1],
		float u3[nx+1][nzdim], float w3[nxdim][nz+1],
		float dx,float dz,float dt,int n)
*/


void advection(int nxdim, int nydim, int nzdim, float p1[nxdim][nydim][nzdim], float p2[nxdim][nydim][nzdim], float p3[nxdim][nydim][nzdim],
float t1[nxdim][nydim][nzdim], float t2[nxdim][nydim][nzdim],
float u1[nxdim+1][nydim][nzdim],float u2[nxdim+1][nydim][nzdim],float u3[nxdim+1][nydim][nzdim],
float v1[nxdim][nydim+1][nzdim],float v2[nxdim][nydim+1][nzdim],float v3[nxdim][nydim+1][nzdim],
float w1[nxdim][nydim][nzdim+1],float w2[nxdim][nydim][nzdim+1],float w3[nxdim][nydim][nzdim+1], float rhobar[nzdim],
int i1,int i2,int j1,int j2,int k1,int k2,float dx,float dy,float dz, float dt, int n)
{

	
	int i,ii,j,k;
/*
	float s11d[nzdim],s21d[nzdim],u1d[nz+1];
	float s11dx[nxdim],s21dx[nxdim],u1dx[nx+1];

	void bc(int nxdim, int nzdim,int nx, int nz, int i1,int i2,int j1,int j2,
		float t1[nxdim][nzdim],float p1[nxdim][nzdim],float u1[nx+1][nzdim], float w1[nxdim][nz+1],
		float dx,float dz,float dt);
	
	void advect1d(int N,int n, float s1[N],float s2[N],float u[n+1],float dt,float dx, char advection_type);
*/




void bc(int nxdim, int nydim, int nzdim, float p1[nxdim][nydim][nzdim], float p2[nxdim][nydim][nzdim], float p3[nxdim][nydim][nzdim],
float t1[nxdim][nydim][nzdim], float t2[nxdim][nydim][nzdim],
float u1[nxdim+1][nydim][nzdim],float u2[nxdim+1][nydim][nzdim],float u3[nxdim+1][nydim][nzdim],
float v1[nxdim][nydim+1][nzdim],float v2[nxdim][nydim+1][nzdim],float v3[nxdim][nydim+1][nzdim],
float w1[nxdim][nydim][nzdim+1],float w2[nxdim][nydim][nzdim+1],float w3[nxdim][nydim][nzdim+1],
int i1,int i2,int j1,int j2,int k1,int k2);
	

	
void advect1d(int N,int i1, int i2, float s1[N],float s2[N],float u[N+1],float dt,float dz);

float s11d[nzdim],s21d[nzdim],u1d[nzdim+1];
float s11dx[nxdim],s21dx[nxdim],u1dx[nxdim+1];
float s11dy[nydim],s21dy[nydim],u1dy[nydim+1];




for(k=k1;k<=k2;k++){
for(j=j1;j<=j2;j++){
	for(ii=0;ii<nxdim;ii++){s11dx[ii]=t1[ii][j][k];}
	for(ii=0;ii<nxdim+1;ii++){u1dx[ii]=u2[ii][j][k];}
	advect1d(nxdim,i1,i2,s11dx,s21dx,u1dx,dt/2.,dx);
	for(i=i1;i<=i2;i++){t2[i][j][k]=s21dx[i];}
}
}

bc(nxdim,nydim, nzdim,p1,p2,p3,t1,t2,u1,u2,u3,v1,v2,v3,w1,w2,w3,i1,i2,j1,j2,k1,k2);

/*THETA ADVECTION*/
	/*z integration*/
	
	
for(k=k1;k<=k2;k++){
for(i = i1; i <= i2; i++){
		for(ii = 0; ii < nydim; ii++){s11dy[ii]=t2[i][ii][k];}
		for(ii = 0; ii < nydim+1; ii++){u1dy[ii] =  v2[i][ii][k];}
	    advect1d(nydim,j1,j2,s11dy,s21dy,u1dy,dt/2.,dy);
		for(j=j1;j<=j2;j++){t2[i][j][k]=s21dy[j];}
	}
}
	
bc(nxdim,nydim, nzdim,p1,p2,p3,t1,t2,u1,u2,u3,v1,v2,v3,w1,w2,w3,i1,i2,j1,j2,k1,k2);


for(j=j1;j<=j2;j++){
for(i = i1; i <= i2; i++){
		for(ii = 0; ii < nzdim; ii++){s11d[ii]=t2[i][j][ii];}
		for(ii = 0; ii < nzdim+1; ii++){u1d[ii] = w2[i][j][ii];}
	    advect1d(nzdim,k1,k2,s11d,s21d,u1d,dt,dz);
		for(k=k1;k<=k2;k++){t2[i][j][k]=s21d[k];}
	}
}
	
bc(nxdim,nydim, nzdim,p1,p2,p3,t1,t2,u1,u2,u3,v1,v2,v3,w1,w2,w3,i1,i2,j1,j2,k1,k2);


for(k=k1;k<=k2;k++){
for(i = i1; i <= i2; i++){
		for(ii = 0; ii < nydim; ii++){s11dy[ii]=t2[i][ii][k];}
		for(ii = 0; ii < nydim+1; ii++){u1dy[ii] = v2[i][ii][k];}
	    advect1d(nydim,j1,j2,s11dy,s21dy,u1dy,dt/2.,dy);
		for(j=j1;j<=j2;j++){t2[i][j][k]=s21dy[j];}
	}
}
	
bc(nxdim,nydim, nzdim,p1,p2,p3,t1,t2,u1,u2,u3,v1,v2,v3,w1,w2,w3,i1,i2,j1,j2,k1,k2);


for(k=k1;k<=k2;k++){
for(j=j1;j<=j2;j++){
	for(ii=0;ii<nxdim;ii++){s11dx[ii]=t2[ii][j][k];}
	for(ii=0;ii<nxdim+1;ii++){u1dx[ii]=u2[ii][j][k];}
	advect1d(nxdim,i1,i2,s11dx,s21dx,u1dx,dt/2.,dx);
	for(i=i1;i<=i2;i++){t2[i][j][k]=s21dx[i];}
}
}

bc(nxdim,nydim, nzdim,p1,p2,p3,t1,t2,u1,u2,u3,v1,v2,v3,w1,w2,w3,i1,i2,j1,j2,k1,k2);


if(n==1){
dt = dt/2.;
}
	


#pragma omp parallel for private(i,j,k)
	
for(i=i1+1;i<=i2;i++){
	for(j=j1;j<=j2;j++){
		for(k=k1;k<=k2;k++){
u3[i][j][k] = u1[i][j][k] \
-2*dt/(4*dx)*( (u2[i+1][j][k]+u2[i][j][k])*(u2[i+1][j][k]-u2[i][j][k])+(u2[i][j][k]+u2[i-1][j][k])*(u2[i][j][k]-u2[i-1][j][k]) ) \
-2*dt/(4*dy)*( (v2[i][j+1][k]+v2[i-1][j+1][k])*(u2[i][j+1][k]-u2[i][j][k])+(v2[i][j][k]+v2[i-1][j][k])*(u2[i][j][k]-u2[i][j-1][k]) ) \
-2*dt/(4*dz)*( (w2[i][j][k+1]+w2[i-1][j][k+1])*(u2[i][j][k+1]-u2[i][j][k])+(w2[i][j][k]+w2[i-1][j][k])*(u2[i][j][k]-u2[i][j][k-1]) );
				}
		}
}


#pragma omp parallel for private(i,j,k)
	
for(i=i1;i<=i2;i++){
	for(j=j1;j<=j2;j++){
		for(k=k1;k<=k2;k++){
v3[i][j][k] = v1[i][j][k]\
-2*dt/(4*dx)*( (u2[i+1][j][k]+u2[i+1][j-1][k])*(v2[i+1][j][k]-v2[i][j][k])+(u2[i][j][k]+u2[i][j-1][k])*(v2[i][j][k]-v2[i-1][j][k]) ) \
-2*dt/(4*dy)*( (v2[i][j+1][k]+v2[i][j][k])*(v2[i][j+1][k]-v2[i][j][k])+(v2[i][j][k]+v2[i][j-1][k])*(v2[i][j][k]-v2[i][j-1][k]) ) \
-2*dt/(4*dz)*( (w2[i][j][k+1]+w2[i][j-1][k+1])*(v2[i][j][k+1]-v2[i][j][k])+(w2[i][j][k]+w2[i][j-1][k])*(v2[i][j][k]-v2[i][j][k-1]) );
				}
		}
}



#pragma omp parallel for private(i,j,k)
	

for(i=i1;i<=i2;i++){
	for(j=j1;j<=j2;j++){
		for(k=k1+1;k<=k2;k++){
w3[i][j][k] = w1[i][j][k] \
-2*dt/(4*dx)*( (u2[i+1][j][k]+u2[i+1][j][k-1])*(w2[i+1][j][k]-w2[i][j][k])+(u2[i][j][k]+u2[i][j][k-1])*(w2[i][j][k]-w2[i-1][j][k]) )\
-2*dt/(4*dy)*( (v2[i][j+1][k]+v2[i][j+1][k-1])*(w2[i][j+1][k]-w2[i][j][k])+(v2[i][j][k]+v2[i][j][k-1])*(w2[i][j][k]-w2[i][j-1][k]) )\
-2*dt/(4*dz)*( (w2[i][j][k+1]+w2[i][j][k])*(w2[i][j][k+1]-w2[i][j][k])+(w2[i][j][k]+w2[i][j][k-1])*(w2[i][j][k]-w2[i][j][k-1]) );
			}
	}
} 

	
	

	return;
}

