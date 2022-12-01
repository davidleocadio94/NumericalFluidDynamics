/*
 *  ATMS 502 / CSE 566 -- Spring, 2017
 *  Demo for pgm1:  Linear and nonlinear advection
 *  =====>>>>> PUT YOUR NAME HERE! and in "name" variable below <<<<<=====
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ncarg/ncargC.h>
#include <ncarg/gks.h>
#define IWTYPE 1
#define WKID   1



main()
{

/*
 * Definitions
 */
#define NX 300 
#define NY 300
#define NZ 75
#define BC_WIDTH 2
#define I1 BC_WIDTH
#define I2 I1+NX-1
#define J1 BC_WIDTH
#define J2 J1+NY-1
#define K1 BC_WIDTH
#define K2 K1+NZ-1
#define NXDIM NX+2*BC_WIDTH
#define NYDIM NY+2*BC_WIDTH
#define NZDIM NZ+2*BC_WIDTH

char *name  = "David Villarreal";
/*char *title = "Program 1"; used elsewhere below */

/* Arrays and other variables */

float p1[NXDIM][NYDIM][NZDIM],p2[NXDIM][NYDIM][NZDIM],p3[NXDIM][NYDIM][NZDIM];
float u1[NXDIM+1][NYDIM][NZDIM],u2[NXDIM+1][NYDIM][NZDIM],u3[NXDIM+1][NYDIM][NZDIM];
float v1[NXDIM][NYDIM+1][NZDIM],v2[NXDIM][NYDIM+1][NZDIM],v3[NXDIM][NYDIM+1][NZDIM];
float w1[NXDIM][NYDIM][NZDIM+1],w2[NXDIM][NYDIM][NZDIM+1],w3[NXDIM][NYDIM][NZDIM+1];
float t1[NXDIM][NYDIM][NZDIM], t2[NXDIM][NYDIM][NZDIM];
int i1,i2,j1,j2,k1,k2,nx,ny,nz;
	
	
	
	
	
	float dt,dtt,courant,c,dx,dy,dz,max,min,EDISS,EDIS,ET;
	float Km, Kt;
	float pi = 4.0*atan(1.0);
	float rhobar[NZDIM];
	char plottitle[20];
	
	int i,j,k,n,nstep,nplot;

	
	dx =50.;
	dy = 50.;
	dz = 50.;
	dt = 0.150; 
	nstep = 600;
	Km = 15.0;
	Kt = 2.0;	




/* Variables to reverse default black/white colors in NCAR Graphics */

	Gcolr_rep rgb1,rgb2;

/* Function prototype declarations */

void ic(int nxdim, int nydim, int nzdim, float p1[nxdim][nydim][nzdim], float p2[nxdim][nydim][nzdim], float p3[nxdim][nydim][nzdim],
float t1[nxdim][nydim][nzdim], float t2[nxdim][nydim][nzdim],
float u1[nxdim+1][nydim][nzdim],float u2[nxdim+1][nydim][nzdim],float u3[nxdim+1][nydim][nzdim],
float v1[nxdim][nydim+1][nzdim],float v2[nxdim][nydim+1][nzdim],float v3[nxdim][nydim+1][nzdim],
float w1[nxdim][nydim][nzdim+1],float w2[nxdim][nydim][nzdim+1],float w3[nxdim][nydim][nzdim+1], float rhobar[nzdim],
int i1,int i2,int j1,int j2,int k1,int k2,float dx,float dy,float dz);



void update(int nxdim, int nydim, int nzdim, float p1[nxdim][nydim][nzdim], float p2[nxdim][nydim][nzdim], float p3[nxdim][nydim][nzdim],
float t1[nxdim][nydim][nzdim], float t2[nxdim][nydim][nzdim],
float u1[nxdim+1][nydim][nzdim],float u2[nxdim+1][nydim][nzdim],float u3[nxdim+1][nydim][nzdim],
float v1[nxdim][nydim+1][nzdim],float v2[nxdim][nydim+1][nzdim],float v3[nxdim][nydim+1][nzdim],
float w1[nxdim][nydim][nzdim+1],float w2[nxdim][nydim][nzdim+1],float w3[nxdim][nydim][nzdim+1]);


	void advect1d(int N,int i1, int i2, float s1[N],float s2[N],float u[N+1],float dt,float dz);



	void error(int nx, int ny, float st[][ny], float s[][ny],float ediss,float edis, float et);
	
void advection(int nxdim, int nydim, int nzdim, float p1[nxdim][nydim][nzdim], float p2[nxdim][nydim][nzdim], float p3[nxdim][nydim][nzdim],
float t1[nxdim][nydim][nzdim], float t2[nxdim][nydim][nzdim],
float u1[nxdim+1][nydim][nzdim],float u2[nxdim+1][nydim][nzdim],float u3[nxdim+1][nydim][nzdim],
float v1[nxdim][nydim+1][nzdim],float v2[nxdim][nydim+1][nzdim],float v3[nxdim][nydim+1][nzdim],
float w1[nxdim][nydim][nzdim+1],float w2[nxdim][nydim][nzdim+1],float w3[nxdim][nydim][nzdim+1], float rhobar[nzdim],
int i1,int i2,int j1,int j2,int k1,int k2,float dx,float dy,float dz, float dt,int n);


		
	
void diffusion(int nxdim, int nydim, int nzdim, float p1[nxdim][nydim][nzdim], float p2[nxdim][nydim][nzdim], float p3[nxdim][nydim][nzdim],
float t1[nxdim][nydim][nzdim], float t2[nxdim][nydim][nzdim],
float u1[nxdim+1][nydim][nzdim],float u2[nxdim+1][nydim][nzdim],float u3[nxdim+1][nydim][nzdim],
float v1[nxdim][nydim+1][nzdim],float v2[nxdim][nydim+1][nzdim],float v3[nxdim][nydim+1][nzdim],
float w1[nxdim][nydim][nzdim+1],float w2[nxdim][nydim][nzdim+1],float w3[nxdim][nydim][nzdim+1], float rhobar[nzdim],
int i1,int i2,int j1,int j2,int k1,int k2,float dx,float dy,float dz, float dt,int Km, int Kt, int n);	


void pgf(int nxdim, int nydim, int nzdim, float p1[nxdim][nydim][nzdim], float p2[nxdim][nydim][nzdim], float p3[nxdim][nydim][nzdim],
float t1[nxdim][nydim][nzdim], float t2[nxdim][nydim][nzdim],
float u1[nxdim+1][nydim][nzdim],float u2[nxdim+1][nydim][nzdim],float u3[nxdim+1][nydim][nzdim],
float v1[nxdim][nydim+1][nzdim],float v2[nxdim][nydim+1][nzdim],float v3[nxdim][nydim+1][nzdim],
float w1[nxdim][nydim][nzdim+1],float w2[nxdim][nydim][nzdim+1],float w3[nxdim][nydim][nzdim+1], float rhobar[nzdim],
int i1,int i2,int j1,int j2,int k1,int k2,float dx,float dy,float dz, float dt,int n);



	
	
void bc(int nxdim, int nydim, int nzdim, float p1[nxdim][nydim][nzdim], float p2[nxdim][nydim][nzdim], float p3[nxdim][nydim][nzdim],
float t1[nxdim][nydim][nzdim], float t2[nxdim][nydim][nzdim],
float u1[nxdim+1][nydim][nzdim],float u2[nxdim+1][nydim][nzdim],float u3[nxdim+1][nydim][nzdim],
float v1[nxdim][nydim+1][nzdim],float v2[nxdim][nydim+1][nzdim],float v3[nxdim][nydim+1][nzdim],
float w1[nxdim][nydim][nzdim+1],float w2[nxdim][nydim][nzdim+1],float w3[nxdim][nydim][nzdim+1],
int i1,int i2,int j1,int j2,int k1,int k2);
	



void putfield(int nx,int ny,int nz,char *fieldname,float datatime,float field[nx][ny][nz]);



	
	
	void contr(int nx,int ny,float splot[nx][ny],float cint,float simtime,
           char *title,int colors,int pltzero,int nestX1,int nestX2,
           int nestY1,int nestY2,char *name);
           
	void sfc(int nx,int ny,int nymax,float splot[nx][ny],float simtime,
         float angh,float angv,char *label,char *name);
         
	
/* Plotting declarations */
	int colors,pltzero;
	char title[25];
	float cint,simtime,angh,angv,contint;


/* Parameters and input .................................... */

	printf("Program #6       Numerical Fluid Dynamics\n\n");
	/*printf("NX=%d, BC_WIDTH=%d, I1=%d, I2=%d, NXDIM=%d\n",
		NX,BC_WIDTH,I1X,I2X,NXDIM);*/
		
	
	
		
		
ic(NXDIM,NYDIM,NZDIM,p1,p2,p3,t1,t2,u1,u2,u3,v1,v2,v3,w1,w2,w3,rhobar,I1,I2,J1,J2,K1,K2,dx,dy,dz);

bc(NXDIM,NYDIM,NZDIM,p1,p2,p3,t1,t2,u1,u2,u3,v1,v2,v3,w1,w2,w3,I1,I2,J1,J2,K1,K2);


float tplot[NX][NY][NZ],uplot[NX][NY][NZ],vplot[NX][NY][NZ],wplot[NX][NY][NZ],pplot[NX][NY][NZ];
int NNX = NX+1;
int NNZ = NZ+1;

for(i=I1;i<=I2;i++){
	for(j=J1;j<=J2;j++){
			for(k=K1;k<=K2;k++){
			tplot[i-I1][j-J1][k-K1]=t1[i][j][k]-300.0;
			}
	}
}

putfield(NX,NY,NZ,"T",0,tplot);


for(i=I1;i<=I2;i++){
	for(j=J1;j<=J2;j++){
			for(k=K1;k<=K2;k++){
			pplot[i-I1][j-J1][k-K1]=p1[i][j][k];
			}
	}
}

putfield(NX,NY,NZ,"P",0,pplot);


for(i=I1;i<=I2;i++){
	for(j=J1;j<=J2;j++){
			for(k=K1;k<=K2;k++){
			uplot[i-I1][j-J1][k-K1]=0.5*(u1[i][j][k]+u1[i+1][j][k]);
			}
	}
}

putfield(NX,NY,NZ,"U",0,uplot);


for(i=I1;i<=I2;i++){
	for(j=J1;j<=J2;j++){
			for(k=K1;k<=K2;k++){
			vplot[i-I1][j-J1][k-K1]=0.5*(v1[i][j][k]+v1[i][j+1][k]);
			}
	}
}

putfield(NX,NY,NZ,"V",0,vplot);


for(i=I1;i<=I2;i++){
	for(j=J1;j<=J2;j++){
			for(k=K1;k<=K2;k++){
			wplot[i-I1][j-J1][k-K1]=0.5*(w1[i][j][k]+w1[i][j][k+1]);
			}
	}
}

putfield(NX,NY,NZ,"W",0,wplot);

float t;
int save_interval=15;
int T = 750;
int ntime = T/dt;
int nsave = (int)rintf(save_interval/dt);


for(n=1;n<=ntime;n++){


dtt = dt;
advection(NXDIM,NYDIM,NZDIM,p1,p2,p3,t1,t2,u1,u2,u3,v1,v2,v3,w1,w2,w3,rhobar,I1,I2,J1,J2,K1,K2,dx,dy,dz,dtt,n);
dtt=dt;

diffusion(NXDIM,NYDIM,NZDIM,p1,p2,p3,t1,t2,u1,u2,u3,v1,v2,v3,w1,w2,w3,rhobar,I1,I2,J1,J2,K1,K2,dx,dy,dz,dtt,Km,Kt,n);
dtt=dt;

pgf(NXDIM,NYDIM,NZDIM,p1,p2,p3,t1,t2,u1,u2,u3,v1,v2,v3,w1,w2,w3,rhobar,I1,I2,J1,J2,K1,K2,dx,dy,dz,dtt,n);

bc(NXDIM,NYDIM,NZDIM,p1,p2,p3,t1,t2,u1,u2,u3,v1,v2,v3,w1,w2,w3,I1,I2,J1,J2,K1,K2);

update(NXDIM,NYDIM,NZDIM,p1,p2,p3,t1,t2,u1,u2,u3,v1,v2,v3,w1,w2,w3);
if(n%nsave==0){
t +=save_interval;

for(i=I1;i<=I2;i++){
	for(j=J1;j<=J2;j++){
			for(k=K1;k<=K2;k++){
			tplot[i-I1][j-J1][k-K1]=t1[i][j][k]-300.0;
			}
	}
}

putfield(NX,NY,NZ,"T",t,tplot);


for(i=I1;i<=I2;i++){
	for(j=J1;j<=J2;j++){
			for(k=K1;k<=K2;k++){
			pplot[i-I1][j-J1][k-K1]=p2[i][j][k];
			}
	}
}

putfield(NX,NY,NZ,"P",t,pplot);


for(i=I1;i<=I2;i++){
	for(j=J1;j<=J2;j++){
			for(k=K1;k<=K2;k++){
			uplot[i-I1][j-J1][k-K1]=0.5*(u2[i][j][k]+u2[i+1][j][k]);
			}
	}
}

putfield(NX,NY,NZ,"U",t,uplot);


for(i=I1;i<=I2;i++){
	for(j=J1;j<=J2;j++){
			for(k=K1;k<=K2;k++){
			vplot[i-I1][j-J1][k-K1]=0.5*(v2[i][j][k]+v2[i][j+1][k]);
			}
	}
}

putfield(NX,NY,NZ,"V",t,vplot);


for(i=I1;i<=I2;i++){
	for(j=J1;j<=J2;j++){
			for(k=K1;k<=K2;k++){
			wplot[i-I1][j-J1][k-K1]=0.5*(w2[i][j][k]+w2[i][j][k+1]);
			}
	}
}

putfield(NX,NY,NZ,"W",t,wplot);


}


}

/*
 * Open the NCAR Graphics package and set colors.
 */
	gopen_gks("stdout",0);
	gopen_ws(WKID, NULL, IWTYPE);
	gactivate_ws(WKID);

	/* omit following four lines to invert black/white colors */
	rgb1.rgb.red = 1.; rgb1.rgb.green = 1.; rgb1.rgb.blue = 1.;
	rgb2.rgb.red = 0.; rgb2.rgb.green = 0.; rgb2.rgb.blue = 0.;
        gset_colr_rep(WKID,0,&rgb1);
        gset_colr_rep(WKID,1,&rgb2);

 /*start program 5*/
 
/*
        putfield(nx,ny,nz,"T",(float)1.0,t);
*/	
 

 

	gdeactivate_ws(WKID);
	gclose_ws(WKID);
	gclose_gks();

	printf("Program finished normally.\n");
	exit;
}
