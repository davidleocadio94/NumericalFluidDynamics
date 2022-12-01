#include <math.h>
/*advect1d(s11d,s21d,u1d,advection_type);*/
void advect1d(int N,int i1, int i2, float s1[N],float s2[N],float u[N+1],float dt,float dz){

int i,j,is;
float courant;
int iflux;
float flux[N+1];

/*


for(i = 0 ; i < n;i++){
u[i] = 0.5*(u[i+1]+u[i]);
}




	for(i=1;i<=N-2;i++){
	courant = u[i]*dt/dz;

	s2[i] = s1[i]-courant/2.*(s1[i+1]-s1[i-1])+pow(courant,2)/2.*\
	    (s1[i+1]-2.*s1[i]+s1[i-1]);

	}




*/

for(i=0;i<N+1;i++){
flux[i]=0.0;
}





for(i=i1;i<=i2+1;i++){

courant = fabs(dt/dz*u[i]);


if(u[i] >= 0.0){
flux[i] = courant*(s1[i-1]+0.5*(1-courant)*0.5*(s1[i]-s1[i-2]));

}else if(u[i]<0.0){
flux[i] = courant*(-s1[i]+0.5*(1-courant)*0.5*(s1[i+1]-s1[i-1]));

}

}





for(i=i1;i<=i2;i++){
s2[i] = s1[i]-(flux[i+1]-flux[i])+dt/dz*s1[i]*(u[i+1]-u[i]);
}




return;

}
