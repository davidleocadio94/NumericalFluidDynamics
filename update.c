

void update(int nxdim, int nydim, int nzdim, float p1[nxdim][nydim][nzdim], float p2[nxdim][nydim][nzdim], float p3[nxdim][nydim][nzdim],
float t1[nxdim][nydim][nzdim], float t2[nxdim][nydim][nzdim],
float u1[nxdim+1][nydim][nzdim],float u2[nxdim+1][nydim][nzdim],float u3[nxdim+1][nydim][nzdim],
float v1[nxdim][nydim+1][nzdim],float v2[nxdim][nydim+1][nzdim],float v3[nxdim][nydim+1][nzdim],
float w1[nxdim][nydim][nzdim+1],float w2[nxdim][nydim][nzdim+1],float w3[nxdim][nydim][nzdim+1])
{

int i,j,k;

/*t and p*/
for(i=0;i<nxdim;i++){
	for(j=0;j<nydim;j++){
		for(k=0;k<nzdim;k++){
	/*t1[i][j]=t2[i][j];*/
	t1[i][j][k]=t2[i][j][k];
	p1[i][j][k]=p2[i][j][k];
	p2[i][j][k]=p3[i][j][k];
		}
	}
}

/*u*/

for(i=0;i<nxdim+1;i++){
	for(j=0;j<nydim;j++){
		for(k=0;k<nzdim;k++){
	u1[i][j][k]=u2[i][j][k];
	u2[i][j][k]=u3[i][j][k];
		}
	}
}

/*w*/

for(i=0;i<nxdim;i++){
	for(j=0;j<nydim+1;j++){
		for(k=0;k<nzdim;k++){
	v1[i][j][k]=v2[i][j][k];
	v2[i][j][k]=v3[i][j][k];
		}
	}
}

for(i=0;i<nxdim;i++){
	for(j=0;j<nydim;j++){
		for(k=0;k<nzdim+1;k++){
	w1[i][j][k]=w2[i][j][k];
	w2[i][j][k]=w3[i][j][k];
		}
	}
}









/*
for(i=0;i<nxdim;i++){
	for(j=0;j<nzdim;j++){
	t2[i][j]=0.0;
	p3[i][j]=0.0;
	}
}



for(i=0;i<nx+1;i++){
	for(j=0;j<nzdim;j++){
	u3[i][j]=0.0;
}
}



for(i=0;i<nxdim;i++){
	for(j=0;j<nz+1;j++){
	w3[i][j]=0.0;
}
}
*/
return;
}

