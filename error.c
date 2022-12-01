


void error(int nx, int ny, float st[][ny], float s[][ny],float ediss,float edis, float et)
 {

int M = (nx)*(ny);
int i,j;
et = 0.;

float sigman2=0.;
float sigmat2=0.;
float meant = 0.;
float meann =0.;
float rho1 = 0.;
float rho2 = 0.;
float rho3 = 0.;
float rho;

for(i = 0; i < nx; i++){
	for(j = 0; j<nx; j++){
	et = et+pow((st[i][j]-s[i][j]),2);
	meann = meann + s[i][j];
	meant = meant + st[i][j];
	}
}

meann = 1/((float)M)*meann;
meant = 1/((float)M)*meant;
et = 1/((float)M)*et;
printf("TOTAL ERROR = %.6f\n",et);

for(i = 0; i < nx; i++){
	for(j = 0; j<ny; j++){
	rho1 = rho1+(s[i][j]-meann)*(st[i][j]-meant);
	rho2 = rho2+pow(s[i][j]-meann,2);
	rho3 = rho3+pow(st[i][j]-meant,2);
	sigman2 = sigman2 + pow(s[i][j]-meann,2);
	sigmat2 = sigmat2+pow(st[i][j]-meant,2);
	}
}


sigman2 = 1/(float)(M)*sigman2;
sigmat2 = 1/(float)(M)*sigmat2;
rho = rho1/(float)sqrt(rho2*rho3);

ediss = pow(sqrt(sigmat2)-sqrt(sigman2),2)+pow(meant-meann,2);
printf("DISSIPATION ERROR =%.5f\n",ediss);
edis = 2*(1-rho)*sqrt(sigmat2)*sqrt(sigman2);
printf("DISPERSION ERROR =%.5f\n",edis);
printf("dispersion+dissipation = %.5f\n",edis+ediss);

return;

}






 