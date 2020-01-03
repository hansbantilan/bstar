
/** I still need to work out 1) score(done), 2) the g_ij               **/
/** 3) la density 4) mass, and  5) la frequency                        **/
/************************************************************************/
/** the expansion employed is phi= phi1 cos(wt) + phi3 cos(3wt)+...    **/
/** A=A0+A2 cos(2wt)+A4 cos(4wt)+... C= C0+C2 cos(2wt)+C4 cos(4wt)+... **/
/** xi=phi' = xi1 cos(wt) + xi3 cos(3wt)+... truncated at the maximum  **/
/** value of N . We identify the vector :                              **/
/** y= (y1,y2,...,yN)= (A0,c0,c2,....,phi1,phi3,...xi1,xi3,...)        **/

#include <stdio.h>
#include <math.h>
#include "nrutil.c"
#include "nrutil.h"
#define PI (acos(-1.0))
#define N 12  /** maximum size of the series **/
#define EPS1 1.0e-7

/** Used in lnsrch.c **/
#define ALF 1.0e-4
#define TOLX 1.0e-7

/** used in newt.c **/
#define MAXITS 200
#define TOLF 1.0e-4
#define TOLMIN 1.0e-6
/**#define TOLX 1.0e-7 **/
#define STPMX 100.0

FILE *fp1;

FILE *fp2;
FILE *fp3;
FILE *fp4;
FILE *fp5;
FILE *fp6;
FILE *fpA;
FILE *fpC;
FILE *fpB;
FILE *fpP;
FILE *fpX;
FILE *ftest;

/** N2 Number of boundary conditions in second x2 **/
int nvar, N2, nnn,nn;
double x1,x2,*fvec;
double A[N+1];
float GAMA, phi1;  /** you are asked this as input **/
double wt1, wt2;


double fminx(double x[])
{
  extern void (*nrfuncv)(int n, double v[], double f[]);
  
  int i;
  double sum;
  
  (*nrfuncv)(nn,x,fvec);
  for (sum=0.0,i=1;i<=nn;i++) sum += SQR(fvec[i]);
  return (0.5*sum);
}



int main(void)
{
  void odeint(double ystart[], int nvar, double x1, double x2, double eps, 
	      double h1,double hmin, int *nok, int *nbad,
	      void (*derivs)(double, double [], double []),
	      void (*rkqs)(double [], double [], int, double *,double, 
			   double, double [], double *, double *, 
			   void (*)(double, double [], double [])));  
  void rkqs(double y[], double dydx[], int n, double *x, 
	    double htry, double eps, double yscal[], 
	    double *hdid, double *hnext,
	    void (*derivs)(double, double [], double []));
  void load(double x1, double v[], double y[]);
  void derivs(double x, double y[], double dydx[]);
  void newt(double x[], int n, int *check,
	    void (*vecfunc)(int, double [], double []));
  void shoot(int n, double v[], double f[]);
  void input4(double v[]); 
  void input6(double v[]); 
  void input8(double v[]);
  void input10(double v[]); 
  void input12(double v[]); 
  
  void printvs(double x2, double v[]);
  void imprimir(double x, double y[]);
  void matriz(double x, double y[]);
  void mass_freq(float xx, double y[]);
  
  int check,i,nok,nbad;
  double q1,*v,h1,hmin,*y,*f;
  float xmax, DELTA, xf, dex, wt_1, wt_2;
  
  printf("You are about to integrate the Eintein-Klein-Gordon\n");
  printf("Write the maximun number in the expansion\n");
  printf("This integer number must be even 4<= N <= 12\n"); 
  scanf("%d",&nnn);
  printf("Now write the input GAMA, and xf: first upper limit\n");
  printf("delx: increment length until reaching xmax\n");
  printf("So you are asked to write GAMA,xf,dex,xmax\n");
  scanf("%f %f %f %f",&GAMA,&xf,&dex,&xmax);

  printf("De el tamano del intervalo en que apareceran los resultados\n");
  printf("Write the interval's size in which the result will appear\n");
  printf("We reccomend a value >= 0.01 not smaller\n");
  scanf("%f",&DELTA);
  printf("We will write the results each %f\n",DELTA);
  printf("Enter two fractions of PI: 0.0,0.25,0.5,0.75, etc you choose\n");
  printf("they will be understood as wt= 0, PI/4, PI/2,3PI/2, etc\n");
  printf("Scalar density rho(x,wt) will be ploted at those values\n\n");
  printf("Enter those two values now: \n");
  scanf("%f %f",&wt_1, &wt_2);
  printf("The scalar density will be plotted at wt= %f PI and %f PI\n",
	 wt_1,wt_2);
  wt1= (double) wt_1*PI;
  wt2= (double) wt_2*PI;
  
  if( (nnn % 2) != 0 || nnn < 2 || nnn > 12) {
    printf("You must have writen an even integer 4<= N <= 12\n");
    printf("Try again, exiting\n"); exit(1);}
  N2=nnn;
  f=dvector(1,N2);
  v=dvector(1,N2);
  
  if(nnn==4) input4(v); if(nnn==6) input6(v); if(nnn==8) input8(v);
  if(nnn==10) input10(v); if(nnn==12) input12(v);
  
  printf("In main phi1 is %f\n",phi1);
  
  nvar= 2 + (3*nnn)/2;
  y=dvector(1,nvar);   

  x1 = 0.0;
  x2= (double) xf;


  
  while( x2 <= xmax){
    newt(v,N2,&check,shoot);	
    if (check) {
      printf("shoot failed, bad initial guess\n");
    } else {
      printvs(x2,v);
    }
    x2 += dex;
  }	
  /*** Here descripancies f[i] are computed and ***/
  /*** shown on screen, they must be small      ***/
 
  shoot(N2,v,f);
  for(i=1;i<=N2;i++) printf("f[%d] = %e\n",i,f[i]);

  

  load(x1,v,y);
  
  fp1= fopen("general_a202.dat","w");
  fp2= fopen("general_c202.dat","w");
  fp3= fopen("general_phi202.dat","w");
  fp4= fopen("general_xi202.dat","w");
  fp5=fopen("general_rho202.dat","w");
  fp6=fopen("general_rho202.dat","w");
  fpA=fopen("general_A202.dat","w");
  fpB=fopen("general_B202.dat","w");
  fpC=fopen("general_C202.dat","w");
  fpP=fopen("general_PHI202.dat","w");
  fpX=fopen("general_XI202.dat","w");
  

  for(i=1;i<=N;i++) A[i]=0.0;
  A[0]=y[1];
  imprimir(0.0,y);
  for(x2= DELTA; x2<= xmax; x2 += DELTA){
    h1=1.0e-5; hmin=0.0;
    odeint(y,nvar,x1,x2,EPS1,h1,hmin,&nok,&nbad,derivs,rkqs);
    matriz(x2,y);
    imprimir(x2,y);
    x1=x2;
  }
 
  mass_freq(xmax,y);
  
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  fclose(fp4);
  fclose(fp5);
  fclose(fp6);
  fclose(fpA);
  fclose(fpB);
  fclose(fpC);
  fclose(fpP);
  fclose(fpX);
  


   if(nnn==2) {
  ftest=fopen("Phidata.dat","w");
  fprintf(ftest,"%d\n %d\n %f\n %f\n %f\n %f\n %f\n %f\n",nnn+2,0,xf,dex,xmax,0.01,0,0.5);
  for(i=1;i<=nnn/2+1;i++) fprintf(ftest,"%e\n",v[i]);
  fprintf(ftest,"%e\n",v[nnn/2+1]/10);
  fprintf(ftest,"%e\n",phi1);
  for(i=nnn/2+2;i<=nnn;i++) fprintf(ftest,"%e\n",v[i]);
  fprintf(ftest,"%e\n",v[nnn]/10);
  fclose(ftest);
    }
    if(nnn==4) {
  ftest=fopen("Phidata.dat","w");
  fprintf(ftest,"%d\n %d\n %f\n %f\n %f\n %f\n %f\n %f\n",nnn+2,0,xf,dex,xmax,0.01,0,0.5);
  for(i=1;i<=nnn/2+1;i++) fprintf(ftest,"%e\n",v[i]);
  fprintf(ftest,"%e\n",v[nnn/2+1]/10);
  fprintf(ftest,"%e\n",phi1);
  for(i=nnn/2+2;i<=nnn;i++) fprintf(ftest,"%e\n",v[i]);
  fprintf(ftest,"%e\n",v[nnn]/10);
  fclose(ftest);
  ftest=fopen("Phi0.2828Deg4","w");
  fprintf(ftest,"%d\n %d\n %f\n %f\n %f\n %f\n %f\n %f\n",nnn,0,xf,dex,xmax,0.01,0,0.5);
  for(i=1;i<=nnn/2+1;i++) fprintf(ftest,"%e\n",v[i]);
  fprintf(ftest,"%e\n",phi1);
  for(i=nnn/2+2;i<=nnn;i++) fprintf(ftest,"%e\n",v[i]);
  fclose(ftest);
    }
    if(nnn==6) {
  ftest=fopen("Phidata.dat","w");
  fprintf(ftest,"%d\n %d\n %f\n %f\n %f\n %f\n %f\n %f\n",nnn+2,0,xf,dex,xmax,0.01,0,0.5);
  for(i=1;i<=nnn/2+1;i++) fprintf(ftest,"%e\n",v[i]);
  fprintf(ftest,"%e\n",v[nnn/2+1]/10);
  fprintf(ftest,"%e\n",phi1);
  for(i=nnn/2+2;i<=nnn;i++) fprintf(ftest,"%e\n",v[i]);
  fprintf(ftest,"%e\n",v[nnn]/10);
  fclose(ftest);
    }
    if(nnn==8) {
  ftest=fopen("Phidata.dat","w");
  fprintf(ftest,"%d\n %d\n %f\n %f\n %f\n %f\n %f\n %f\n",nnn+2,0,xf,dex,xmax,0.01,0,0.5);
  for(i=1;i<=nnn/2+1;i++) fprintf(ftest,"%e\n",v[i]);
  fprintf(ftest,"%e\n",v[nnn/2+1]/10);
  fprintf(ftest,"%e\n",phi1);
  for(i=nnn/2+2;i<=nnn;i++) fprintf(ftest,"%e\n",v[i]);
  fprintf(ftest,"%e\n",v[nnn]/10);
  fclose(ftest);
    }
    if(nnn==10) {
  ftest=fopen("Phidata.dat","w");
  fprintf(ftest,"%d\n %d\n %f\n %f\n %f\n %f\n %f\n %f\n",nnn+2,0,xf,dex,xmax,0.01,0,0.5);
  for(i=1;i<=nnn/2+1;i++) fprintf(ftest,"%e\n",v[i]);
  fprintf(ftest,"%e\n",v[nnn/2+1]/10);
  fprintf(ftest,"%e\n",phi1);
  for(i=nnn/2+2;i<=nnn;i++) fprintf(ftest,"%e\n",v[i]);
  fprintf(ftest,"%e\n",v[nnn]/10);
  fclose(ftest);
    }
    if(nnn==12) {
  ftest=fopen("Phidata.dat","w");
  fprintf(ftest,"%d\n %d\n %f\n %f\n %f\n %f\n %f\n %f\n",nnn,0,xf,dex,xmax,0.01,0,0.5);
  for(i=1;i<=nnn/2+1;i++) fprintf(ftest,"%e\n",v[i]);
  fprintf(ftest,"%e\n",phi1);
  for(i=nnn/2+2;i<=nnn;i++) fprintf(ftest,"%e\n",v[i]);
  fclose(ftest);
  ftest=fopen("Phi0.2828Deg4","w");
  fprintf(ftest,"%d\n %d\n %f\n %f\n %f\n %f\n %f\n %f\n",4,0,xf,dex,xmax,0.01,0,0.5);
  for(i=1;i<=3;i++) fprintf(ftest,"%e\n",v[i]);
  fprintf(ftest,"%e\n",phi1);
  fprintf(ftest,"%e\n",v[4]);
  fclose(ftest);
    }
  
  free_dvector(v,1,N2);
  free_dvector(y,1,nvar);
  return 0;
}

void mass_freq(float xx,double y[])
{
  extern int nnn;
  extern double A[N+1];
  int k,j,jmax;
  double a[N+1],c[N+1],phi[N+1],xi[N+1],freq,mass,sum_a,sum_c;
  
  jmax = nnn/2;
  
  for(j=0;j<=N;j++) a[j]=c[j]=phi[j]=xi[j]=0.0;
  a[0]=y[1];
  for(j=1;j<=jmax;j++) a[2*j]=A[2*j];
  for(j=0,k=2;j <= jmax; j++, k +=1) c[2*j]= y[k];
  
  /** computing mass and frequency **/
  for(sum_a=0.0, sum_c=0.0, k=0;k<=nnn;k += 2){
    sum_a += a[k];
    sum_c += c[k];
  }
  mass= 0.5*xx*(1.0-1.0/sum_a);
  freq= sqrt(sum_c)/sum_a;     

  printf("The mass is: %e, and the frequency = %e\n",mass,freq);
}

void load(double x1, double v[], double y[])
{
  extern float phi1;  
  /*** there are other parts where phi1 changes **/
  /** and I should check if here it has the right value **/
  
  
  if(nnn==2){
    y[1]=1.0; 
    y[2]= v[1]; y[3]= v[2]; 
    y[4]= phi1; 
    y[5]=0.0;
  }
  
  if(nnn==4){
    y[1]=1.0; 
    y[2]= v[1]; y[3]= v[2]; y[4]= v[3]; 
    y[5]= phi1; 
    y[6]= v[4]; 
    y[7]= y[8]= 0.0;
  }
  
  if(nnn==6){
    y[1]=1.0; 
    y[2]= v[1]; y[3]= v[2]; y[4]= v[3]; y[5]= v[4]; 
    y[6]= phi1; 
    y[7]= v[5]; y[8]= v[6];
    y[9]=y[10]=y[11]=0.0;
  }
  
  if(nnn==8){
    y[1]=1.0; 
    y[2]= v[1]; y[3]= v[2] ; y[4]= v[3]; y[5]= v[4]; y[6]=v[5]; 
    y[7]= phi1; 
    y[8]= v[6]; y[9]= v[7]; y[10]= v[8];
    y[11]=y[12]=y[13]=y[14]=0.0;
  }
  
  if(nnn==10){
    y[1]=1.0; 
    y[2]= v[1]; y[3]= v[2]; y[4]= v[3]; y[5]= v[4]; y[6]=v[5]; y[7]=v[6];
    y[8]= phi1; 
    y[9]= v[7]; y[10]= v[8]; y[11]= v[9]; y[12]= v[10];
    y[13]=y[14]=y[15]=y[16]=y[17]=0.0;
  }
  
  if(nnn==12){
    y[1]=1.0; 
    y[2]=v[1];y[3]=v[2];y[4]=v[3];y[5]=v[4];y[6]=v[5];y[7]=v[6];y[8]=v[7];
    y[9]= phi1; 
    y[10]=v[8];y[11]=v[9]; y[12]=v[10]; y[13]=v[11]; y[14]=v[12];
    y[15]=y[16]=y[17]=y[18]=y[19]=y[20]=0.0;
  }
  
}

void score(double xf, double y[], double f[])
{
    extern int nnn;
    int jmax, j, k;
    double c[N+1],phi[N+1],xi[N+1];
    
    jmax= nnn/2;

    for(j=0,k=2;j <= jmax; j++, k +=1) c[2*j]= y[k];
    for(j=1;j <= jmax; j++, k +=1) phi[2*j-1]= y[k];
    for(j=1;j <= jmax; j++, k +=1) xi[2*j-1]= y[k];

    if(nnn==2) {
	f[1]= phi[1]; f[2]= c[2];
    }

    if(nnn==4){
	f[1]= phi[1]; f[2]= phi[3]; f[3]= c[2]; f[4]=c[4];
    }

    if(nnn==6){
	f[1]= phi[1]; f[2]= phi[3]; f[3]= phi[5]; 
	f[4]=c[2]; f[5]=c[4]; f[6]=c[6];
    }

    if(nnn==8){
	f[1]= phi[1]; f[2]= phi[3]; f[3]= phi[5]; f[4]= phi[7]; 
	f[5]=c[2]; f[6]=c[4]; f[7]=c[6]; f[8]=c[8];     
    }

    if(nnn==10){
	f[1]= phi[1];f[2]= phi[3];f[3]= phi[5];f[4]= phi[7];f[5]=phi[9]; 
	f[6]=c[2]; f[7]=c[4]; f[8]=c[6]; f[9]=c[8]; f[10]=c[10];     
    }

    if(nnn==12){
	f[1]= phi[1]; f[2]= phi[3]; f[3]= phi[5]; f[4]= phi[7]; 
	f[5]=phi[9]; f[6]=phi[11]; 
	f[7]=c[2]; f[8]=c[4]; f[9]=c[6]; f[10]=c[8]; 
	f[11]=c[10]; f[12]=c[12];     
    }
   
}

void derivs(double x, double y[], double dydx[])
{
    int k,j, jmax;
    void phidot2(double phi[], double t[]);
    void phi2d(double phi[], double phidd[]);
    void cdot_phidot(double c[],double phi[], double cgot[]);
    void nonlinear(double a[],double c[], double pr[]); 
    void cuadrado(double a[],double pr[]);
    void matriz(double x, double yout[]);

    /** y = (a0,c,phi,xi) **/
    double a[N+1],c[N+1], phi[N+1], xi[N+1],phidd[N+1],t[N+1];
    double ac[N+1], xi2[N+1], a2[N+1], phi2[N+1],axi[N+1],cgot[N+1];
    double uno[N+1],dos[N+1], tres[N+1],cinco[N+1], cinco_prim[N+1];
    double seis[N+1],siete[N+1],ocho[N+1], phi3[N+1], phi4[N+1];
    double f1[N+1],f2[N+1],f3[N+1],f4[N+1], siete_prim[N+1];
    double tres_prim[N+1], ocho_prim[N+1];
    extern double A[N+1];
    extern float GAMA;
    extern int nnn;

    jmax= nnn/2;
    matriz(x,y);

    for(j=0;j<=N;j++) a[j]=c[j]=phi[j]=xi[j]=0.0;

    for(j=0;j <= jmax;j++) a[2*j]= A[2*j];
    for(j=0,k=2;j <= jmax; j++, k +=1) c[2*j]= y[k];
    for(j=1;j <= jmax; j++, k +=1) phi[2*j-1]= y[k];
    for(j=1;j <= jmax; j++, k +=1) xi[2*j-1]= y[k];

    if( fabs(x) <= 1.0e-15 ) {
      cuadrado(phi,phi2);
      nonlinear(phi,phi2,phi3); 	
      phi2d(phi,phidd);
      nonlinear(c,phidd,seis);
      cdot_phidot(c,phi,cgot);
      nonlinear(a,phi,siete);  
      nonlinear(a,phi3,siete_prim);
      
      for(k=1; k<=nnn-1; k += 2) 
	f4[k]= seis[k] + 0.5*cgot[k] + siete[k] + GAMA*siete_prim[k];
      
      dydx[1]= 0.0;
      for(k=2; k<=(nnn/2)+2; k +=1) dydx[k]=0.0;	
      for(k=(nnn/2)+3; k<= nnn+2; k +=1) dydx[k]=0.0;
      for(k=nnn+3, j=1; k<= 3*(nnn/2)+2; k+=1, j +=2) {
	if(j>(nnn-1)) {
	  printf("something wrong in derivs exiting\n"); 
	  exit(1);}
	dydx[k]= f4[j];
      }       
    }
    else {
      phidot2(phi,t);
      phi2d(phi,phidd); /** d^2(phi)/dt^2 **/
      cuadrado(a,a2);
      cuadrado(xi,xi2);
      cuadrado(phi,phi2);
      cuadrado(phi2,phi4); /* */
      nonlinear(phi,phi2,phi3); /** **/
      nonlinear(a,c,ac);
      nonlinear(a,xi,axi);
      nonlinear(ac,t,uno);
      nonlinear(a,xi2,dos);
      nonlinear(a2,phi2,tres);
      nonlinear(ac,phi2,cinco);
      nonlinear(c,phidd,seis); 
      cdot_phidot(c,phi,cgot);
      nonlinear(a,phi,siete);  
      nonlinear(axi,phi2,ocho);
      /*** new elements for the quartic term ***/
      nonlinear(a2,phi4,tres_prim);
      nonlinear(ac,phi4,cinco_prim);
      nonlinear(a,phi3,siete_prim);
      nonlinear(axi,phi4,ocho_prim);
      /*This is the tt component of Einstein equation  */
      /*dAdx = 0.5*x* ( A*C*Dphidt^2 + A*dphidx^2 + A^2*phi^2  +  0   */	 
      f1[0]=   0.5*x* (    uno[0]    +   dos[0]   +  tres[0]    +  0.5*GAMA*tres_prim[0])+
      /*    A/x - A^2/x    */ 
	 a[0]/x - a2[0]/x ;
      /* This is the rr compenent of Einstein equation */
      /*                  dC/dx= (2.0/x)*( C  - A*C  */
      for(k=0;k<=nnn;k++) f2[k]= (2.0/x)*(c[k]-ac[k])+ 
      /*		    x*A*phi^2*C	 + 0 		*/
			    x*cinco[k]   + 0.5*GAMA*x*cinco_prim[k];
      /* This is the Klein Gordon equation              */
      for(k=1;k<=nnn-1;k++){
	f3[k]= xi[k];
      /*d^2phidx^2 =  C * d^2phi/dt^2 +  0.5*C*phi   + A*phi     + 0    */	
	f4[k]      =  seis[k]         +  0.5*cgot[k] + siete[k]  + GAMA*siete_prim[k]
      /*  - dphidx/x + 0.5*x*A*dphidx*phi^2  +             0            - A*dphidx/x */
	  - xi[k]/x  + 0.5*x*ocho[k]         +  0.25*GAMA*x*ocho_prim[k]- axi[k]/x;
      }
      
      dydx[1]=f1[0];
      
      for(k=2, j=0; k<=(nnn/2)+2; k +=1, j+=2) {dydx[k]=f2[j];	
      if(j>(nnn)) {printf("1 something wrong in derivs exiting\n");exit(1);}}
      
      for(k=(nnn/2)+3,j=1; k<=nnn+2; k +=1,j +=2) {dydx[k]=f3[j];
      if(j>(nnn-1)){printf("2 something wrong in derivs exiting\n");exit(1);}}
      
      for(k=nnn+3,j=1;k<=3*(nnn/2)+2;k +=1,j +=2) {dydx[k]= f4[j];
      if(j>(nnn-1)){printf("3 something wrong in derivs exiting\n");exit(1);}}  
    }
}


/** Product of two sums with n=0,....N **/

void nonlinear(double a[],double c[], double pr[])
{
    extern int nnn;    
    int j,p,k;
    double sum1,sum2;
    
    for(k=0;k<=nnn;k++){
	for(j=0,sum1=0.0,sum2=0.0;j<= nnn;j++){
	    for(p=0;p<= nnn;p++){
		if(abs(j-p)==k) sum1 += a[j]*c[p];
		if((j+p)==k) sum2 += a[j]*c[p];
	    }      
	}
	pr[k]= 0.5*sum1 + 0.5*sum2;
    }
}


void cuadrado(double a[],double pr[])
{
    extern int nnn;
    int j,p,k;
    double sum1,sum2;
    
    for(k=0;k<=nnn;k++){
	for(j=0,sum1=0.0,sum2=0.0;j<=nnn;j++){
	    for(p=0;p<= nnn;p++){
		if(abs(j-p)==k) sum1 += a[j]*a[p];
		if((j+p)==k) sum2 += a[j]*a[p];
	    }      
	}
	pr[k]= 0.5*sum1 + 0.5*sum2;
    }
}

void phi2d(double phi[], double phidd[])
{
  extern int nnn;
  int j;
    
    phidd[0]=0.0;
    for(j=1;j<=nnn;j++)
	phidd[j]= -j*j*phi[j];
}

/** la sumatoria de phi empieza desde n=1 **/
 
void phidot2(double phi[], double t[])
{
    extern int nnn;
    int j,p,k;
    double sum1,sum2;
    
    for(j=1,sum1=0.0;j<=nnn;j++) sum1 += j*j*phi[j]*phi[j];
    t[0]= 0.5*sum1;
    
    for(k=1;k<=nnn;k++){
	for(j=1,sum1=0.0,sum2=0.0;j<=nnn;j++){
	    for(p=1;p<= nnn;p++){
		if(abs(j-p)==k) sum1 += j*p*phi[j]*phi[p];
		if( (j+p)== k ) sum2 += j*p*phi[j]*phi[p];
	    }      
	}
	t[k]= 0.5*sum1 - 0.5*sum2;
    }
}

void cdot_phidot(double c[],double phi[],double cgot[])
{
    extern int nnn;
    int j,p,k;
    double sum1,sum2;
    
    for(j=1,sum1=0.0;j<=nnn;j++) sum1 += j*j*c[j]*phi[j];
    cgot[0]= 0.5*sum1;
    
    for(k=1;k<=nnn;k++){
	for(j=1,sum1=0.0,sum2=0.0;j<= nnn;j++){
	    for(p=1;p<=nnn;p++){
		if(abs(j-p)==k) sum1 += j*p*c[j]*phi[p];
		if( (j+p)== k ) sum2 += j*p*c[j]*phi[p];
	    }      
	}
	cgot[k]= 0.5*sum1 - 0.5*sum2;
    }
}

void lnsrch(int n, double xold[], double fold, double g[], double p[], 
double x[],double *f, double stpmax, int *check, double (*func)(double []))
{
    int i;
    double a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,sum,temp,
	test,tmplam;
    *check=0;
    for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
    sum=sqrt(sum);
    if (sum > stpmax)
	for (i=1;i<=n;i++) p[i] *= stpmax/sum;
    for (slope=0.0,i=1;i<=n;i++)
	slope += g[i]*p[i];
    if (slope >= 0.0) nrerror("Roundoff problem in lnsrch.");
    test=0.0;
    for (i=1;i<=n;i++) {
	temp=fabs(p[i])/FMAX(fabs(xold[i]),1.0);
	if (temp > test) test=temp;
    }
    alamin=TOLX/test;
    alam=1.0;
    for (;;) {
	for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];
	*f=(*func)(x);
	if (alam < alamin) {
	    for (i=1;i<=n;i++) x[i]=xold[i];
	    *check=1;
	    return;
	} else if (*f <= fold+ALF*alam*slope) return;
	else {
	    if (alam == 1.0)
		tmplam = -slope/(2.0*(*f-fold-slope));
	    else {
		rhs1 = *f-fold-alam*slope;
		rhs2=f2-fold-alam2*slope;
		a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
		b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
		if (a == 0.0) tmplam = -slope/(2.0*b);
		else {
		    disc=b*b-3.0*a*slope;
		    if (disc < 0.0) tmplam=0.5*alam;
		    else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
		    else tmplam=-slope/(b+sqrt(disc));
		}
		if (tmplam > 0.5*alam)
		    tmplam=0.5*alam;
	    }
	}
	alam2=alam;
	f2 = *f;
	alam=FMAX(tmplam,0.1*alam);
    }
}


void (*nrfuncv)(int n, double v[], double f[]);

/***********************
#define FREERETURN {free_vector(fvec,1,n);free_vector(xold,1,n);
 free_vector(p,1,n);free_vector(g,1,n);free_matrix(fjac,1,n,1,n);
 free_ivector(indx,1,n);return;}********/

void newt(double x[], int n, int *check,
	  void (*vecfunc)(int, double [], double []))
{
    void fdjac(int n, double x[], double fvec[], double **df,
	       void (*vecfunc)(int, double [], double []));
    void lnsrch(int n, double xold[], double fold, double g[], double p[], 
		double x[], double *f, double stpmax, int *check, 
		double (*func)(double []));
    void lubksb(double **a, int n, int *indx, double b[]);
    void ludcmp(double **a, int n, int *indx, double *d);
    int i,its,j,*indx;
    double d,den,f,fold,stpmax,sum,temp,test,**fjac,*g,*p,*xold;

    indx=ivector(1,n);
    fjac=dmatrix(1,n,1,n);
    g=dvector(1,n);
    p=dvector(1,n);
    xold=dvector(1,n);
    fvec=dvector(1,n);
    nn=n;
    nrfuncv=vecfunc;
    f=fminx(x); 
    test=0.0;

    for (i=1;i<=n;i++)
	if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
    if (test < 0.01*TOLF) {
	*check=0;     
	/*FREERETURN*/
	{free_dvector(fvec,1,n);free_dvector(xold,1,n);
	free_dvector(p,1,n);free_dvector(g,1,n);free_dmatrix(fjac,1,n,1,n);
	free_ivector(indx,1,n);return;}
    }
    for (sum=0.0,i=1;i<=n;i++) sum += SQR(x[i]);
    stpmax=STPMX*FMAX(sqrt(sum),(double)n);
    for (its=1;its<=MAXITS;its++) {
	fdjac(n,x,fvec,fjac,vecfunc);
	for (i=1;i<=n;i++) {
	    for (sum=0.0,j=1;j<=n;j++) sum += fjac[j][i]*fvec[j];
	    g[i]=sum;
	}
	for (i=1;i<=n;i++) xold[i]=x[i];
	fold=f;
	for (i=1;i<=n;i++) p[i] = -fvec[i];
	ludcmp(fjac,n,indx,&d);
	lubksb(fjac,n,indx,p);
	lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,fminx);
	test=0.0;
	for (i=1;i<=n;i++)
	    if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
	if (test < TOLF) {
	    *check=0;     
	    /*FREERETURN*/
	   { free_dvector(fvec,1,n);free_dvector(xold,1,n);
	    free_dvector(p,1,n);free_dvector(g,1,n);free_dmatrix(fjac,1,n,1,n);
	    free_ivector(indx,1,n);return;}

		}
	if (*check) {
	    test=0.0;
	    den=FMAX(f,0.5*n);
	    for (i=1;i<=n;i++) {
		temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
		if (temp > test) test=temp;
	    }
	    *check=(test < TOLMIN ? 1 : 0);
	    /*FREERETURN*/  
	    {free_dvector(fvec,1,n);free_dvector(xold,1,n);
	    free_dvector(p,1,n);free_dvector(g,1,n);free_dmatrix(fjac,1,n,1,n);
	    free_ivector(indx,1,n);return;}
	    
		}
	test=0.0;
	for (i=1;i<=n;i++) {
	    temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
	    if (temp > test) test=temp;
	}
	if (test < TOLX) 
	    /*FREERETURN*/  
	    {free_dvector(fvec,1,n);free_dvector(xold,1,n);
	    free_dvector(p,1,n);free_dvector(g,1,n);free_dmatrix(fjac,1,n,1,n);
	    free_ivector(indx,1,n);return;}
    }
    nrerror("MAXITS exceeded in newt");
}


#define EPS 1.0e-4
void fdjac(int n, double x[], double fvec[], double **df,
	void (*vecfunc)(int, double [], double []))
{
    int i,j;
    double h,temp,*f;
    
    f=dvector(1,n);
    for (j=1;j<=n;j++) {
	temp=x[j];
	h=EPS*fabs(temp);
	if (h == 0.0) h=EPS;
	x[j]=temp+h;
	h=x[j]-temp;
	(*vecfunc)(n,x,f);
	x[j]=temp;
	for (i=1;i<=n;i++) df[i][j]=(f[i]-fvec[i])/h;
    }
    free_dvector(f,1,n);
}
#undef EPS


extern void (*nrfuncv)(int n, double v[], double f[]);



#define EPS 1.0e-6

int kmax,kount;
double *xp,**yp,dxsav;

void shoot(int n, double v[], double f[])
{
    void derivs(double x, double y[], double dydx[]);
    void load(double x1, double v[], double y[]);
    void odeint(double ystart[], int nvar, double x1, double x2,
		double eps, double h1, double hmin, int *nok, int *nbad,
		void (*derivs)(double, double [], double []),
		void (*rkqs)(double [], double [], int, double *, double, 
			     double, double [], double *, double *, 
			     void (*)(double, double [], double [])));
    void rkqs(double y[], double dydx[], int n, double *x,
		double htry, double eps, double yscal[], double *hdid, 
	      double *hnext,
	      void (*derivs)(double, double [], double []));
    void score(double xf, double y[], double f[]);
    int i,nbad,nok;
    double h1,hmin=0.0,*y;

    y=dvector(1,nvar);
    kmax=0;
    h1=(x2-x1)/100.0;
    load(x1,v,y);
    odeint(y,nvar,x1,x2,EPS,h1,hmin,&nok,&nbad,derivs,rkqs);
    score(x2,y,f); 
    free_dvector(y,1,nvar);
}
#undef EPS


void lubksb(double **a, int n, int *indx, double b[])
{
    int i,ii=0,ip,j;
    double sum;
    
    for (i=1;i<=n;i++) {
	ip=indx[i];
	sum=b[ip];
	b[ip]=b[i];
	if (ii)
	    for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
	else if (sum) ii=i;
	b[i]=sum;
    }
    for (i=n;i>=1;i--) {
	sum=b[i];
	for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
	b[i]=sum/a[i][i];
    }
}


#define TINY 1.0e-20
void ludcmp(double **a, int n, int *indx, double *d)
{
    int i,imax,j,k;
    double big,dum,sum,temp;
    double *vv;
    
    vv=dvector(1,n);
    *d=1.0;
    for (i=1;i<=n;i++) {
	big=0.0;
	for (j=1;j<=n;j++)
	    if ((temp=fabs(a[i][j])) > big) big=temp;
	if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
	vv[i]=1.0/big;
    }
    for (j=1;j<=n;j++) {
	for (i=1;i<j;i++) {
	    sum=a[i][j];
	    for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
	    a[i][j]=sum;
	}
	big=0.0;
	for (i=j;i<=n;i++) {
	    sum=a[i][j];
	    for (k=1;k<j;k++)
		sum -= a[i][k]*a[k][j];
	    a[i][j]=sum;
	    if ( (dum=vv[i]*fabs(sum)) >= big) {
		big=dum;
		imax=i;
	    }
	}
	if (j != imax) {
	    for (k=1;k<=n;k++) {
		dum=a[imax][k];
		a[imax][k]=a[j][k];
		a[j][k]=dum;
	    }
	    *d = -(*d);
	    vv[imax]=vv[j];
	}
	indx[j]=imax;
	if (a[j][j] == 0.0) a[j][j]=TINY;
	if (j != n) {
	    dum=1.0/(a[j][j]);
	    for (i=j+1;i<=n;i++) a[i][j] *= dum;
	}
    }
    free_dvector(vv,1,n);
}
#undef TINY


#define MAXSTP 10000
#define TINY 1.0e-30

extern int kmax,kount;
extern double *xp,**yp,dxsav;

void odeint(double ystart[], int nvar, double x1, double x2, double eps, 
	    double h1, double hmin, int *nok, int *nbad,
	void (*derivs)(double, double [], double []),
	void (*rkqs)(double [], double [], int, double *, double, 
		     double, double [], double *, double *, 
		     void (*)(double, double [], double [])))
{
    int nstp,i;
    double xsav,x,hnext,hdid,h;
    double *yscal,*y,*dydx;
    
    yscal=dvector(1,nvar);
    y=dvector(1,nvar);
    dydx=dvector(1,nvar);
    x=x1;
    h=SIGN(h1,x2-x1);
    *nok = (*nbad) = kount = 0;
    for (i=1;i<=nvar;i++) y[i]=ystart[i];
    if (kmax > 0) xsav=x-dxsav*2.0;
    for (nstp=1;nstp<=MAXSTP;nstp++) {
	(*derivs)(x,y,dydx);
	for (i=1;i<=nvar;i++)
	    yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
	if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
	    xp[++kount]=x;
	    for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
	    xsav=x;
	}
	if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
	(*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
	if (hdid == h) ++(*nok); else ++(*nbad);
	if ((x-x2)*(x2-x1) >= 0.0) {
	    for (i=1;i<=nvar;i++) ystart[i]=y[i];
	    if (kmax) {
		xp[++kount]=x;
		for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
	    }
	    free_dvector(dydx,1,nvar);
	    free_dvector(y,1,nvar);
	    free_dvector(yscal,1,nvar);
	    return;
	}
	if (fabs(hnext) <= hmin) nrerror("Step size too small in odeint");
	h=hnext;
    }
    nrerror("Too many steps in routine odeint");
}
#undef MAXSTP
#undef TINY
#undef NRANSI


void rk4(double y[], double dydx[], int n, double x, double h, double yout[],
	void (*derivs)(double, double [], double []))
{
    int i;
    double xh,hh,h6,*dym,*dyt,*yt;
    
    dym=dvector(1,n);
    dyt=dvector(1,n);
    yt=dvector(1,n);
    hh=h*0.5;
    h6=h/6.0;
    xh=x+hh;
    for (i=1;i<=n;i++) yt[i]=y[i]+hh*dydx[i];
    (*derivs)(xh,yt,dyt);
    for (i=1;i<=n;i++) yt[i]=y[i]+hh*dyt[i];
    (*derivs)(xh,yt,dym);
    for (i=1;i<=n;i++) {
	yt[i]=y[i]+h*dym[i];
	dym[i] += dyt[i];
    }
    (*derivs)(x+h,yt,dyt);
    for (i=1;i<=n;i++)
	yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
    free_dvector(yt,1,n);
    free_dvector(dyt,1,n);
    free_dvector(dym,1,n);
}


void rkck(double y[], double dydx[], int n, double x, double h, double yout[],
	double yerr[], void (*derivs)(double, double [], double []))
{
    int i;
    static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
	b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
	b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
	b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
	b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
	c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
	dc5 = -277.00/14336.0;
    double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
	dc4=c4-13525.0/55296.0,dc6=c6-0.25;
    double *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;
    
    ak2=dvector(1,n);
    ak3=dvector(1,n);
    ak4=dvector(1,n);
    ak5=dvector(1,n);
    ak6=dvector(1,n);
    ytemp=dvector(1,n);
    for (i=1;i<=n;i++)
	ytemp[i]=y[i]+b21*h*dydx[i];
    (*derivs)(x+a2*h,ytemp,ak2);
    for (i=1;i<=n;i++)
	ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
    (*derivs)(x+a3*h,ytemp,ak3);
    for (i=1;i<=n;i++)
	ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
    (*derivs)(x+a4*h,ytemp,ak4);
    for (i=1;i<=n;i++)
	ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
    (*derivs)(x+a5*h,ytemp,ak5);
    for (i=1;i<=n;i++)
	ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
    (*derivs)(x+a6*h,ytemp,ak6);
    for (i=1;i<=n;i++)
	yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
    for (i=1;i<=n;i++)
	yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
    free_dvector(ytemp,1,n);
    free_dvector(ak6,1,n);
    free_dvector(ak5,1,n);
    free_dvector(ak4,1,n);
    free_dvector(ak3,1,n);
    free_dvector(ak2,1,n);
}


#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []))
{
    void rkck(double y[], double dydx[], int n, double x, double h,
	      double yout[], double yerr[], 
	      void (*derivs)(double, double [], double []));
    int i;
    double errmax,h,htemp,xnew,*yerr,*ytemp;
    
    yerr=dvector(1,n);
    ytemp=dvector(1,n);
    h=htry;
    for (;;) {	
	rkck(y,dydx,n,*x,h,ytemp,yerr,derivs);
	errmax=0.0;
	for (i=1;i<=n;i++) errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
	errmax /= eps;
	if (errmax <= 1.0) break;
	htemp=SAFETY*h*pow(errmax,PSHRNK);
	h=(h >= 0.0 ? FMAX(htemp,0.1*h) : FMIN(htemp,0.1*h));
	xnew=(*x)+h;
	if (xnew == *x) nrerror("stepsize underflow in rkqs");
    }
    if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
    else *hnext=5.0*h;
    *x += (*hdid=h);
    for (i=1;i<=n;i++) y[i]=ytemp[i];
    free_dvector(ytemp,1,n);
    free_dvector(yerr,1,n);
}

#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON


void input4(double v[])
{

  int k;
  float c0, c2, c4, phi3;

  printf("Write the initial values of c0,c2,c4\n");
  scanf("%f %f %f",&c0,&c2,&c4);
  printf("Write the initial values of  phi1,phi3\n");
  scanf("%f %f",&phi1,&phi3);
  printf("You wrote c0=%1.6f c2=%1.6f c4=%1.6f\n",c0,c2,c4);
  printf("phi1= %1.6f phi3= %1.6f\n",phi1,phi3); 
    
  v[1]=c0; v[2]=c2; v[3]=c4; 
  v[4]= phi3;
  
  for(k=1;k<=4;k++) printf("v[%d]= %e\n",k,v[k]);

}

void input6(double v[])
{
    float c0, c2, c4, c6, phi3, phi5;
    
    printf("Write the initial values of  c0,c2,c4,c6\n");
    scanf("%f %f %f %f",&c0,&c2,&c4,&c6);
    printf("Write the initial values of  phi1,phi3,phi5\n");
    scanf("%f %f %f",&phi1,&phi3,&phi5);
    printf("You wrote c0=%1.6f c2=%1.6f c4=%1.6f c6=%1.6f\n",
	   c0,c2,c4,c6);
    printf("phi1= %1.6f phi3= %1.6f phi5= %1.6f\n",phi1,phi3,phi5);
    
    v[1]=c0; v[2]=c2; v[3]=c4; v[4]= c6; 
    v[5]= phi3; v[6]=phi5;
}

void input8(double v[])
{
    float c0,c2,c4,c6,c8,phi3,phi5,phi7;

    printf("Write the initial values of  c0,c2,c4,c6,c8\n");
    scanf("%f %f %f %f %f",&c0,&c2,&c4,&c6,&c8);
    printf("Write the initial values of  phi1,phi3,phi5,phi7\n");
    scanf("%f %f %f %f",&phi1,&phi3,&phi5,&phi7);
    printf("You wrote c0=%1.6f c2=%1.6f c4=%1.6f c6=%1.6f c8=%1.6f\n",
	   c0,c2,c4,c6,c8);
    printf("phi1= %1.6f phi3= %1.6f phi5= %1.6f phi7=%1.6f\n",
	   phi1,phi3,phi5,phi7);
        
    v[1]=c0; v[2]=c2; v[3]=c4; v[4]= c6; v[5]=c8; 
    v[6]=phi3; v[7]= phi5; v[8]= phi7;
}

void input10(double v[])
{
    float c0,c2,c4,c6,c8,c10,phi3,phi5,phi7,phi9;
    
    printf("Write the initial values of  c0,c2,c4,c6,c8,c10\n");
    scanf("%f %f %f %f %f %f",&c0,&c2,&c4,&c6,&c8,&c10);
    printf("Write the initial values of  phi1,phi3,phi5,phi7,phi9\n"); 
    scanf("%f %f %f %f %f",&phi1,&phi3,&phi5,&phi7,&phi9);
    printf("You wrote c0=%1.6f c2=%1.6f c4=%1.6f c6=%1.6f c8=%1.6f c10=%1.6f\n",
	   c0,c2,c4,c6,c8,c10);
    printf("phi1= %1.6f phi3= %1.6f phi5= %1.6f phi7= %1.6f phi9= %1.6f\n",
	   phi1,phi3,phi5,phi7,phi9);
        
    v[1]=c0; v[2]=c2; v[3]=c4; v[4]= c6; v[5]=c8; v[6]= c10;
    v[7]=phi3; v[8]= phi5; v[9]= phi7; v[10]= phi9;

}

void input12(double v[])
{
    float c0,c2,c4,c6,c8,c10,c12,phi3,phi5,phi7,phi9,phi11;

    printf("Write the initial values of  c0,c2,c4,c6,c8,c10,c12\n");
    scanf("%f %f %f %f %f %f %f",&c0,&c2,&c4,&c6,&c8,&c10,&c12);
    printf("Write the initial values of  phi1,phi3,phi5,phi7,phi9\n"); 
    scanf("%f %f %f %f %f %f",&phi1,&phi3,&phi5,&phi7,&phi9,&phi11);
    printf("You wrote c0=%1.6f c2=%1.6f c4=%1.6f c6=%1.6f\n",
	   c0,c2,c4,c6); 
    printf("c8=%1.6f c10=%1.6f c12=%1.6f\n",c8,c10,c12);	   
    printf("phi1= %1.6f phi3= %1.6f phi5= %1.6f phi7= %1.6f\n",
	   phi1,phi3,phi5,phi7);
    printf("phi9= %1.6f phi11= %1.6f\n",phi9,phi11);
        
    v[1]=c0; v[2]=c2; v[3]=c4; v[4]= c6; v[5]=c8; v[6]= c10; v[7]= c12;
    v[8]=phi3; v[9]= phi5; v[10]= phi7; v[11]= phi9; v[12]= phi11;
}

void printvs(double x2, double v[])
{
    if(nnn==2) {
	printf("x=%f c0=%e c2=%e\n",x2,v[1],v[2]);
    }
    if(nnn==4) {
	printf("x=%f c0=%e %e %e\n",x2,v[1],v[2],v[3]);
	printf("%e\n",v[4]);
    }
    if(nnn==6) {
      printf("x=%f %e %e %e %e\n",x2,v[1],v[2],v[3],v[4]);
	printf("phi3=%e %e\n",v[5],v[6]);
    }
    if(nnn==8) {
	printf("x=%f %e %e %e %e %e\n",
	       x2,v[1],v[2],v[3],v[4],v[5]);
	printf("phi3=%e %e %e\n",v[6],v[7],v[8]);
    }
    if(nnn==10) {
	printf("x=%f c0=%e c2=%e c4=%e c6=%e\n",
	       x2,v[1],v[2],v[3],v[4]);
	printf("c8=%e c10=%e phi3=%e phi5=%e\n",
	       v[5],v[6],v[7],v[8]);
	printf("phi7=%e phi9=%e\n",v[9],v[10]);
    }
    if(nnn==12) {
	printf("x=%f c0=%e c2=%e c4=%e c6=%e\n",
	       x2,v[1],v[2],v[3],v[4]);
	printf("c8=%e c10=%e c12=%e phi3=%e\n",
	       v[5],v[6],v[7],v[8]);
	printf("phi5=%e phi7=%e phi9=%e phi11=%e\n",
	       v[9],v[10],v[11],v[12]);
    }

}
/*** Here A[j] j=1,2,..,N are computed **/

void matriz(double x, double yout[])
{
  extern int nvar;
  extern double A[N+1];
  int n,j,k,p,la,*indx, jmax;
  double s1,s2,s3,s4;
  double **aa,*b,d;
  double a[N+1],c[N+1],phi[N+1],xi[N+1];
  void ludcmp(double **aa, int n, int *indx, double *d);
  void lubksb(double **aa,int n,int *indx,double b[]);
  
  jmax= nnn/2; 
  indx=ivector(1,jmax);
  b=dvector(1,jmax);
  aa=dmatrix(1,jmax,1,jmax);
  
  if( fabs(x) <= 1.0e-8 ) {
    for(j=1; j<= N; j++) A[j]=0.0;
    A[0] = 1.0;
  }
  else {
    A[0]=a[0]=yout[1];  
    for(j=0, k=2; j <= jmax; j++, k +=1) c[2*j]= yout[k];
    for(j=1;j <= jmax; j++, k +=1) phi[2*j-1]= yout[k];
    for(j=1;j <= jmax; j++, k +=1) xi[2*j-1]= yout[k];
    
    /** Computing matrix C **/
    for(k=1;k<= jmax;k++){
      for(n=1;n<=jmax;n++){
	for(p=1,s1=0.0,s2=0.0,s3=0.0,s4=0.0;p<=jmax;p++){
	  for(la=1;la<=jmax;la++){
	    if(2*p-1-abs(2*n-2*la+1)== 2*k) 
	      s1 += (2*p-1)*phi[2*p-1]*xi[2*la-1];
	    if(2*p-1-abs(2*n-2*la+1)== -2*k) 
	      s1 += -(2*p-1)*phi[2*p-1]*xi[2*la-1];  	    
	    if(2*p-1+abs(2*n-2*la+1)== 2*k) 
	      s2 += (2*p-1)*phi[2*p-1]*xi[2*la-1];    	    
	    if(p-(n+la)==k) 
	      s3 += (2*p-1)*phi[2*p-1]*xi[2*la-1];  	     
	    if(p-(n+la)== -k) 
	      s3 += -(2*p-1)*phi[2*p-1]*xi[2*la-1];	     
	    if(p+n+la-1==k) 
	      s4 += (2*p-1)*phi[2*p-1]*xi[2*la-1];	     
	  }
	}
	aa[k][n]= 0.25*x*(s1+s2+s3+s4);
      }
    }
    
    for(n=1;n<=jmax;n++) aa[n][n] += -2.0*n;
    
    for(k=1;k<=jmax;k++) {
      for(p=1,s1=0.0,s2=0.0;p<=jmax;p++){
	for(la=1;la<=jmax;la++){
	  if(p+la-1==k) s1 += (2*p-1)*phi[2*p-1]*xi[2*la-1];
	  if(p-la==k) s2 += (2*p-1)*phi[2*p-1]*xi[2*la-1];
	  if(p-la== -k) s2 += - (2*p-1)*phi[2*p-1]*xi[2*la-1];
	}
      }
      b[k]= -0.5*x*(s1+s2)*A[0]; 
    }
    
    ludcmp(aa,jmax,indx,&d);
    lubksb(aa,jmax,indx,b);
    
    for(k=1;k<=jmax;k++) A[2*k]= b[k];
  }
  
  free_dmatrix(aa,1,jmax,1,jmax);
  free_dvector(b,1,jmax);
  free_ivector(indx,1,jmax);
}

// Printing
void imprimir(double x1, double y[])
{
    extern int nnn;
    extern double A[N+1], wt1, wt2;
    int k,j,jmax;
    double a[N+1],c[N+1],phi[N+1],xi[N+1];
    double rho1,rho2, sumaA, sumaB, sumaC,sumaPHI,sumaXI;
    void densidad(double *rho,double a[],double c[], double phi[],
		  double xi[],double wt);
    
    jmax = nnn/2;

    for(j=0;j<=N;j++) {a[j]=c[j]=phi[j]=xi[j]=0.0;}
    for(j=1;j<=jmax;j++) a[2*j]=A[2*j];
    a[0]=y[1];
    for(j=0,k=2;j <= jmax; j++, k +=1) c[2*j]= y[k];
    for(j=1;j <= jmax; j++, k +=1) phi[2*j-1]= y[k];
    for(j=1;j <= jmax; j++, k +=1) xi[2*j-1]= y[k];
    

    for(sumaA=a[0],j=1;j<=jmax;j++) sumaA += a[2*j]; 
    for(sumaC=0.0,j=0;j<=jmax;j++) sumaC += c[2*j];
    for(sumaPHI=0.0,j=1;j<=jmax;j++) sumaPHI += phi[2*j-1];
    for(sumaXI=0.0,j=1;j<=jmax;j++) sumaXI += xi[2*j-1];


    sumaB= sumaA/sumaC;

    if(nnn==2){
	if(fabs(x1)<=1.0e-15){
	    densidad(&rho1,a,c,phi,xi,wt1);
	    densidad(&rho2,a,c,phi,xi,wt2);
	    fprintf(fp5,"%e\t%e\t%e\n",x1,rho1,rho2);
	    fprintf(fp1,"%e\t%e\t%e\n",x1,y[1],0.0);
	    fprintf(fp2,"%e\t%e\t%e\n",x1,c[0],c[2]);
	    fprintf(fp3,"%e\t%e\n",x1,phi[1]);
	    fprintf(fp4,"%e\t%e\n",x1,xi[1]);
	}
	else {
	    densidad(&rho1,a,c,phi,xi,wt1);
	    densidad(&rho2,a,c,phi,xi,wt2);
	    fprintf(fp5,"%e\t%e\t%e\n",x1,rho1,rho2);
	    fprintf(fp1,"%e\t%e\t%e\n",x1,y[1],A[2]);
	    fprintf(fp2,"%e\t%e\t%e\n",x1,c[0],c[2]);
	    fprintf(fp3,"%e\t%e\n",x1,phi[1]);
	    fprintf(fp4,"%e\t%e\n",x1,xi[1]);
	}
    }

    if(nnn==4){
	if(fabs(x1)<=1.0e-15){
	    densidad(&rho1,a,c,phi,xi,wt1);
	    densidad(&rho2,a,c,phi,xi,wt2);
	    fprintf(fp5,"%e\t%e\t%e\n",x1,rho1,rho2);
	    fprintf(fp1,"%e\t%e\t%e\t%e\n",x1,y[1],0.0,0.0);
	    fprintf(fp2,"%e\t%e\t%e\t%e\n",x1,c[0],c[2],c[4]);
	    fprintf(fp3,"%e\t%e\t%e\n",x1,phi[1],phi[3]);
	    fprintf(fp4,"%e\t%e\t%e\n",x1,xi[1],xi[3]);
	    fprintf(fpA,"%e\t%e\n",x1,sumaA);
	    fprintf(fpC,"%e\t%e\n",x1,sumaC);
	    fprintf(fpB,"%e\t%e\n",x1,sumaB);
	    fprintf(fpP,"%e\t%e\n",x1,sumaPHI);
	    fprintf(fpX,"%e\t%e\n",x1,sumaXI);
	}
	else{
	    densidad(&rho1,a,c,phi,xi,wt1);
	    densidad(&rho2,a,c,phi,xi,wt2);
	    fprintf(fp5,"%e\t%e\t%e\n",x1,rho1,rho2);
	    fprintf(fp1,"%e\t%e\t%e\t%e\n",x1,y[1],A[2],A[4]);
	    fprintf(fp2,"%e\t%e\t%e\t%e\n",x1,c[0],c[2],c[4]);
	    fprintf(fp3,"%e\t%e\t%e\n",x1,phi[1],phi[3]);
	    fprintf(fp4,"%e\t%e\t%e\n",x1,xi[1],xi[3]);
	    fprintf(fpA,"%e\t%e\n",x1,sumaA);
	    fprintf(fpB,"%e\t%e\n",x1,sumaB);
	    fprintf(fpC,"%e\t%e\n",x1,sumaC);
	    fprintf(fpP,"%e\t%e\n",x1,sumaPHI);
	    fprintf(fpX,"%e\t%e\n",x1,sumaXI);
	}
    }

    if(nnn==6){
	if(fabs(x1)<=1.0e-15){
	    densidad(&rho1,a,c,phi,xi,wt1);
	    densidad(&rho2,a,c,phi,xi,wt2);
	    fprintf(fp5,"%e\t%e\t%e\n",x1,rho1,rho2);
	    fprintf(fp1,"%e\t%e\t%e\t%e\t%e\n",x1,y[1],0.0,0.0,0.0);
	    fprintf(fp2,"%e\t%e\t%e\t%e\t%e\n",x1,c[0],c[2],c[4],c[6]);
	    fprintf(fp3,"%e\t%e\t%e\t%e\n",x1,phi[1],phi[3],phi[5]);
	    fprintf(fp4,"%e\t%e\t%e\t%e\n",x1,xi[1],xi[3],xi[5]);
	}
	else {
	    densidad(&rho1,a,c,phi,xi,wt1);
	    densidad(&rho2,a,c,phi,xi,wt2);
	    fprintf(fp5,"%e\t%e\t%e\n",x1,rho1,rho2);
	    fprintf(fp1,"%e\t%e\t%e\t%e\t%e\n",x1,y[1],A[2],A[4],A[6]);
	    fprintf(fp2,"%e\t%e\t%e\t%e\t%e\n",x1,c[0],c[2],c[4],c[6]);
	    fprintf(fp3,"%e\t%e\t%e\t%e\n",x1,phi[1],phi[3],phi[5]);
	    fprintf(fp4,"%e\t%e\t%e\t%e\n",x1,xi[1],xi[3],xi[5]);
	}
    }

    if(nnn==8){
	if(fabs(x1)<=1.0e-15){
	    densidad(&rho1,a,c,phi,xi,wt1);
	    densidad(&rho2,a,c,phi,xi,wt2);
	    fprintf(fp5,"%e\t%e\t%e\n",x1,rho1,rho2);
	    fprintf(fp1,"%e\t%e\t%e\t%e\t%e\t%e\n",x1,y[1],0.0,0.0,0.0,0.0);
	    fprintf(fp2,"%e\t%e\t%e\t%e\t%e\t%e\n",x1,c[0],c[2],c[4],c[6],c[8]);
	    fprintf(fp3,"%e\t%e\t%e\t%e\t%e\n",x1,phi[1],phi[3],phi[5],phi[7]);
	    fprintf(fp4,"%e\t%e\t%e\t%e\t%e\n",x1,xi[1],xi[3],xi[5],xi[7]);
	}
	else {
	    densidad(&rho1,a,c,phi,xi,wt1);
	    densidad(&rho2,a,c,phi,xi,wt2);
	    fprintf(fp5,"%e\t%e\t%e\n",x1,rho1,rho2);
	    fprintf(fp1,"%e\t%e\t%e\t%e\t%e\t%e\n",x1,y[1],A[2],A[4],A[6],A[8]);
	    fprintf(fp2,"%e\t%e\t%e\t%e\t%e\t%e\n",x1,c[0],c[2],c[4],c[6],c[8]);
	    fprintf(fp3,"%e\t%e\t%e\t%e\t%e\n",x1,phi[1],phi[3],phi[5],phi[7]);
	    fprintf(fp4,"%e\t%e\t%e\t%e\t%e\n",x1,xi[1],xi[3],xi[5],xi[7]);
	}
    }

    if(nnn==10){
	if(fabs(x1)<=1.0e-15){
	densidad(&rho1,a,c,phi,xi,wt1);
	densidad(&rho2,a,c,phi,xi,wt2);
	fprintf(fp5,"%e\t%e\t%e\n",x1,rho1,rho2);
	fprintf(fp1,"%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
		x1,y[1],0.0,0.0,0.0,0.0,0.0);
	fprintf(fp2,"%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
		x1,c[0],c[2],c[4],c[6],c[8],c[10]);
	fprintf(fp3,"%e\t%e\t%e\t%e\t%e\t%e\n",
		x1,phi[1],phi[3],phi[5],phi[7],phi[9]);
	fprintf(fp4,"%e\t%e\t%e\t%e\t%e\t%e\n",
		x1,xi[1],xi[3],xi[5],xi[7],xi[9]);
	}
	else {
	    densidad(&rho1,a,c,phi,xi,wt1);
	    densidad(&rho2,a,c,phi,xi,wt2);
	    fprintf(fp5,"%e\t%e\t%e\n",x1,rho1,rho2);
	    fprintf(fp1,"%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
		    x1,y[1],A[2],A[4],A[6],A[8],A[10]);
	    fprintf(fp2,"%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
		    x1,c[0],c[2],c[4],c[6],c[8],c[10]);
	    fprintf(fp3,"%e\t%e\t%e\t%e\t%e\t%e\n",
		    x1,phi[1],phi[3],phi[5],phi[7],phi[9]);
	    fprintf(fp4,"%e\t%e\t%e\t%e\t%e\t%e\n",
		    x1,xi[1],xi[3],xi[5],xi[7],xi[9]);
	}
    }

    if(nnn==12){
	if(fabs(x1)<=1.0e-15){
	    densidad(&rho1,a,c,phi,xi,wt1);
	    densidad(&rho2,a,c,phi,xi,wt2);
	    fprintf(fp5,"%e\t%e\t%e\n",x1,rho1,rho2);
	    fprintf(fp1,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
		    x1,y[1],0.0,0.0,0.0,0.0,0.0,0.0);
	    fprintf(fp2,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
		    x1,c[0],c[2],c[4],c[6],c[8],c[10],c[12]);
	    fprintf(fp3,"%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
		    x1,phi[1],phi[3],phi[5],phi[7],phi[9],phi[11]);
	    fprintf(fp4,"%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
		    x1,xi[1],xi[3],xi[5],xi[7],xi[9],xi[11]);
	}
	else {
	    densidad(&rho1,a,c,phi,xi,wt1);
	    densidad(&rho2,a,c,phi,xi,wt2);
	    fprintf(fp5,"%e\t%e\t%e\n",x1,rho1,rho2);
	    fprintf(fp1,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
		    x1,y[1],A[2],A[4],A[6],A[8],A[10],A[12]);
	    fprintf(fp2,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
		    x1,c[0],c[2],c[4],c[6],c[8],c[10],c[12]);
	    fprintf(fp3,"%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
		    x1,phi[1],phi[3],phi[5],phi[7],phi[9],phi[11]);
	    fprintf(fp4,"%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
		    x1,xi[1],xi[3],xi[5],xi[7],xi[9],xi[11]);   
	}
    }
}

void densidad(double *rho, double a[], double c[], double phi[],
	      double xi[], double wt)
{
    extern int nnn;
    int k;
    double uno,dos,tres,sum_phi2,sum_phiprim2;
    double sum_c,sum_a;

    for(sum_c=0.0,sum_a=0.0,k=0;k<=nnn;k +=2){
	sum_a += a[k]*cos(k*wt); 
	sum_c += c[k]*cos(k*wt);
    }
    
    if(sum_a<1.0e-15){printf("A^-1(x,t) goes to infinity, exiting\n"); 
	exit(1);} 
    
    for(sum_phi2=0.0,sum_phiprim2=0.0,k=1;k<=nnn;k += 2){
	sum_phi2 += phi[k]*k*sin(k*wt);
	sum_phiprim2 += xi[k]*cos(k*wt);
    }
    
    sum_phi2 = sum_phi2*sum_phi2;
    sum_phiprim2 = sum_phiprim2*sum_phiprim2;
    
    uno= sum_c*sum_phi2/sum_a;
    dos= sum_phiprim2/sum_a;
    
    for(tres=0.0,k=1;k<=nnn;k += 2)
	tres += phi[k]*cos(k*wt);
    
    *rho= 0.5*(uno+dos+tres*tres+0.5*GAMA*tres*tres*tres*tres);
}






