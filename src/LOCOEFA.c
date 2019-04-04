#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define N 50
#define DEBUG 0
#define MAXERROR 1e-13

double T;
int K;
double area;

typedef struct datapoints {
  double x;
  double y;
  double t;
  double deltax;
  double deltay;
  double deltat;
  double sumdeltaxj;
  double sumdeltayj;
  double xi;
  double epsilon;
} Datapoints;

Datapoints *contour;

typedef struct modes {
  double alpha;
  double beta;
  double gamma;
  double delta;
  double tau;
  double alphaprime;
  double gammaprime;
  double rho;
  double alphastar;
  double betastar;
  double gammastar;
  double deltastar;
  double r;
  double a;
  double b;
  double c;
  double d;
  double aprime;
  double bprime;
  double cprime;
  double dprime;
  double phi;
  double theta;
  double lambda1;
  double lambda2;
  double lambda21;
  double lambda12;
  double lambdaplus;
  double lambdaminus;
  double zetaplus;
  double zetaminus;
  double locooffseta;
  double locooffsetc;
  double locolambdaplus;
  double locolambdaminus;
  double locozetaplus;
  double locozetaminus;
  double locoL;
  double locoaplus;
  double locobplus;
  double lococplus;
  double locodplus;
  double locoaminus;
  double locobminus;
  double lococminus;
  double locodminus;
  } Modes;

Modes mode[N+2];

void ReadData(double *x,int *n,double *out)
{
  int i;
  int lines;

  lines=n[0]+1;
  //printf("datapoints: %d\n",lines);
  
  if(((contour=(Datapoints *)calloc(lines,sizeof(Datapoints)))==NULL)) {
    fprintf(stderr,"Error in memory allocation\n");
    exit(EXIT_FAILURE);
  }
  
  K=0;
  area=x[3*(lines-1)];
  for(i=0;i<lines-1;i++) {
    contour[K].x=x[i+lines-1];
    contour[K].y=x[i+2*(lines-1)];
    //check if not accidentally same point; if so, ignore
    if(K==0 || fabs(contour[K].x-contour[K-1].x)>MAXERROR ||
       fabs(contour[K].y-contour[K-1].y)>MAXERROR)
      K++;
  }
  if(K && (fabs(contour[0].x-contour[K-1].x)>MAXERROR ||
	   fabs(contour[0].y-contour[K-1].y)>MAXERROR)) {
    //if last point not equal to first point, add one point equal to first point
    contour[K].x=contour[0].x;
    contour[K].y=contour[0].y;
  }
  else
    K--;

  if(K<2) {
    fprintf(stderr,"Error: contour contains only one point\n");
    exit(EXIT_FAILURE);
  }
  
  if(DEBUG) {
    fprintf(stdout,"\n\ndata points:\n============\n\n");
    for(i=0;i<=K;i++)
      fprintf(stdout,"point %d:\tx=%lf\ty=%lf\n",i,contour[i].x,contour[i].y);
  }
}

void CalculateEFACoefficients()
{
  int i,j;

  //Below eq.5: Deltax, DeltaY, DeltaT
  contour[0].t=0.;
  for(i=1;i<=K;i++) {
    contour[i].deltax=contour[i].x-contour[i-1].x;
    contour[i].deltay=contour[i].y-contour[i-1].y;
    contour[i].deltat=hypot(contour[i].deltax,contour[i].deltay);
    contour[i].t=contour[i-1].t+contour[i].deltat;
  }
  
  //Below eq.5: T
  T=contour[K].t;

  if(T<MAXERROR) {
    fprintf(stderr,"Error: contour of length zero\n");
    exit(EXIT_FAILURE);
  }

  //Below eq. 7: sumDeltaxj, sumDeltayj, xi, epsilon
  contour[1].xi=0.;
  contour[1].epsilon=0.;
  contour[1].sumdeltaxj=0.;
  contour[1].sumdeltayj=0.;
  for(i=2;i<=K;i++) {
    contour[i].sumdeltaxj=contour[i-1].sumdeltaxj+contour[i-1].deltax;
    contour[i].sumdeltayj=contour[i-1].sumdeltayj+contour[i-1].deltay;
    contour[i].xi=contour[i].sumdeltaxj-contour[i].deltax/contour[i].deltat*contour[i-1].t;
    contour[i].epsilon=contour[i].sumdeltayj-contour[i].deltay/contour[i].deltat*contour[i-1].t;
  }

  //Equation 7: alpha0, gamma0
  for(j=0;j<=N;j++) {
    mode[j].alpha=0.;
    mode[j].beta=0.;
    mode[j].gamma=0.;
    mode[j].delta=0.;
  }
  mode[0].alpha=contour[0].x;
  mode[0].gamma=contour[0].y;
  for(i=1;i<=K;i++) {
    mode[0].alpha+=(contour[i].deltax/(2.*contour[i].deltat)*(contour[i].t*contour[i].t-contour[i-1].t*contour[i-1].t)+contour[i].xi*(contour[i].t-contour[i-1].t))/T;
    mode[0].gamma+=(contour[i].deltay/(2.*contour[i].deltat)*(contour[i].t*contour[i].t-contour[i-1].t*contour[i-1].t)+contour[i].epsilon*(contour[i].t-contour[i-1].t))/T;

    //Equation 6: alpha, beta, gamma, delta
    for(j=1;j<=N;j++) {
      mode[j].alpha+=contour[i].deltax/contour[i].deltat*(cos(2.*(double)j*M_PI*contour[i].t/T)-cos(2*(double)j*M_PI*contour[i-1].t/T));
      mode[j].beta+=contour[i].deltax/contour[i].deltat*(sin(2.*(double)j*M_PI*contour[i].t/T)-sin(2*(double)j*M_PI*contour[i-1].t/T));
      mode[j].gamma+=contour[i].deltay/contour[i].deltat*(cos(2.*(double)j*M_PI*contour[i].t/T)-cos(2*(double)j*M_PI*contour[i-1].t/T));
      mode[j].delta+=contour[i].deltay/contour[i].deltat*(sin(2.*(double)j*M_PI*contour[i].t/T)-sin(2*(double)j*M_PI*contour[i-1].t/T));
    }
  }
  for(j=1;j<=N;j++) {
    mode[j].alpha*=T/(2.*(double)(j*j)*M_PI*M_PI);
    mode[j].beta*=T/(2.*(double)(j*j)*M_PI*M_PI);
    mode[j].gamma*=T/(2.*(double)(j*j)*M_PI*M_PI);
    mode[j].delta*=T/(2.*(double)(j*j)*M_PI*M_PI);
  }

  if(DEBUG) {
    fprintf(stdout,"\n\nEFA coefficients:\n=================\n\n");
    for(j=0;j<=N;j++) {
      fprintf(stdout,"mode %d:\n",j);
      fprintf(stdout,"(%g\t%g)\n",mode[j].alpha,mode[j].beta);
      fprintf(stdout,"(%g\t%g)\n\n",mode[j].gamma,mode[j].delta);
    }
  }
}    

void CalculateLOCOCoefficients(double *out)
{
  int i;

  //Equation 14: tau1
  mode[1].tau=0.5*atan2(2.*(mode[1].alpha*mode[1].beta+mode[1].gamma*mode[1].delta),mode[1].alpha*mode[1].alpha+mode[1].gamma*mode[1].gamma-mode[1].beta*mode[1].beta-mode[1].delta*mode[1].delta); 

  //Below eq. 15: alpha1prime, gamma1prime
  mode[1].alphaprime=mode[1].alpha*cos(mode[1].tau)+mode[1].beta*sin(mode[1].tau);
  mode[1].gammaprime=mode[1].gamma*cos(mode[1].tau)+mode[1].delta*sin(mode[1].tau);

  //Equation 16: rho
  mode[1].rho=atan2(mode[1].gammaprime,mode[1].alphaprime);

  //Equation 17: tau1
  if(mode[1].rho<0.)
    mode[1].tau+=M_PI;

  //Equation 18: alphastar, betastar, gammastar, deltastar
  for(i=1;i<=N;i++) {
    mode[i].alphastar=mode[i].alpha*cos((double)i*mode[1].tau)+mode[i].beta*sin((double)i*mode[1].tau);
    mode[i].betastar=-mode[i].alpha*sin((double)i*mode[1].tau)+mode[i].beta*cos((double)i*mode[1].tau);
    mode[i].gammastar=mode[i].gamma*cos((double)i*mode[1].tau)+mode[i].delta*sin((double)i*mode[1].tau);
    mode[i].deltastar=-mode[i].gamma*sin((double)i*mode[1].tau)+mode[i].delta*cos((double)i*mode[1].tau);
  }

  //Equation 9: r
  mode[1].r=mode[1].alphastar*mode[1].deltastar-mode[1].betastar*mode[1].gammastar;

  //Equation 19: betastar, deltastar
  if(mode[1].r<0.) {
    for(i=1;i<=N;i++) {
      mode[i].betastar=-mode[i].betastar;
      mode[i].deltastar=-mode[i].deltastar;
    }
  }

  //Equation 20: a, b, c, d
  mode[0].a=mode[0].alpha;
  mode[0].c=mode[0].gamma;

  for(i=1;i<=N;i++) {
    mode[i].a=mode[i].alphastar;
    mode[i].b=mode[i].betastar;
    mode[i].c=mode[i].gammastar;
    mode[i].d=mode[i].deltastar;
    }

  if(DEBUG) {
    fprintf(stdout,"\n\nmodified EFA coefficients:\n==========================\n\n");
    for(i=0;i<=N;i++) {
      fprintf(stdout,"mode %d:\n",i);
      fprintf(stdout,"(%g\t%g)\n",mode[i].a,mode[i].b);
      fprintf(stdout,"(%g\t%g)\n\n",mode[i].c,mode[i].d);
    }
  }
  
  if(DEBUG)
    fprintf(stdout,"\n\nLambda matrices:\n================\n\n");

  for(i=1;i<=N;i++) {
    //Equation 26: phi
    mode[i].phi=0.5*atan2(2.*(mode[i].a*mode[i].b+mode[i].c*mode[i].d),mode[i].a*mode[i].a+mode[i].c*mode[i].c-mode[i].b*mode[i].b-mode[i].d*mode[i].d); 

    //Below eq. 27: aprime, bprime, cprime, dprime
    mode[i].aprime=mode[i].a*cos(mode[i].phi)+mode[i].b*sin(mode[i].phi);
    mode[i].bprime=-mode[i].a*sin(mode[i].phi)+mode[i].b*cos(mode[i].phi);
    mode[i].cprime=mode[i].c*cos(mode[i].phi)+mode[i].d*sin(mode[i].phi);
    mode[i].dprime=-mode[i].c*sin(mode[i].phi)+mode[i].d*cos(mode[i].phi);

    //Equation 27: theta
    mode[i].theta=atan2(mode[i].cprime,mode[i].aprime);

    //Equation 25: Lambda
    mode[i].lambda1=cos(mode[i].theta)*mode[i].aprime+sin(mode[i].theta)*mode[i].cprime;
    mode[i].lambda12=cos(mode[i].theta)*mode[i].bprime+sin(mode[i].theta)*mode[i].dprime;
    mode[i].lambda21=-sin(mode[i].theta)*mode[i].aprime+cos(mode[i].theta)*mode[i].cprime;
    mode[i].lambda2=-sin(mode[i].theta)*mode[i].bprime+cos(mode[i].theta)*mode[i].dprime;

    if(DEBUG /*|| fabs(mode[i].lambda12)>MAXERROR || fabs(mode[i].lambda21)>MAXERROR || mode[i].lambda1<0 || mode[i].lambda1 < fabs(mode[i].lambda2)*/) {
      if(fabs(mode[i].lambda12)>MAXERROR || fabs(mode[i].lambda21)>MAXERROR)
	fprintf(stderr,"Warning: off-diagonal Lambda matrix unequal to zero:\n");
      if(mode[i].lambda1<0)
	fprintf(stderr,"Warning: lambda1 negative:\n");
      if(mode[i].lambda1 < fabs(mode[i].lambda2))
	fprintf(stderr,"Warning: lambda1 < |lambda2|:\n");
      fprintf(stdout,"mode %d:\n",i);
      fprintf(stdout,"(%g\t%g)\n",mode[i].lambda1,mode[i].lambda12);
      fprintf(stdout,"(%g\t%g)\n\n",mode[i].lambda21,mode[i].lambda2);
    }

    //Equation 32: lambdaplus, lambdaminus 
    mode[i].lambdaplus=(mode[i].lambda1+mode[i].lambda2)/2.;
    mode[i].lambdaminus=(mode[i].lambda1-mode[i].lambda2)/2.;

    //Below eq. 37: zetaplus, zetaminus
    mode[i].zetaplus=mode[i].theta-mode[i].phi;
    mode[i].zetaminus=-mode[i].theta-mode[i].phi;
  }

  //Below eq. 39: A0
  mode[0].locooffseta=mode[0].a;
  mode[0].locooffsetc=mode[0].c;

  if(DEBUG) {
    fprintf(stdout,"\n\noffset:\n===============\n\n");
    fprintf(stdout,"LOCO-EFA A0 offset:\ta=%g\tc=%g\n",mode[0].locooffseta,mode[0].locooffsetc);
  }

  //Below eq. 41: A+(l=0)
  mode[0].locolambdaplus=mode[2].lambdaplus;
  mode[0].locozetaplus=mode[2].zetaplus;

  //Below eq. 41: A+(l=1)
  mode[1].locolambdaplus=mode[1].lambdaplus;
  mode[1].locozetaplus=mode[1].zetaplus;

  //Below eq. 41: A+(l>1)
  for(i=2;i<=N-1;i++) {
    mode[i].locolambdaplus=mode[i+1].lambdaplus;
    mode[i].locozetaplus=mode[i+1].zetaplus;
  }

  //Below eq. 41: A-(l>0)
  for(i=2;i<=N+1;i++) {
    mode[i].locolambdaminus=mode[i-1].lambdaminus;
    mode[i].locozetaminus=mode[i-1].zetaminus;
  }

  if(1) {
    if(DEBUG)
      fprintf(stdout,"\n\nLn quadruplets:\n===============\n\n");
    for(i=0;i<=N+1;i++) {
      if(DEBUG) {
	fprintf(stdout,"LOCO-EFA mode %d:\tlambdaplus=%g\tlambdaminus=%g\tzetaplus=%g\tzetaminus=%g\n",i,mode[i].locolambdaplus,mode[i].locolambdaminus,mode[i].locozetaplus,mode[i].locozetaminus);
      }
      out[i]=i;
      out[i+2*(N+2)]=mode[i].locolambdaplus;
      out[i+3*(N+2)]=mode[i].locolambdaminus;
      out[i+4*(N+2)]=mode[i].locozetaplus;
      out[i+5*(N+2)]=mode[i].locozetaminus;
      out[i+6*(N+2)]=mode[i].locolambdaplus/(mode[i].locolambdaplus+mode[i].locolambdaminus+DBL_EPSILON);
    }
  }

  //Equation 38: Lambda*Zeta
  for(i=0;i<=N+1;i++) {
    mode[i].locoaplus=mode[i].locolambdaplus*cos(mode[i].locozetaplus);
    mode[i].locobplus=-mode[i].locolambdaplus*sin(mode[i].locozetaplus);
    mode[i].lococplus=mode[i].locolambdaplus*sin(mode[i].locozetaplus);
    mode[i].locodplus=mode[i].locolambdaplus*cos(mode[i].locozetaplus);
    mode[i].locoaminus=mode[i].locolambdaminus*cos(mode[i].locozetaminus);
    mode[i].locobminus=-mode[i].locolambdaminus*sin(mode[i].locozetaminus);
    mode[i].lococminus=-mode[i].locolambdaminus*sin(mode[i].locozetaminus);
    mode[i].locodminus=-mode[i].locolambdaminus*cos(mode[i].locozetaminus);
  }

  if(DEBUG) {
    fprintf(stdout,"\n\nLOCO coefficients:\n==================\n\n");
    for(i=0;i<=N+1;i++) {
      fprintf(stdout,"mode %d, Aplus:\n",i);
      fprintf(stdout,"(%g\t%g)\n",mode[i].locoaplus,mode[i].locobplus);
      fprintf(stdout,"(%g\t%g)\n",mode[i].lococplus,mode[i].locodplus);
      fprintf(stdout,"mode %d, Aminus:\n",i);
      fprintf(stdout,"(%g\t%g)\n",mode[i].locoaminus,mode[i].locobminus);
      fprintf(stdout,"(%g\t%g)\n",mode[i].lococminus,mode[i].locodminus);
    }
  }

  //Equation 47: L
  for(i=1;i<=N+1;i++) {
    mode[i].locoL=sqrt(mode[i].locolambdaplus*mode[i].locolambdaplus+mode[i].locolambdaminus*mode[i].locolambdaminus+2.*mode[i].locolambdaplus*mode[i].locolambdaminus*cos(mode[i].locozetaplus-mode[i].locozetaminus-2.*mode[1].locozetaplus));
    //the next line applies area normalisation, if wanted
    mode[i].locoL/=sqrt(area);
  }
  
  if(1) {
    if(DEBUG)
      fprintf(stdout,"\n\nLn scalar:\n==========\n\n");
    for(i=0;i<=N+1;i++) {
      if(DEBUG)
	fprintf(stdout,"LOCO-EFA mode %d:\tLn=%g\n",i,mode[i].locoL);
      out[i+N+2]=mode[i].locoL;
    }
  }
}

int LOCOEFA_(double *x,int *n,double *out)
{
  ReadData(x,n,out);
  CalculateEFACoefficients();
  CalculateLOCOCoefficients(out);
  return 0;
}  
