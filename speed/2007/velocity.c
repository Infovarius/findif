//----------------------- Calculation of known velocity profile  ----------//
#define LEVEL extern
#include "head.h"
#include <math.h>
#define Vmax 0.849931

double ch(double x)
{
 return((exp(x)+exp(-x))/2);
}

double th(double x)
{
 return( (exp(x)-exp(-x))/(exp(x)+exp(-x)) );
}

double fv2(double tcur,double a,double amp_beg,double amp_end,double tmax)
{
 return( amp_beg+ (amp_end-amp_beg)*(th(a*(tcur-tmax))+1)/2 );
}

double vmod(double t)
{
double a1=30,amp1=0,amp2=1,t1=0.03,vlim=0.1,dec=3,tbr=0.2,tpereg=0.3;
 if(TDV) return fv2(t, 10., fv2(t,a1,amp1,amp2,t1) - fv2(0,a1,amp1,amp2,t1), fv2(t,dec,1,vlim,tpereg), tbr);
    else return Vmax;
}

double chimod(double t,double phi)
{
 if(TDV) return( (th(-2*exp(-50*t)*(phi-40*M_PI*t))+1) / (th(80*exp(-50*t)*M_PI*t)+1) );
    else return chi;
}

double Gartmann(double rho,double rho0,double ksi)
{
 return( (ch(ksi)-ch(ksi*rho/rho0))/(ch(ksi)-1) );
}

double vfi_given(double t,double rho,double rho0)
{
const double ksi=Hksi;
 if(ksi>0) return( vmod(t)*Gartmann(rho,rho0,ksi)/Vmax);
 if(ksi==0) return( vmod(t)*Gartmann(rho,rho0,0.0001)/Vmax);
 return 1.;
}

double vtheta_given(double t,double rho,double rho0,double phi)
{
return( rho*chimod(t,phi)*vfi_given(t,rho,rho0)/rho0 );
}

void divertor(double t_cur)
{
  /*============================ divertor =================================*/
  // divertors are off till t sec
/*  if(t_cur<00)
  if(n[2]==0)
  for(i=0;i<m1;i++)
    for(j=ghost;j<=ghost;j++)
      for(k=0;k<m3;k++)
         {
	 vrho = f[3][i][j][k]*costh[i][k]+f[1][i][j][k]*sinth[i][k];
         vth  = -f[3][i][j][k]*sinth[i][k]+f[1][i][j][k]*costh[i][k];
         vphi = sqrt(pow(f[2][i][j][k],2)+vth*vth);           //sqrt(vfi*vfi+vth*vth)
         f[1][i][j][k] = vrho*sinth[i][k]+vphi*costh[i][k]*sin(chi[i][k]);
         f[2][i][j][k] = vphi*cos(chi[i][k]);
         f[3][i][j][k] = vrho*costh[i][k]-vphi*sinth[i][k]*sin(chi[i][k]);
         }*/
   return;
}
