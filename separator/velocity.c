//----------------------- Calculation of known velocity profile  ----------//
//#define LEVEL extern
//#include "head.h"
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
// return fv2(t, 10., fv2(t,a1,amp1,amp2,t1) - fv2(0,a1,amp1,amp2,t1), fv2(t,dec,1,vlim,tpereg), tbr);
 return Vmax;
}

double chimod(double t,double phi)
{
// return( (th(-2*exp(-50*t)*(phi-40*M_PI*t))+1) / (th(80*exp(-50*t)*M_PI*t)+1) );
 return 1.;
}

double Gartmann(double rho,double rho0,double ksi)
{
 return( (ch(ksi)-ch(ksi*rho/rho0))/(ch(ksi)-1) );
}

double vfi_given(double t,double rho,double rho0)
{
const double ksi=0.001;
 return( vmod(t)*Gartmann(rho,rho0,ksi)/Vmax);
// return 1.;
}

double vtheta_given(double t,double rho,double rho0,double phi)
{
return( rho*chimod(t,phi)*vfi_given(t,rho,rho0)/rho0 );
}
