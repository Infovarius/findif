//----------------------- Calculation of PDE right part  ----------//
#define LEVEL extern
#include "head.h"

void pde(double t, double ****f, double ****df)
{
   int i,k,l,m;
   double dv1[3][3],dv2[3][3],dp1[3],dp2[3],dn1[3];

   boundary_conditions(f);

   for(i=ghost;i<mm1;i++)
   for(k=ghost;k<mm3;k++) {
      for(l=0;l<3;l++ for2D(l))
      for(m=0;m<3;m++ for2D(m)) {
         dv1[l][m]=dr(f[l],i,0,k,m+1,0,dx[m],ghost, approx);
         dv2[l][m]=dr(f[l],i,0,k,m+1,1,dx[m]*dx[m],ghost, approx);
      }
      for(m=0;m<3;m++ for2D(m)) {
         dp1[m]=dr(f[3],i,0,k,m+1,0,dx[m],ghost, approx);
         dn1[m]=dr(nut,i,0,k,m+1,0,dx[m],ghost, approx);
      }

      for(l=0;l<3;l++ for2D(l))
       df[l][i][0][k]= (dv2[l][0] + dv2[l][2])*nut[i][0][k]
                     - dp1[l] + (dn1[0]-f[0][i][0][k])*dv1[l][0]
                              + (dn1[2]-f[2][i][0][k])*dv1[l][2];
      df[3][i][0][k]= (-(dv1[0][0] + dv1[2][2]))/Gamma;
   }

   return;
}

void nut_by_flux(int ind) //calculating nu_turbulent by velocity fluctuations
{
double flux, maxflux = 0;
double maschtab = 10000;
int i,j,k,kol = 0;
for(i=0;i<=Nx+1;i++)
    for(j=0;j<=Ny+1;j++)
        for(k=0;k<=Nz+1;k++)
            {
            flux = 0;
            if(i>0)
            	{flux += modul(f[0][i][j][k]-f[0][i-1][j][k],f[1][i][j][k]-f[1][i-1][j][k],f[2][i][j][k]-f[2][i-1][j][k]);
            	kol++;
               }
            if(i<Nx+1)
            	{flux += modul(f[0][i][j][k]-f[0][i+1][j][k],f[1][i][j][k]-f[1][i+1][j][k],vz[ind][i][j][k]-f[2][i+1][j][k]);
            	kol++;
               }
            if(j>0)
            	{flux += modul(f[0][i][j][k]-f[0][i][j-1][k],f[1][i][j][k]-f[1][i][j-1][k],f[2][i][j][k]-f[2][i][j-1][k]);
            	kol++;
               }
            if(j<Ny+1)
	            {flux += modul(f[0][i][j][k]-f[0][i][j+1][k],f[1][i][j][k]-f[1][i][j+1][k],f[2][i][j][k]-f[2][i][j+1][k]);
            	kol++;
               }
            if(k>0)
            	{flux += modul(f[0][i][j][k]-f[0][i][j][k-1],f[1][i][j][k]-f[1][i][j][k-1],f[2][i][j][k]-f[2][i][j][k-1]);
            	kol++;
               }
            if(k<Nz+1)
            	{flux += modul(f[0][i][j][k]-f[0][i][j][k+1],f[1][i][j][k]-f[1][i][j][k+1],f[2][i][j][k]-f[2][i][j][k+1]);
            	kol++;
               }
            flux /= kol;
            if(flux>maxflux) maxflux = flux;
            nut[i][j][k] = (1. + maschtab * flux)/Re;
            }
}

void  boundary_conditions(double ****f)
{
   int i, k, l, g;

   // stream surfaces
   for(l=0;l<nvar;l++ for2D(l))
   for(k=ghost;k<mm3;k++)
   for(g=0;g<ghost;g++)
   {
   //periodic for velocities and gradient-periodic for pressure
      f[l][g][0][k] = f[l][n1+g][0][k] - ((l==3)?(p2-p1):0);
      f[l][mm1+g][0][k] = f[l][ghost+g][0][k] + ((l==3)?(p2-p1):0);
   }

/*   // vertical surfaces
   for(l=0;l<nvar;l++)
   for(i=ghost;i<mm1;i++)
   for(k=ghost;k<mm3;k++)
   for(g=0;g<ghost;g++)
   {
   //periodic conditions for velocities and pressure
      f[l][i][g][k] = f[l][i][n2+g][k];
      f[l][i][mm2+g][k] = f[l][i][ghost+g][k];
   }*/

   // on horizontal surfaces

   for(l=0;l<nvar;l++ for2D(l))
   for(i=ghost;i<mm1;i++)
   for(g=0;g<ghost;g++)
   {
   //sticking for velocities and free conditions for pressure
      f[l][i][0][g] = f[l][i][0][2*ghost-1-g]*((l==3)?1:-1);
      f[l][i][0][mm3+g] = f[l][i][0][mm3-1-g]*((l==3)?1:-1);
   }

   return;
}

void  init_conditions(double ****f)
{
   int i,k,l;
   double Noise=0.01, Noise1=0.;
   double k1,k2,k3;

   k1=2*M_PI/l1;  k3=M_PI/l3;
   for(i=0;i<m1;i++)
   for(k=0;k<m3;k++) {
        f[0][i][0][k]=coordin(k,2)*(l3-coordin(k,2))*4/l3/l3
          - Noise1*3*k3*pow(sin(k3*coordin(k,2) ),2.)*cos(k3*coordin(k,2))*sin(k1*coordin(i,0))
                     + Noise*((double)rand()-RAND_MAX/2)/RAND_MAX;
        f[2][i][0][k]=Noise1* pow(sin(k3*coordin(k,2)),3.) * k1*cos(k1*coordin(i,0))
                     + Noise*((double)rand()-RAND_MAX/2)/RAND_MAX;
        f[3][i][0][k]=p1+(i-0.5)*(p2-p1)/n1;
        nut[i][0][k]=1./Re;
   }

   return;
}

static double kf3[2][3][3]={{{-3./2.0, 2.0, -1./2.0}, {-1./2.0, 0.0, 1./2.0}, {1./2.0, -2.0, 3./2.0}},
                     {{1.0, -2.0, 1.0}, {1.0, -2.0, 1.0}, {1.0, -2.0, 1.0}}};
static double kf5[2][5][5]={{{-25./12.0, 4.0, -3.0, 4./3.0, -1./4.0}, {-1./4.0, -5./6.0, 3./2.0, -1./2.0, 1./12.0},
							{1./12.0, -2./3.0, 0.0, 2./3.0, -1./12.0}, {-1./12.0, 1./2.0, -3./2.0, 5./6.0, 1./4.0},
							{1./4.0, -4./3.0, 3.0, -4.0, 25./12.0}}, {{35./12.0, -26./3.0, 19./2.0, -14./3.0, 11./12.0},
                     {11./12.0, -5./3.0, 1./2.0, 1./3.0, -1./12.0}, {-1./12.0, 4./3.0, -5./2.0, 4./3.0, -1./12.0},
                     {-1./12.0, 1./3.0, 1./2.0, -5./3.0, 11./12.0}, {11./12.0, -14./3.0, 19./2.0, -26./3.0, 35./12.0}}};
static double kf7[2][7][7]={{{-49./20.0, 6.0, -15./2.0, 20./3.0, -15./4.0, 6./5.0, -1./6.0},	{-1./6.0, -77./60.0, 5./2.0, -5./3.0, 5./6.0,-1./4.0, 1./30.0},
                     {1./30.0, -2./5.0, -7./12.0, 4./3.0, -1./2.0, 2./15.0, -1./60.0}, {-1./60.0, 3./20.0, -3./4.0, 0.0, 3./4.0, -3./20.0, 1./60.0},
                     {1./60.0, -2./15.0, 1./2.0, -4./3.0, 7./12.0, 2./5.0, -1./30.0}, {-1./30.0, 1./4.0, -5./6.0, 5./3.0, -5./2.0, 77./60.0, 1./6.0},
							{1./6.0, -6./5.0, 15./4.0, -20./3.0, 15./2.0, -6.0, 49./20.0}}, {{203./45.0, -87./5.0, 117./4.0, -254./9.0, 33./2.0, -27./5.0, 137./180.0},
							{137./180.0, -49./60.0, -17./12.0, 47./18.0, -19./12.0, 31./60.0, -13./180.0}, {-13./180.0, 19./15.0, -7./3.0, 10./9.0, 1./12.0, -1./15.0, 1./90.0},
                     {1./90.0, -3./20.0, 3./2.0, -49./18.0, 3./2.0, -3./20.0, 1./90.0}, {1./90.0, -1./15.0, 1./12.0, 10./9.0, -7./3.0, 19./15.0, -13./180.0},
                     {-13./180.0, 31./60.0, -19./12.0, 47./18.0, -17./12.0, -49./60.0, 137./180.0}, {137./180.0, -27./5.0, 33./2.0, -254./9.0, 117./4.0, -87./5.0, 203./45.0}}};

double dr(double ***m, int ii, int jj, int kk, int dr, int or, double dx, int sh,  int sm)
/*        matrix     , point                 , direct, order , differ   , shift , sample */
/*                                           , 1,2,3 ,  0,1    dx,dx^2  , 0-left , 3,5,7 */
{
double tmp=0.0;
int i;

switch (sm*dr) {
	case 3 : for(i=0; i<sm; i++) tmp += m[ii+i-sh][jj][kk]*kf3[or][sh][i]; break;
  //	case 6 : for(i=0; i<sm; i++) tmp += m[ii][jj+i-sh][kk]*kf3[or][sh][i]; break;
	case 9 : for(i=0; i<sm; i++) tmp += m[ii][jj][kk+i-sh]*kf3[or][sh][i]; break;
	case 5 : for(i=0; i<sm; i++) tmp += m[ii+i-sh][jj][kk]*kf5[or][sh][i]; break;
  //	case 10: for(i=0; i<sm; i++) tmp += m[ii][jj+i-sh][kk]*kf5[or][sh][i]; break;
	case 15: for(i=0; i<sm; i++) tmp += m[ii][jj][kk+i-sh]*kf5[or][sh][i]; break;
	case 7 : for(i=0; i<sm; i++) tmp += m[ii+i-sh][jj][kk]*kf7[or][sh][i]; break;
  //	case 14: for(i=0; i<sm; i++) tmp += m[ii][jj+i-sh][kk]*kf7[or][sh][i]; break;
	case 21: for(i=0; i<sm; i++) tmp += m[ii][jj][kk+i-sh]*kf7[or][sh][i]; break;
	default :
    	nrerror("\nNO SUCH SAMPLE for derivative Bye ...");
	}
return(tmp/dx);
}

double coordin(int i, int dir)
                      //0-x,1-y,2-z
{
 return dir==1? 0 : dx[dir]*((double)(i-ghost)+0.5);
};

