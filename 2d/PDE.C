//----------------------- Calculation of PDE right part  ----------//
#define LEVEL extern
#include "head.h"

void pde(double t, double ****f, double ****df)
{
   int i,j,k,l,m;
   double dv1[3][3],dv2[3][3],dp1[3],dp2[3],dn1[3];

   boudary_conditions(f);

   for(i=ghost;i<mm1;i++)
   for(j=ghost;j<mm2;j++)
   for(k=ghost;k<mm3;k++) {
      for(l=0;l<3;l++)
      for(m=0;m<3;m++) {
         if(l==1 || m==1) continue;
         dv1[l][m]=dr(f[l],i,j,k,m+1,0,dx[m],ghost, approx);
         dv2[l][m]=dr(f[l],i,j,k,m+1,1,dx[m]*dx[m],ghost, approx);
      }
      for(m=0;m<3;m++) {
         if(l==1 || m==1) continue;
         dp1[m]=dr(f[3],i,j,k,m+1,0,dx[m],ghost, approx);
         dp2[m]=dr(f[3],i,j,k,m+1,1,dx[m]*dx[m],ghost, approx);
         dn1[m]=dr(nut,i,j,k,m+1,0,dx[m],ghost, approx);
      }

      for(l=0;l<3;l++) {
               if(l==1 || m==1) continue;
       df[l][i][j][k]= (dv2[l][0]  + dv2[l][2])*nut[i][j][k]
                     - dp1[l] + (dn1[0]-f[0][i][j][k])*dv1[l][0]
                              + (dn1[2]-f[2][i][j][k])*dv1[l][2];
      }
      df[3][i][j][k]= (-(dv1[0][0]  + dv1[2][2]) +
                       0*(dp2[0] + dp2[2]))/Gamma;
   }

   return;
}

void  boudary_conditions(double ****f)
{
   int i, j, k, l, g;

   // 1 boundary
   for(l=0;l<nvar;l++)
   for(j=ghost;j<mm2;j++)
   for(k=ghost;k<mm3;k++)
   for(g=0;g<ghost;g++)
   {
      f[l][g][j][k] = f[l][n1+g][j][k] - ((l==3)?(p2-p1):0);
      f[l][n1+ghost+g][j][k] = f[l][ghost+g][j][k] + ((l==3)?(p2-p1):0);
   }


   // 3 boundary
   for(l=0;l<nvar;l++)
   for(i=ghost;i<mm1;i++)
   for(j=ghost;j<mm2;j++)
   for(g=0;g<ghost;g++)
   {
      f[l][i][j][g] = f[l][i][j][2*ghost-1-g]*((l==3)?1:-1);
      f[l][i][j][n3+ghost+g] = f[l][i][j][n3+ghost-1-g]*((l==3)?1:-1);
   }

   return;
}

void  init_conditions(double ****f)
{
   int i,j,k,l;
   double Noise=0., Noise1=0.01;

//   for(l=0;l<nvar;l++)
   for(i=0;i<m1;i++)
   for(j=0;j<m2;j++)
   for(k=0;k<m3;k++) {
        f[0][i][j][k]=coordin(k,2)*(l3-coordin(k,2))*4
                - Noise1*(l3-2*coordin(k,2))*sin(coordin(i,0))
                + Noise*((double)rand()-RAND_MAX/2)/RAND_MAX;
        f[2][i][j][k]=Noise1*coordin(k,2)*(l3-coordin(k,2))*cos(coordin(i,0))
                + Noise*((double)rand()-RAND_MAX/2)/RAND_MAX;
        f[3][i][j][k]=p1+(i-0.5)*(p2-p1)/n1;
        nut[i][j][k]=1./Re;
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
	case 6 : for(i=0; i<sm; i++) tmp += m[ii][jj+i-sh][kk]*kf3[or][sh][i]; break;
	case 9 : for(i=0; i<sm; i++) tmp += m[ii][jj][kk+i-sh]*kf3[or][sh][i]; break;
	case 5 : for(i=0; i<sm; i++) tmp += m[ii+i-sh][jj][kk]*kf5[or][sh][i]; break;
	case 10: for(i=0; i<sm; i++) tmp += m[ii][jj+i-sh][kk]*kf5[or][sh][i]; break;
	case 15: for(i=0; i<sm; i++) tmp += m[ii][jj][kk+i-sh]*kf5[or][sh][i]; break;
	case 7 : for(i=0; i<sm; i++) tmp += m[ii+i-sh][jj][kk]*kf7[or][sh][i]; break;
	case 14: for(i=0; i<sm; i++) tmp += m[ii][jj+i-sh][kk]*kf7[or][sh][i]; break;
	case 21: for(i=0; i<sm; i++) tmp += m[ii][jj][kk+i-sh]*kf7[or][sh][i]; break;
	default :
    	nrerror("\nNO SUCH SAMPLE for derivative Bye ...");
	}
return(tmp/dx);
}

double coordin(int i, int dir)
{
 return dx[dir]*((double)(i-ghost)+0.5);
};

