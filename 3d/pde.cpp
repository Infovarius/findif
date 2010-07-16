//----------------------- Calculation of PDE right part  ----------//
#define LEVEL extern
#include "head.h"

inline double norma(double a,double b,double c,int order)
{
if(order==2) return ((a)*(a)+(b)*(b)+(c)*(c));
   else  return  pow(((a)*(a)+(b)*(b)+(c)*(c)),order/2.);
}

void pde(double t, double ****f, double ****df)
{
   int i,j,k,l,m;
   double dv1[3][3],dv2[3][3],dp1[3],dp2[3],dn1[3];

   boundary_conditions(f);

   for(i=ghost;i<mm1;i++)
   for(j=ghost;j<mm2;j++)
   for(k=ghost;k<mm3;k++) {
      for(l=0;l<3;l++)
      for(m=0;m<3;m++) {
         dv1[l][m]=dr(f[l],i,j,k,m+1,0,dx[m],ghost, approx);
         dv2[l][m]=dr(f[l],i,j,k,m+1,1,dx[m]*dx[m],ghost, approx);
      }
      for(m=0;m<3;m++) {
         dp1[m]=dr(f[3],i,j,k,m+1,0,dx[m],ghost, approx);
         dn1[m]=dr(nut,i,j,k,m+1,0,dx[m],ghost, approx);
      }

      for(l=0;l<3;l++)
       df[l][i][j][k]= (dv2[l][0] + dv2[l][1] + dv2[l][2])*nut[i][j][k]
                     - dp1[l] + (dn1[0]-f[0][i][j][k])*dv1[l][0]
                              + (dn1[1]-f[1][i][j][k])*dv1[l][1]
                              + (dn1[2]-f[2][i][j][k])*dv1[l][2];
      df[3][i][j][k]= (-(dv1[0][0] + dv1[1][1] + dv1[2][2]))/Gamma;
   }

   return;
}

double deviation(double ****f,int i,int j,int k)
{
const size_okr=min(1,ghost);
double flux;
int kol=0,l;
   for(l=1;l<=size_okr;l++)
            {
            flux = 0;
            if(i>0)
               {flux += norma(f[0][i][j][k]-f[0][i-l][j][k],f[1][i][j][k]-f[1][i-l][j][k],f[2][i][j][k]-f[2][i-l][j][k],2);
            	kol++;
               }
            if(i<m1)
               {flux += norma(f[0][i][j][k]-f[0][i+l][j][k],f[1][i][j][k]-f[1][i+l][j][k],f[2][i][j][k]-f[2][i+l][j][k],2);
            	kol++;
               }
            if(j>0)
               {flux += norma(f[0][i][j][k]-f[0][i][j-l][k],f[1][i][j][k]-f[1][i][j-l][k],f[2][i][j][k]-f[2][i][j-l][k],2);
            	kol++;
               }
            if(j<m2)
               {flux += norma(f[0][i][j][k]-f[0][i][j+l][k],f[1][i][j][k]-f[1][i][j+l][k],f[2][i][j][k]-f[2][i][j+l][k],2);
            	kol++;
               }
            if(k>0)
               {flux += norma(f[0][i][j][k]-f[0][i][j][k-l],f[1][i][j][k]-f[1][i][j][k-l],f[2][i][j][k]-f[2][i][j][k-l],2);
                kol++;
               }
            if(k<m3)
               {flux += norma(f[0][i][j][k]-f[0][i][j][k+l],f[1][i][j][k]-f[1][i][j][k+l],f[2][i][j][k]-f[2][i][j][k+l],2);
                kol++;
               }
           };
   flux /= kol;
return(flux);
}

void nut_by_flux(double ****f,double dt) //calculating nu_turbulent by velocity fluctuations
{
double maschtab = 100000;
int i,j,k,l;
double koef;
struct_func(f,2,2,3);
/*for(i=0;i<n3;i++)
    {
    koef=sqrt(s_func[i][0]/(pow(sha[i][1],2.)+pow(shb[i][1],2.)));
    sha[i][1] *= koef;
    shb[i][1] *= koef;
    koef=sqrt(s_func[i][1]/(pow(sha[i][0],2.)+pow(shb[i][0],2.)));
    sha[i][0] *= koef;
    shb[i][0] *= koef;
    }*/
/*clrscr();
for(j=0;j<n3;j++)
    {
    printf("%lf  %lf",s_func[j][0],s_func[j][1]);
    double en;
    for (i=0,en=0; i<=Ns; i++)
      en+=sha[j][i]*sha[j][i]+shb[j][i]*shb[j][i];
    printf("   totEn=%lf\n",en);
    }  */
time_step_shell(dt);
for(k=0;k<n3;k++)
   {
   double tmp = maschtab*pow(
           (nl[2]*s_func[k][0] + nl[1]*s_func[k][1] + nl[0]*s_func[k][2])*pow(dx[2],4),
                         1./3);
   for(i=ghost;i<mm1;i++)
       for(j=ghost;j<mm2;j++)
            nut[i][j][k+ghost] = (1. + tmp)/Re;
   }
}

void  boundary_conditions(double ****f)
{
   int i, j, k, l, g;

   // stream surfaces
   for(l=0;l<nvar;l++)
   for(j=ghost;j<mm2;j++)
   for(k=ghost;k<mm3;k++)
   for(g=0;g<ghost;g++)
   {
   //periodic for velocities and gradient-periodic for pressure
      f[l][g][j][k] = f[l][n1+g][j][k] - ((l==3)?(p2-p1):0);
      f[l][mm1+g][j][k] = f[l][ghost+g][j][k] + ((l==3)?(p2-p1):0);
   }

   // vertical surfaces
   for(l=0;l<nvar;l++)
   for(i=ghost;i<mm1;i++)
   for(k=ghost;k<mm3;k++)
   for(g=0;g<ghost;g++)
   {
   //periodic conditions for velocities and pressure
      f[l][i][g][k] = f[l][i][n2+g][k];
      f[l][i][mm2+g][k] = f[l][i][ghost+g][k];
   }

   // on horizontal surfaces

   for(l=0;l<nvar;l++)
   for(i=ghost;i<mm1;i++)
   for(j=ghost;j<mm2;j++)
   for(g=0;g<ghost;g++)
   {
   //sticking for velocities and free conditions for pressure
      f[l][i][j][g] = f[l][i][j][2*ghost-1-g]*((l==3)?1:-1);
      f[l][i][j][mm3+g] = f[l][i][j][mm3-1-g]*((l==3)?1:-1);
   }

   return;
}

void  init_conditions(double ****f,double Re)
{
   int i,j,k,l;
   double Noise=0.3, Noise1=0;
//   double k1,k2,k3;

//   k1=2*M_PI/l1;  k3=M_PI/l3;

   for(i=0;i<m1;i++)
   for(j=0;j<m2;j++)
   for(k=0;k<m3;k++) {
        f[0][i][j][k]=(1+Noise*((double)rand()-RAND_MAX/2)/RAND_MAX)*
                       coordin(k,2)*(l3-coordin(k,2))*4/l3/l3;
        f[1][i][j][k]=Noise1*cos(2*M_PI*coordin(j,1)/l2)*cos(2*M_PI*coordin(k,2)/l3)
                      + Noise*((double)rand()-RAND_MAX/2)/RAND_MAX*
                       coordin(k,2)*(l3-coordin(k,2))*4/l3/l3;
        f[2][i][j][k]=Noise1*sin(2*M_PI*coordin(j,1)/l2)*sin(2*M_PI*coordin(k,2)/l3)
                      + Noise*((double)rand()-RAND_MAX/2)/RAND_MAX*
                       coordin(k,2)*(l3-coordin(k,2))*4/l3/l3;
        f[3][i][j][k]=p1+(i-0.5)*(p2-p1)/n1;
        nut[i][j][k]=(
        (0.39+14.8*exp(-2.13*pow(2*coordin(k,2)-l3,2)))*1
                +1)/Re;
   }
//   struct_func(f,2,2,3);
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
    	nrerror("\nNO SUCH SAMPLE for derivative Bye ...",0);
	}
return(tmp/dx);
}

double coordin(int i, int dir)
                      //0-x,1-y,2-z
{
 return dx[dir]*((double)(i-ghost)+0.5);
};

