//----------------------- Adopting Runge-Cutta method ----------//
#define LEVEL extern
#include "head.h"

#define Safety       (double)0.9
#define Pgrow        (double)-0.2
#define Pshrink      (double)-0.25
#define ErrCon       (double)1.89e-4
#define ErrScale     (double)1e4
#define MinScale     (double)1e-12


void timestep(double ****f, double ****df, double ***nut, double t, double ****fout,
              double dttry, double *dtdid, double *dtnext)
{
   double dt, err, errs;

   dt=dttry;
   for (;;)
   {
      err = rkck(f, df, nut, t, dt, fout);
      if (VarStep==0 || count%VarStep!=0) {
            *dtdid = dt;
            *dtnext = dt;
            return;
            }

      MPI_Allreduce(&err, &errs, 1, MPI_DOUBLE , MPI_MAX, MPI_COMM_WORLD);
      err=errs;

      err *=ErrScale;
      if (err<=1.0) break;
      dt = max(Safety * dt * pow(err, Pshrink), dt*(double)0.1);
      if (t+dt == t)
	    {
	     dump(f,nut,t,count);
             nrerror("Stepsize underflow in rk\n",t,count);
            }
   }
   *dtdid = dt;
   if (err>ErrCon)
      *dtnext = Safety * dt * pow(err, Pgrow);
   else
      *dtnext = dt * (double)5.0;
}

double rkck(double ****f, double ****df1, double ***nut, double t, double dt, double ****fout)
{
static double   a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
		b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
		b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
		b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
		b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
		c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
		dc5 = -277.00/14336.0;
double          dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
		dc4=c4-13525.0/55296.0,dc6=c6-0.25;
   int i,j,k,l;
   double err, sca, err1;
enter++;
   /*2rd step*/
   for(l=0;l<nvar;l++)
   for(i=ghost;i<mm1;i++)
   for(j=ghost;j<mm2;j++)
   for(k=ghost;k<mm3;k++)
//     if(isType(node[i][k],NodeFluid))
      fout[l][i][j][k] = f[l][i][j][k]+dt*b21*df1[l][i][j][k];
   pde(t+a2*dt,fout, df2);

   /*3rd step*/
   for(l=0;l<nvar;l++)
   for(i=ghost;i<mm1;i++)
   for(j=ghost;j<mm2;j++)
   for(k=ghost;k<mm3;k++)
//     if(isType(node[i][k],NodeFluid))
      fout[l][i][j][k] = f[l][i][j][k]+dt*(b31*df1[l][i][j][k]+
                                           b32*df2[l][i][j][k]);
   pde(t+a3*dt,fout, df3);

   /*4th step*/
   for(l=0;l<nvar;l++)
   for(i=ghost;i<mm1;i++)
   for(j=ghost;j<mm2;j++)
   for(k=ghost;k<mm3;k++)
//     if(isType(node[i][k],NodeFluid))
      fout[l][i][j][k] = f[l][i][j][k]+dt*(b41*df1[l][i][j][k]+
                                           b42*df2[l][i][j][k]+
                                           b43*df3[l][i][j][k]);
   pde(t+a4*dt,fout, df4);

   /*5th step*/
   for(l=0;l<nvar;l++)
   for(i=ghost;i<mm1;i++)
   for(j=ghost;j<mm2;j++)
   for(k=ghost;k<mm3;k++)
//     if(isType(node[i][k],NodeFluid))
      fout[l][i][j][k] = f[l][i][j][k]+dt*(b51*df1[l][i][j][k]+
                                           b52*df2[l][i][j][k]+
                                           b53*df3[l][i][j][k]+
                                           b54*df4[l][i][j][k]);
   pde(t+a5*dt,fout, df5);

   /*6th step*/
   for(l=0;l<nvar;l++)
   for(i=ghost;i<mm1;i++)
   for(j=ghost;j<mm2;j++)
   for(k=ghost;k<mm3;k++)
//     if(isType(node[i][k],NodeFluid))
      fout[l][i][j][k] = f[l][i][j][k]+dt*(b61*df1[l][i][j][k]+
                                           b62*df2[l][i][j][k]+
                                           b63*df3[l][i][j][k]+
                                           b64*df4[l][i][j][k]+
                                           b65*df5[l][i][j][k]);
   pde(t+a6*dt,fout, df2);

 /*calculating output matrix and error value*/
   err = 0.0;
   for(l=0;l<nvar;l++)
   for(i=ghost;i<mm1;i++)
   for(j=ghost;j<mm2;j++)
   for(k=ghost;k<mm3;k++)
//     if(isType(node[i][k],NodeFluid) && !isType(node[i][k],NodeClued))
   {
      fout[l][i][j][k] = f[l][i][j][k]+dt*(c1*df1[l][i][j][k]+
					   c3*df3[l][i][j][k]+
					   c4*df4[l][i][j][k]+
					   c6*df2[l][i][j][k]);
      sca = fabs(f[l][i][j][k])+fabs(dt*df1[l][i][j][k])+MinScale;
      err1 = fabs(dt*(    dc1*df1[l][i][j][k]+
			  dc3*df3[l][i][j][k]+
			  dc4*df4[l][i][j][k]+
			  dc5*df5[l][i][j][k]+
			  dc6*df2[l][i][j][k]))/sca;
      err = max(err, err1);
   }
   return err;
}
