#define LEVEL extern
#include "head.h"

double dtfict=1e-5,dvr;

/*void main()
{
   double dp1;
   Gamma=1.;*/

/*   n1=20;
   n2=1;
   n3=20;
   m1 = n1+2*ghost;
   m2 = n2+2*ghost;
   m3 = n3+2*ghost;
   mm1 = ghost+n1;
   mm2 = ghost+n2;
   mm3 = ghost+n3;*/
/*   l1=3;
   l2=1;
   l3=1;
   dx[0]=l1/n1;
   dx[1]=l2/n2;
   dx[2]=l3/n3;
   t_cur = 0;*/
/*   f  =alloc_mem_4f(nvar, m1, m2, m3);   //f[3]-pressure,f[0..2]-v(vector)
   f1 =alloc_mem_4f(nvar, m1, m2, m3);
   df =alloc_mem_4f(nvar, m1, m2, m3);
   diver=alloc_mem_3f(m1, m2, m3);*/

/*   double maxp;
   dvr = divergention(f);
   init_conditions(f,10000);
   dump("vv.dat",f[0],0);
   dump("vv.dat",f[2],1);
   boundary_conditions(f);
   dump("vv.dat",f[0],0);
   dump("vv.dat",f[2],1);
   maxp = MaxArray(&f[3][0][0][0],m1*m2*m3);
   dvr = divergention(f);
   dump("diver.dat",diver,0);
   p_relax(f);
   for(int i=ghost;i<mm1;i++)
   for(int j=ghost;j<=ghost;j++)
   for(int k=ghost;k<mm3;k++)
      for(int m=0;m<3;m++ for2D(m)) {
         dp1=dr(p[3],i,j,k,m+1,0,dx[m],ghost, approx);
         f[m][i][j][k] -= dp1;
                               }

   dump("vv.dat",f[0],0);
   dump("vv.dat",f[2],1);
   dvr = divergention(f);
   dump("diver.dat",diver,1);*/

/*   free_mem_4f(f  ,nvar, m1, m2, m3);
   free_mem_4f(f1 ,nvar, m1, m2, m3);
   free_mem_4f(df ,nvar, m1, m2, m3);
   free_mem_3f(diver, m1, m2, m3);*/
//}

double divergention(double ****f)
{
double sumdiv=0;
int i,j,k,m;
boundary_conditions(f);
   for(i=ghost;i<mm1;i++)
   for(j=ghost;j<=ghost;j++)
   for(k=ghost;k<mm3;k++)   {
      diver[i][j][k] = 0;
      for(m=0;m<3;m++ for2D(m))
         diver[i][j][k] += dr(f[m],i,j,k,m+1,0,dx[m],ghost, approx);
//                          +dr(df[m],i,j,k,m+1,0,dx[m],ghost, approx);
//    if(fabs(diver[i][j][k])>maxdiv)
      sumdiv += fabs(diver[i][j][k]);
                             }
return(sumdiv);
}

double MaxArray(double *f,long kol)
{
 double dop=0;
 while(kol>0)
  {
  if (fabs(*f) > dop) dop = fabs(*f);
  f++;
  kol--;
  }
 return dop;
}

void p_relax(double ****f,double ****df)
{
int i,j,k,m;
double  divpar=1e-6,MinDiv=1e-3, maxp, maxd;
double dttry, dtdid, dtnext;

pres_step = 0;
   for(i=0;i<m1;i++)
   for(j=0;j<m2;j++)
   for(k=0;k<m3;k++)
      for(m=0;m<nvar;m++)
       prel[m][i][j][k] = f[m][i][j][k];
//dump("p.dat",p[3],pres_step+1);
 dtnext= 1e-4;
  do {
   dvr = divergention(f);
   maxp = MaxArray(&prel1[3][0][0][0],m1*m2*m3);
   maxd=pde_prelax(t_cur, prel, df);
   dttry = dtnext;
   timestep(prel,df,t_cur,prel1,dttry,dtdid,dtnext,pde_prelax);
   t_cur += dtdid;

   for(i=0;i<m1;i++)
   for(j=0;j<m2;j++)
   for(k=0;k<m3;k++)
       prel[3][i][j][k] = prel1[3][i][j][k];
   boundary_conditions(prel);
   for(int i=ghost;i<mm1;i++)
   for(int j=ghost;j<=ghost;j++)
   for(int k=ghost;k<mm3;k++)
      for(int m=0;m<3;m++ for2D(m))
         prel[m][i][j][k] = f[m][i][j][k]-dr(prel[3],i,j,k,m+1,0,dx[m],ghost, approx);
   dvr = divergention(prel);
   maxp = MaxArray(&prel[3][0][0][0],m1*m2*m3);
   if(pres_step%OutStep==0)
      {
//        dump("vv.dat",p[0],0);
//        dump("vv.dat",p[2],1);
//        dump("diver.dat",diver,1);
//        dump("p.dat",p[3],pres_step);
     clrscr();
     printf("T=%0.4lf; dt=%0.10lf\n",t_cur,dtdid);
     printf("Divergention on %d step = %0.10lf\n",pres_step,dvr/n1/n3);
     printf("Difference between steps=%0.10lf\n",maxd*dtdid);
     printf("Maximum pressure=%0.10lf\n",maxp);
      }

   pres_step++;
   } while (/*maxd*dtdid>divpar &&*/ dvr/n1/n2>MinDiv);
for(i=0;i<m1;i++)
 for(j=0;j<m2;j++)
   for(k=0;k<m3;k++)
      for(m=0;m<3;m++)
       df[m][i][j][k] = prel[m][i][j][k] - f[m][i][j][k];
}

double pde_prelax(double t,double ****f,double ****df)
{
int i,j,k,m;
double dp2[3];
double maxd=0;
boundary_conditions(f);
   for(i=ghost;i<mm1;i++)
   for(j=ghost;j<=ghost;j++)
   for(k=ghost;k<mm3;k++)
     {
     for(m=0;m<3;m++ for2D(m))
         {
         dp2[m]=dr(f[3],i,j,k,m+1,1,dx[m]*dx[m],ghost, approx);
//         if(dp2[m]>maxd) maxd = dp2[m];
         }
     df[3][i][j][k]= dp2[0]+dp2[2]-diver[i][j][k]/Gamma;
     if(fabs(df[3][i][j][k])>maxd) maxd = fabs(df[3][i][j][k]);
     }
return maxd;
}
