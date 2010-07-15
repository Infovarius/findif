#define LEVEL

#include "head.h"

int main(int argc, char** argv)
{
   double t_cur, Ttot, dttry, dtdid, dtnext;
   int i,j,k,l, count;
   double temp, mf, mda, mdr, en, div, see[10][10];

   Re=10000;
   Gamma=1e-1;
   l1=3;
   l2=1;
   l3=1;

   nvar=4;
   n1=30;
   n2=1;
   n3=30;
   approx=7;
   ghost=(approx-1)/2;
   dx[0]=l1/n1;
   dx[1]=l2/n2;
   dx[2]=l3/n3;
   p1 = 8*l1/(l3*Re) ; p2 = 0;


   m1 = n1+2*ghost;
   m2 = n2+2*ghost;
   m3 = n3+2*ghost;
   mm1 = ghost+n1;
   mm2 = ghost+n2;
   mm3 = ghost+n3;


   f  =alloc_mem_4f(nvar, m1, m2, m3);
   f1 =alloc_mem_4f(nvar, m1, m2, m3);
   df =alloc_mem_4f(nvar, m1, m2, m3);
   df2=alloc_mem_4f(nvar, m1, m2, m3);
   df3=alloc_mem_4f(nvar, m1, m2, m3);
   df4=alloc_mem_4f(nvar, m1, m2, m3);
   df5=alloc_mem_4f(nvar, m1, m2, m3);
   nut=alloc_mem_3f(m1, m2, m3);

   init_conditions(f);
   boudary_conditions(f);

         div=0;
          for(i=ghost;i<mm1;i++)
          for(j=ghost;j<mm2;j++)
          for(k=ghost;k<mm3;k++) {
           temp=0;
           for(l=0;l<3;l++)
            temp+=dr(f[l],i,j,k,l+1,0,dx[l],ghost, approx);
           see[i-ghost][j-ghost]=temp;
           if (fabs(temp)>div) div=fabs(temp);
          }
         printf("%e %e %d %e\n", t_cur, dtdid, count, div);


   dtnext=1e-3;
   t_cur=0; Ttot=1000;
   count=0;

   while (t_cur < Ttot) {            /*----- MAIN ITERATIONS ----*/
        pde(t_cur, f, df);
        dttry=dtnext;
        timestep(f, df, t_cur, f1, dttry, &dtdid, &dtnext);
        t_cur+=dtdid;
        if (count%20==0) {
         clrscr();
         div=0;
         boudary_conditions(f1);
          for(i=ghost;i<mm1;i++)
          for(j=ghost;j<mm2;j++)
          for(k=ghost;k<mm3;k++) {
           temp=0;
           for(l=0;l<3;l++)
            temp+=dr(f1[l],i,j,k,l+1,0,dx[l],ghost, approx);
//            see[i-ghost][j-ghost]=temp;
           if (fabs(temp)>div) div=fabs(temp);
          }
         printf("%e %e %d %e\n", t_cur, dtdid, count, div);
        for(l=0;l<nvar;l++) {
          mf=mda=mdr=en=0;
          for(i=ghost;i<mm1;i++)
          for(j=ghost;j<mm2;j++)
          for(k=ghost;k<mm3;k++) {
            en+=f1[l][i][j][k]*f1[l][i][j][k];
            temp=fabs(f[l][i][j][k]-f1[l][i][j][k]);
            if (temp>mda) mda=temp;
            temp/=(f1[l][i][j][k]+1.e-100);
            if (temp>mdr) mdr=temp;
            if (fabs(f1[l][i][j][k])>mf) mf=fabs(f1[l][i][j][k]);
          }
         printf("%d %e %e %e %e\n",l, mf, mda, mdr, en);
        }
           for(k=0;k<m3;k++) printf("%e\n",f1[0][5][3][k]);
        }
        for(l=0;l<nvar;l++)
        for(i=ghost;i<mm1;i++)
        for(j=ghost;j<mm2;j++)
        for(k=ghost;k<mm3;k++)
           f[l][i][j][k]=f1[l][i][j][k];

        count++;
   }

   free_mem_4f(f  ,nvar, m1, m2, m3);
   free_mem_4f(f1 ,nvar, m1, m2, m3);
   free_mem_4f(df ,nvar, m1, m2, m3);
   free_mem_4f(df2,nvar, m1, m2, m3);
   free_mem_4f(df3,nvar, m1, m2, m3);
   free_mem_4f(df4,nvar, m1, m2, m3);
   free_mem_4f(df5,nvar, m1, m2, m3);
   free_mem_3f(nut, m1, m2, m3);

return 0;
}


