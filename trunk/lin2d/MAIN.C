#define LEVEL

//#include <conio.h>
#include "head.h"

#define EACH 100
int main(int argc, char** argv)
{
   double t_cur, Ttot, dttry, dtdid, dtnext;
   int i,k,l;

   Re=10000;
   Gamma=1e-3;
   l1=3;
   l2=1;
   l3=1;

   nvar=4;
   n1=10;
   n2=1;
   n3=10;
   approx=7;                      //derivatives approximation order
   ghost=(approx-1)/2;            //radius of approx sample
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

   t_cur=0; Ttot=1000;
   count=0; enter = 0;

   f  =alloc_mem_4f(nvar, m1, m2, m3);   //f[3]-pressure,f[0,2]-v(vector)
   f1 =alloc_mem_4f(nvar, m1, m2, m3);
   df =alloc_mem_4f(nvar, m1, m2, m3);
   df2=alloc_mem_4f(nvar, m1, m2, m3);
   df3=alloc_mem_4f(nvar, m1, m2, m3);
   df4=alloc_mem_4f(nvar, m1, m2, m3);
   df5=alloc_mem_4f(nvar, m1, m2, m3);
   nut=alloc_mem_3f(m1, m2, m3);

   init_conditions(f);
   printing(f,dtdid,t_cur,count);

   dtnext=1e-3;
   dump(f,t_cur,count);
   while (t_cur < Ttot) {            /*----- MAIN ITERATIONS ----*/
        pde(t_cur, f, df);
        dttry=dtnext;
        timestep(f, df, t_cur, f1, dttry, &dtdid, &dtnext);
        t_cur+=dtdid;
        if (count%EACH==0)
            printing(f1,dtdid,t_cur,count);
        for(l=0;l<nvar;l++)
        for(i=ghost;i<mm1;i++)
        for(k=ghost;k<mm3;k++)
           f[l][i][0][k]=f1[l][i][0][k];

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


