#define LEVEL

#include "head.h"
#include <time.h>

int main(int argc, char** argv)
{
   double dttry, dtdid, dtnext;
   int i,j,k,l;

   time(&time_begin);
   nmessage("work has begun",0);

   OutStep = 10*(CheckStep=10);
#include "init_vars.h"


   Re=1000000;
   Gamma=1e-4;
   l1=3;
   l2=1;
   l3=1;

   nvar=4;
   n1=16;
   n2=16;
   n3=16;
   Ns=15;
   approx=7;                      //derivatives approximation order
   ghost=(approx-1)/2;            //radius of approx sample
   dx[0]=l1/n1;
   dx[1]=l2/n2;
   dx[2]=l3/n3;
   p1 = 80*l1/(l3*Re) ; p2 = 0;


   m1 = n1+2*ghost;
   m2 = n2+2*ghost;
   m3 = n3+2*ghost;
   mm1 = ghost+n1;
   mm2 = ghost+n2;
   mm3 = ghost+n3;

   t_cur=0; Ttot=1000;
   count=0; enter = 0;

   s_func = alloc_mem_2f(n3+2,kol_masht);
   f  =alloc_mem_4f(nvar, m1, m2, m3);   //f[3]-pressure,f[0..2]-v(vector)
   f1 =alloc_mem_4f(nvar, m1, m2, m3);
   df =alloc_mem_4f(nvar, m1, m2, m3);
   df2=alloc_mem_4f(nvar, m1, m2, m3);
   df3=alloc_mem_4f(nvar, m1, m2, m3);
   df4=alloc_mem_4f(nvar, m1, m2, m3);
   df5=alloc_mem_4f(nvar, m1, m2, m3);
   nut=alloc_mem_3f(m1, m2, m3);
   init_shell();

   fileopen("error.err",0);
   init_conditions(f,Re);
   if(argc>1)
           {
           FILE *inp=fileopen(argv[1],-1);
           float tmpd;
           char tmpc;
           for(l=0;l<nvar;l++)
              {
              do fscanf(inp,"%c",&tmpc); while (tmpc!='{');
              for(i=0;i<m1;i++)
                 {
                 do fscanf(inp,"%c",&tmpc); while (tmpc!='{');
                 for(j=0;j<m2;j++)
                    {
                    do fscanf(inp,"%c",&tmpc); while (tmpc!='{');
                    for(k=0;k<m3;k++)
                       {
                       tmpc=fscanf(inp,"%g",&tmpd);
                       fscanf(inp,"%c",&tmpc);
                       f[l][i][j][k]=tmpd;
                       }
                    fscanf(inp,"%c",&tmpc);
                    }
                 fscanf(inp,"%c",&tmpc);
                 }
              fscanf(inp,"%c",&tmpc);
              }
           }
   PulsEn=check(f);
   boundary_conditions(f);
   printing(f,0,t_cur,count,PulsEn);

   dtnext=1e-3;
   dump(f,t_cur,count);

/*------------------------ MAIN ITERATIONS -------------------------*/
   while (t_cur < Ttot && !razlet) {
        pde(t_cur, f, df);
        dttry=dtnext;
        timestep(f, df, t_cur, f1, dttry, &dtdid, &dtnext);
        nut_by_flux(f,dtdid);
        t_cur+=dtdid;
        count++;
        if (count%CheckStep==0)
            PulsEn=check(f);
        if (count%OutStep==0)
            {
            if (count%CheckStep!=0)
                PulsEn=check(f);
            printing(f1,dtdid,t_cur,count,PulsEn);
            };
        for(l=0;l<nvar;l++)
        for(i=ghost;i<mm1;i++)
        for(j=ghost;j<mm2;j++)
        for(k=ghost;k<mm3;k++)
           f[l][i][j][k]=f1[l][i][j][k];
        if(kbhit())
             {
                switch (getch()) {
                        case 'd' : dump(f,t_cur,count);  break;
                        case 'q' : { dump(f,t_cur,count); 
                                     nrerror("You asked to exit. Here you are...",t_cur);
                                    }
                        }
              }
   }

   dump(f,t_cur,count); 
//   free_mem_2f(s_func,n3+2,kol_masht);
   free_mem_4f(f  ,nvar, m1, m2, m3);
   free_mem_4f(f1 ,nvar, m1, m2, m3);
   free_mem_4f(df ,nvar, m1, m2, m3);
   free_mem_4f(df2,nvar, m1, m2, m3);
   free_mem_4f(df3,nvar, m1, m2, m3);
   free_mem_4f(df4,nvar, m1, m2, m3);
   free_mem_4f(df5,nvar, m1, m2, m3);
   free_mem_3f(nut, m1, m2, m3);
   if(t_cur>Ttot&&!razlet) nmessage("work is succesfully done",t_cur);
       else nrerror("this is break of scheme",t_cur);

return 0;
}


