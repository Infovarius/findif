#define LEVEL

#include "head.h"

int main(int argc, char** argv)
{
   double dttry, dtdid, dtnext;
   int i,j,k,l,i2,j2,k2;
   double ****ft;

 /* Initialize MPI */
 MPI_Init(&argc,&argv);
 /*  Get rank of this process and process group size */
 MPI_Comm_rank(MPI_COMM_WORLD,&rank);
 MPI_Comm_size(MPI_COMM_WORLD,&size);           putlog("main:reach this ",numlog++);

  srand(rank);
  NameMessageFile = "message.dat";
  NameErrorFile = "error.err";
  NameNuFile = "nut.dat";
  NameVFile  = "vv.dat";
  //NameDumpFile = "dump.dat";
  NameEnergyFile = "energy.dat";
  KaskadVarFile = "kaskvar.dat";
                                       numlog=0;
   init_param(argc,argv,&dtnext);       // initilization of parameters
   Gamma=1e-3;

   ghost=(approx-1)/2;            //radius of approx sample
   dx[0]=l1/N1;
   dx[1]=l2/N2;
   dx[2]=l3/N3;
   p1 = 4*l1/(l3*Re) ; p2 = 0;

   t_cur=0;
   count=0; enter = 0;

   nmessage("work has begun",0,0);

   sprintf(fname,"%s",(argc>1)? argv[1] : argv[0]);
   NameSnapFile = fname;
   sprintf(NameDumpFile ,"%s.dmp",fname);
   sprintf(NameStatFile,"%s_%d.sta",fname,rank);

   init_parallel();

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

   time_begin = MPI_Wtime();

   Master fileopen(NameErrorFile,0);

   init_conditions();

   boundary_conditions(f, nut);                           putlog("main:reach this ",numlog++);
//   dump(f,t_cur,count);                                 putlog("main:reach this ",numlog++);

   if(CheckStep!=0) check(f);                   putlog("main:reach this ",numlog++);
   if (OutStep!=0) printing(f,0,t_cur,count,PulsEnergy);   //in parallel version better not to do

/*------------------------ MAIN ITERATIONS -------------------------*/
   while ((Ttot==0 || t_cur < Ttot) && !razlet) {
        pde(t_cur, f, df);
        dttry=dtnext;
        timestep(f, df, nut, t_cur, f1, dttry, &dtdid, &dtnext);
        nut_by_flux(f,dtdid);
        t_cur+=dtdid;
        count++;
        if(t_cur >= Ttot && Ttot>0) break;
        if (CheckStep!=0 && count%CheckStep==0)
            {
	    boundary_conditions(f1,nut);
            check(f1);
            }
        if (OutStep!=0 && count%OutStep==0)
            {
            if (CheckStep!=0 && count%CheckStep!=0)
                {
		boundary_conditions(f1,nut);
                check(f1);
                }
	      else boundary_conditions(f1,nut);
            printing(f1,dtdid,t_cur,count,PulsEnergy);
            }
        if (SnapStep!=0 && count%SnapStep==0)
            snapshot(f,nut,t_cur,count);
	ft = f;  f = f1;  f1 = ft;
/*        if(kbhit())
             {
                switch (getch()) {
                        case 'd' : dump(f,t_cur,count);  break;
                        case 'q' : { dump(f,t_cur,count);
                                     MPI_Finalize();
                                     nrerror("You asked to exit. Here you are...",t_cur);
                                    }
                        }
              }*/
   }

   dump(f,nut,t_cur,count);
//   free_mem_2f(s_func,n3+2,kol_masht);
   free_mem_4f(f  ,nvar, m1, m2, m3);
   free_mem_4f(f1 ,nvar, m1, m2, m3);
   free_mem_4f(df ,nvar, m1, m2, m3);
   free_mem_4f(df2,nvar, m1, m2, m3);
   free_mem_4f(df3,nvar, m1, m2, m3);
   free_mem_4f(df4,nvar, m1, m2, m3);
   free_mem_4f(df5,nvar, m1, m2, m3);
   free_mem_3f(nut, m1, m2, m3);
   if(t_cur>=Ttot&&!razlet) nmessage("work is succesfully done",t_cur,count);
       else nrerror("this is break of scheme",t_cur,count);
   MPI_Finalize();
   nmessage("mpi_finalize is done",t_cur,count);
return 0;
}
