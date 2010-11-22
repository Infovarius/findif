#define LEVEL

#include "head.h"

int main(int argc, char** argv)
{
   double dttry, dtdid, dtnext;
   int i,j,k,l,i2,j2,k2;                           FILE *fd;

 /* Initialize MPI */
 MPI_Init(&argc,&argv);
 /*  Get rank of this process and process group size */
 MPI_Comm_rank(MPI_COMM_WORLD,&rank);
 MPI_Comm_size(MPI_COMM_WORLD,&size);           putlog("main:reach this ",numlog++);

  NameMessageFile = "message.dat";
  NameErrorFile = "error.dat";
  NameNuFile = "nut.dat";
  NameVFile  = "vv.dat";
  //NameDumpFile = "dump.dat";
  NameEnergyFile = "energy.dat";
  KaskadVarFile = "kaskvar.dat";
                                       numlog=0;
   init_param(argc,argv,&dtnext);       // initilization of parameters
   Gamma=1e-3;

   ghost=(approx-1)/2;            //radius of approx sample
   dx[0]=2*R/N1;
   dx[1]=lfi/N2;
   dx[2]=2*R/N3;
   omega=1;

   t_cur=0;
   count=0; enter = 0;

   nmessage("work has begun",0);

   sprintf(fname,"%s",(argc>1)? argv[1] : argv[0]);
   NameSnapFile = fname;
   sprintf(NameDumpFile ,"%s.dmp",fname);
   sprintf(fname_stat,"%s_%d.sta",fname,rank);

   init_parallel();

   s_func = alloc_mem_2f(n3+2,kol_masht);
   f  =alloc_mem_4f(nvar, m1, m2, m3);   //f[3]-pressure,f[0..2]-v(vector)
   f1 =alloc_mem_4f(nvar, m1, m2, m3);
   df =alloc_mem_4f(nvar, m1, m2, m3);
   df2=alloc_mem_4f(nvar, m1, m2, m3);
   df3=alloc_mem_4f(nvar, m1, m2, m3);
   df4=alloc_mem_4f(nvar, m1, m2, m3);
   df5=alloc_mem_4f(nvar, m1, m2, m3);
   r_1=alloc_mem_1f(m1);                  // r^(-1)
   r_2=alloc_mem_1f(m1);                  // r^(-2)
   node=alloc_mem_2f(m1, m3);         // kind of nodes
   refr=alloc_mem_2f(m1, m3);         // reflection of nodes relative circle
   refz=alloc_mem_2f(m1, m3);
   nut=alloc_mem_3f(m1, m2, m3);
   averf = alloc_mem_3f(3,m1,m3);

   init_shell();

   time_begin = MPI_Wtime();

   Master fileopen("error.err",0);

   init_conditions(f,Re);
   for(i=0;i<m1;i++) { r_1[i]=1./coordin(i,0); r_2[i]=r_1[i]*r_1[i]; }
//--------------------------------------
    if(rank!=0) MPI_Recv("",0,MPI_CHAR,rank-1,1,MPI_COMM_WORLD,statuses);

 fd=fileopen("debug",rank);

 Master nmessage("debug has been started",t_cur);

 print_array2d(fd,node,0,m1,0,m3);
 print_array2d(fd,refr,0,m1,0,m3);
 print_array2d(fd,refz,0,m1,0,m3);
 fclose(fd);

 if(rank!=size-1) MPI_Send("",0,MPI_CHAR,rank+1,1,MPI_COMM_WORLD);
             else nmessage("dump is done",t_cur);
//--------------------------------------

   boundary_conditions(f);                             putlog("main:reach this ",numlog++);
   dump(f,t_cur,count);                                 putlog("main:reach this ",numlog++);

   if(CheckStep!=0) check(f);                            putlog("main:reach this ",numlog++);
   if (OutStep!=0) printing(f,0,t_cur,count,PulsEnergy);   //in parallel version better not to do

/*------------------------ MAIN ITERATIONS -------------------------*/
   while (t_cur < Ttot && !razlet) {
        pde(t_cur, f, df);
        dttry=dtnext;
        timestep(f, df, t_cur, f1, dttry, &dtdid, &dtnext);
        nut_by_flux(f,dtdid);
        t_cur+=dtdid;
        count++;
        if (CheckStep!=0 && count%CheckStep==0)
            {
            boundary_conditions(f1);
            check(f);
            }
        if (OutStep!=0 && count%OutStep==0)
            {
            if (CheckStep!=0 && count%CheckStep!=0)
                {
                boundary_conditions(f1);
                check(f);
                }
            printing(f1,dtdid,t_cur,count,PulsEnergy);
            }
        if (SnapStep!=0 && count%SnapStep==0)
            snapshot(f,t_cur,count);
        for(l=0;l<nvar;l++)
        for(i=0;i<m1;i++)
        for(j=0;j<m2;j++)
        for(k=0;k<m3;k++)
           f[l][i][j][k]=f1[l][i][j][k];
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
   free_mem_2f(refr,m1, m3);
   free_mem_2f(refz,m1, m3);
   free_mem_2f(node,m1, m3);
   free_mem_3f(averf,3,m1,m3);

   if(t_cur>=Ttot&&!razlet) nmessage("work is succesfully done",t_cur);
       else nrerror("this is break of scheme",t_cur);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   nmessage("mpi_finalize is done",t_cur);
return 0;
}
