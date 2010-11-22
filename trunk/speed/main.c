#define LEVEL

#include "head.h"

int main(int argc, char** argv)
{
   double dttry, dtdid, dtnext;
   double vrho,vphi,vth, rho, r1, phi1, z1;              // for given velocity profile
   int i,j,k,l,i2,j2,k2;
   FILE *fd;
   int strl;

 /* Initialize MPI */
 MPI_Init(&argc,&argv);
 /*  Get rank of this process and process group size */
 MPI_Comm_rank(MPI_COMM_WORLD,&rank);
 MPI_Comm_size(MPI_COMM_WORLD,&size);

start_tick(0,"init");
  NameMessageFile = "message.dat";
  NameErrorFile = "error.err";
  NameNuFile = "nut.dat";
  NameVFile  = "vv.dat";
  NameEnergyFile = "energy.dat";
  KaskadVarFile = "kaskvar.dat";
  fname = (argc>1)? argv[1] : argv[0];
  NameSnapFile = fname;
  strl=strlen(fname);
  NameCPFile = (char *)calloc(strl+4,sizeof(char));
  sprintf(NameCPFile,"%s.cp",fname);   NameCPFile[strl+3] = 0;
  NameDumpFile =(char *)calloc(strl+5,sizeof(char));
  sprintf(NameDumpFile ,"%s.dmp",fname); NameDumpFile[strl+4] = 0;
  NameStatFile =(char *)calloc(strl+9,sizeof(char));
  sprintf(NameStatFile,"%s_%d.sta",fname,rank); NameStatFile[strl+8] = 0;

   Master nmessage("--------------------------------------------------------------------------",0,0);
   init_param(argc,argv,&dtnext);       // initialization of parameters
   Gamma=1e-4;
   ghost=(approx-1)/2;                  //radius of approx sample
   MPI_Barrier(MPI_COMM_WORLD);
   t_cur=0;
   count=0; enter = 0;

/* ---------------------- initialization of arrays --------------------- */
   goon = ((fd=fopen(NameCPFile,"r"))>0);
   if(goon)
      { do fscanf(fd,"%s\n",NameInitFile); while (!feof(fd));
        goon = strcmp(NameInitFile,"END");
      }
   if(goon) {if(init_data()) nrerror("error of reading initial arrays",0,0);}
       else { init_parallel();  operate_memory(1);}
   fclose(fd);

   dx[0]=2*R/N1;
   dx[1]=lfi/N2;
   dx[2]=2*R/N3;
   init_conditions();
   init_timers();

//--------------------------------------
  if(!goon) Master {
      fd=fileopen("coord",rank);
      for(i=0;i<N1+2*ghost;i++) fprintf(fd,"%e ",coordin(i,0));
      fprintf(fd,"\n");
      for(i=0;i<N2+2*ghost;i++) fprintf(fd,"%e ",coordin(i,1));
      fprintf(fd,"\n");
      for(i=0;i<N3+2*ghost;i++) fprintf(fd,"%e ",coordin(i,2));
      fclose(fd);
     }
//--------------------------------------
 if(!goon) {
    if(rank!=0) MPI_Recv("",0,MPI_CHAR,rank-1,1,MPI_COMM_WORLD,statuses);

 fd=fileopen("node",rank);

 Master nmessage("nodes outputting has been started",t_cur,count);

 print_array2d(fd,node,0,m1,0,m3);
 print_array2d(fd,refr_f,0,m1,0,m3);
 print_array2d(fd,refz_f,0,m1,0,m3);
 print_array2d(fd,refr_m,0,m1,0,m3);
 print_array2d(fd,refz_m,0,m1,0,m3);
 fclose(fd);

 if(rank!=size-1) MPI_Send("",0,MPI_CHAR,rank+1,1,MPI_COMM_WORLD);
             else nmessage("nodes has been dumped",t_cur,count);
            }
//--------------------------------------

   boundary_conditions(f);

   if(!goon)  dump(f,eta,t_cur,count);

   time_begin = MPI_Wtime();
   if(!goon) Master nmessage("work has begun",0,0);
       else Master nmessage("work continued",t_cur,count);
   Master fileopen(NameErrorFile,abs(goon));

   if(CheckStep!=0) check(f);
   if(!goon) if (OutStep!=0) printing(f,0,t_cur,count,PulsEnergy);
finish_tick(0);

/*------------------------ MAIN ITERATIONS -------------------------*/
   while (t_cur < Ttot && !razlet) {
 start_tick(10,"");
       for(i=0;i<m1;i++)
         for(j=0;j<m2;j++)
            for(k=0;k<m3;k++)
            if(isType(node[i][k],NodeFluid))
            {
            r1 = coordin(i,0);  /*phi1 = coordin(j,1);*/  z1 = coordin(k,2);
            rho=sqrt(pow(r1-rc,2) + z1*z1);
            vrho = 0;
            if(rho>=Rfl) continue;
            vth  = vtheta_given(t_cur,rho,Rfl,phi1);
            vphi = vfi_given(t_cur,rho,Rfl);
/*            vth = 0;
            vphi = vfi_given(0.2,rho,Rfl);      */
            f[1][i][j][k] = vrho*sinth[i][k]+vth*costh[i][k];
            f[2][i][j][k] = vphi;
            f[3][i][j][k] = vrho*costh[i][k]-vth*sinth[i][k];
            }
    finish_tick(10);   start_tick(11,"");
        start_tick(7,"outpde");  pde(t_cur, f, df); finish_tick(7);
        dttry=dtnext;
    finish_tick(11);   start_tick(12,"");
        timestep(f, df, t_cur, f1, dttry, &dtdid, &dtnext);
    finish_tick(12);   start_tick(13,"");
        nut_by_flux(f,dtdid);
        t_cur+=dtdid;
        count++;
        if (CheckStep!=0 && count%CheckStep==0)
            {
            boundary_conditions(f1);
            check(f1);
            }
        if (OutStep!=0 && count%OutStep==0)
            {
            if (CheckStep!=0 && count%CheckStep!=0)
                {
                boundary_conditions(f1);
                check(f1);
                }
              else boundary_conditions(f1);
            printing(f1,dtdid,t_cur,count,PulsEnergy);
            }
        if (SnapStep!=0 && count%SnapStep==0)
            snapshot(f1,eta,t_cur,count);
    finish_tick(13);   start_tick(14,"exch");
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
    finish_tick(14); 
   }

   printing(f1,dtdid,t_cur,count,PulsEnergy);
   snapshot(f,eta,t_cur,count);
   if(rank==size-1) add_control_point("END");
   print_CPU_usage();

   operate_memory(-1);
//   Master fclose(NameErrorFile);

   if(t_cur>=Ttot&&!razlet) nmessage("work is succesfully done",t_cur,count);
       else nrerror("this is break of scheme",t_cur,count);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   nmessage("mpi_finalize is done",t_cur,count);
return 0;
}
