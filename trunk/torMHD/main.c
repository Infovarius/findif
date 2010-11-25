#define LEVEL

#include "head.h"
#include "time.h"
#define Tunit 0.36915

int main(int argc, char** argv)
{
   double dttry, dtdid, dtnext, tmp,tmpT;
   int i,j,k,l,tmpC;
   int outed;
   FILE *fd, *ferror;
   int strl;
   double ****ft;

 /* Initialize MPI */
 MPI_Init(&argc,&argv);
 /*  Get rank of this process and process group size */
 MPI_Comm_rank(MPI_COMM_WORLD,&rank);
 MPI_Comm_size(MPI_COMM_WORLD,&size);

  srand(rank);
  NameMessageFile = "message.dat";
  NameErrorFile = "error.err";
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

   Master nmessage("--------------------------------------------------------------------------",-1,-1);

/* ---------------------- reading files and arrays --------------------- */
   goon = ((fd=fopen(NameCPFile,"r+"))>0);
   if(fd==NULL) { //putlog("File of cp were not opened",goon);
                  Master if((fd=fopen(NameCPFile,"w+"))!=NULL) ;//putlog("File cp was successfully created",1);
                }
           else ;//putlog("File of control points opened=",(long)fd);
   if(goon)
      { do fscanf(fd,"%s\n",NameInitFile); while (!feof(fd));
        goon = strcmp(NameInitFile,"END");
      }

   init_param(argc,argv,&dtnext);       // initialization of parameters
   Gamma=1e-4;
   ghost=(approx-1)/2;                  //radius of approx sample
   t_cur=0;
   count=0; enter = 0;
   timeE1=0;

   if(goon) {if(init_data()) nrerror("error of reading initial arrays",-1,-1);}
       else { init_parallel();  operate_memory(1);}
   fileclose(fd);

   dx[0]=2*R/N1;
   dx[1]=lfi/N2;
   dx[2]=2*R/N3;
   p1 = 4/Re;

   init_conditions();

//--------------------------------------
  Master {
      fd=fileopen("coord",rank);
      for(i=0;i<N1+2*ghost;i++) fprintf(fd,"%e ",coordin(i,0));
      fprintf(fd,"\n");
      for(i=0;i<N2+2*ghost;i++) fprintf(fd,"%e ",coordin(i,1));
      fprintf(fd,"\n");
      for(i=0;i<N3+2*ghost;i++) fprintf(fd,"%e ",coordin(i,2));
      fileclose(fd);
     }
//--------------------------------------
 {
    if(rank!=0) MPI_Recv(fname,0,MPI_CHAR,rank-1,1,MPI_COMM_WORLD,statuses);

 fd=fileopen("node",rank);

 Master nmessage("nodes outputting has been started",t_cur,count);

 print_array2i(fd,node,0,m1,0,m3);
 print_array2d(fd,refr_f,0,m1,0,m3);
 print_array2d(fd,refz_f,0,m1,0,m3);
 print_array2d(fd,refr_m,0,m1,0,m3);
 print_array2d(fd,refz_m,0,m1,0,m3);
 fileclose(fd);

 if(rank!=size-1) MPI_Send(fname,0,MPI_CHAR,rank+1,1,MPI_COMM_WORLD);
             else nmessage("nodes has been dumped",t_cur,count);
            }
//--------------------------------------

   boundary_conditions(f,nut);

//   if(!goon) dump(f,eta,t_cur,count);
   time_begin = MPI_Wtime();
   if(!goon) Master nmessage("work has begun",0,0);
       else Master nmessage("work continued",t_cur,count);
   Master ferror = fileopen(NameErrorFile,abs(goon));

   if(CheckStep!=0) check(f);
   if (OutStep!=0) printing(f,0,t_cur,count,PulsEnergy);

/*------------------------ MAIN ITERATIONS -------------------------*/
   while ((Ttot==0 || t_cur < Ttot) && !razlet) {
        pde(t_cur, f, df);
        dttry=dtnext;
        timestep(f, df, t_cur, f1, dttry, &dtdid, &dtnext);
  //      nut_by_flux(f,dtdid);
        t_cur+=dtdid;
        count++;
        if(Ttot!=0 && t_cur >= Ttot) break;
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
                };
//              else boundary_conditions(f1,nut);
            time1=MPI_Wtime();
            printing(f1,dtdid,t_cur,count,PulsEnergy);
            time0=MPI_Wtime();
            timeE1=0;
            }
        outed = 0;
        if (SnapStep!=0 && count%SnapStep==0)
            { snapshot(f1,eta,t_cur,count); outed=1; }
        if (SnapDelta>5*dtdid && floor((t_cur-dtdid)/SnapDelta)<floor(t_cur/SnapDelta) && !outed)
            { snapshot(f1,eta,t_cur,count); outed=1; }
        ft = f;  f = f1;  f1 = ft;
        if (count%100==0) {
//          MPI_Barrier(MPI_COMM_WORLD);
          tmp=Rm;
          init_param(argc,argv,&dttry);
          Rm=tmp;
//          MPI_Barrier(MPI_COMM_WORLD);
        }

        if (ChangeParamTime!=0 && floor((t_cur-dtdid)/ChangeParamTime)<floor(t_cur/ChangeParamTime))
		{
		//Rm = floor(t_cur/ChangeParamTime+0.5)*DeltaParam-190;
		if(!outed) { snapshot(f,eta,t_cur,count); outed = 1;}
//                MPI_Barrier(MPI_COMM_WORLD);
		tmp = (Rm += DeltaParam);
                tmpC = count;  tmpT = t_cur;
		goon = ((fd=fopen(NameCPFile,"r+"))>0);
		if(goon)
			{ do fscanf(fd,"%s\n",NameInitFile); while (!feof(fd));
			goon = strcmp(NameInitFile,"END");
			}
		init_param(argc,argv,&dtnext);       // initialization of parameters
		ghost=(approx-1)/2;                  //radius of approx sample
		dx[0]=2*R/N1;
		dx[1]=lfi/N2;
		dx[2]=2*R/N3;

		fileclose(fd);
		if(goon) {if(init_data()) nrerror("error of reading initial arrays",-1,-1);}

                count = tmpC;  t_cur = tmpT;  Rm = tmp;
                p1 = 4/Re;
		init_conditions();
		goon = 1;
		Master nmessage("Rm was changed to",Rm,count);
   	   	}
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
   } // end while

   printing(f1,dtdid,t_cur,count,PulsEnergy);
   if(!outed) snapshot(f,eta,t_cur,count);
   if(rank==0) add_control_point("END");

   Master fileclose(ferror);

   if(t_cur>=Ttot&&!razlet) nmessage("work is succesfully done",t_cur,count);
       else nrerror("this is break of scheme",t_cur,count);
   MPI_Barrier(MPI_COMM_WORLD);
   operate_memory(-1);
   MPI_Finalize();
   nmessage("mpi_finalize is done",t_cur,count);
return 0;
}
