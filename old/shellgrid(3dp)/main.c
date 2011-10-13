#define LEVEL

#include "head.h"

int main(int argc, char** argv)
{
   double dttry, dtdid, dtnext;
   int i,j,k,l,i2,j2,k2;
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
   NameSnapFile = fname;
  NameStatFile =(char *)calloc(strl+9,sizeof(char));
  sprintf(NameStatFile,"%s_%d.sta",fname,rank); NameStatFile[strl+8] = 0;

   numlog = 1;
   Master nmessage("--------------------------------------------------------------------------",-1,-1);

/* ---------------------- reading files and arrays --------------------- */
   goon = ((fd=fopen(NameCPFile,"r+"))>0);
   if(fd==NULL) { putlog("File of cp were not opened",goon);
                  Master if((fd=fopen(NameCPFile,"w+"))!=NULL) putlog("File cp was successfully created",1);
                }
           else putlog("File of control points opened=",(long)fd);
	if(goon)
		{
		do fscanf(fd,"%s\n",NameInitFile); while (!feof(fd));
		goon = strcmp(NameInitFile,"END") || strlen(NameInitFile)==0;
		fileclose(fd);
		}
   init_param(argc,argv,&dtnext);       // initilization of parameters
   Gamma=1e-3;

   ghost=(approx-1)/2;            //radius of approx sample
   dx[0]=l1/N1;
   dx[1]=l2/N2;
   dx[2]=l3/N3;
   p1 = 8*l1/(l3*Re) ; p2 = 0;

   t_cur=0;
   count=0; enter = 0;

   if(goon) {if(init_data()) nrerror("error of reading initial arrays",-1,-1);}
       else { init_parallel();  operate_memory(1);}

   init_conditions();

//   init_shell();
//--------------------------------------
   boundary_conditions(f,nut);

   if(!goon)  dump(f,nut,t_cur,count);

   time_begin = MPI_Wtime();
   if(!goon) Master nmessage("work has begun",0,0);
       else Master nmessage("work continued",t_cur,count);
   Master ferror = fileopen(NameErrorFile,abs(goon));

   if(CheckStep!=0) check(f);
   if(!goon) if (OutStep!=0) printing(f,0,t_cur,count,PulsEnergy);   //in parallel version better not to do

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
   if(rank==0) add_control_point("END");

   Master fileclose(ferror);

   MPI_Barrier(MPI_COMM_WORLD);
   operate_memory(-1);
   if(t_cur>=Ttot&&!razlet) nmessage("work is succesfully done",t_cur,count);
       else nrerror("this is break of scheme",t_cur,count);
   MPI_Finalize();
   nmessage("mpi_finalize is done",t_cur,count);
return 0;
}
