#define LEVEL

#include "head.h"

int main(int argc, char** argv)
{
   double dttry, dtdid, dtnext;
   int i,j,k,l,i2,j2,k2;
   FILE *fd;
   int strl;
   double ChangeParamTime = 3, DeltaParam = 10;         // for iteration on parameters

 /* Initialize MPI */
 MPI_Init(&argc,&argv);
 /*  Get rank of this process and process group size */
 MPI_Comm_rank(MPI_COMM_WORLD,&rank);
 MPI_Comm_size(MPI_COMM_WORLD,&size);

// MPI_Barrier(MPI_COMM_WORLD);

  srand(rank);
  NameMessageFile = "message.dat";
  NameErrorFile = "error.err";
  NameVFile  = "vv.dat";
  NameEnergyFile = "energy.dat";
  MainFlowFile = "torstab.main";
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
   init_param(argc,argv,&dtnext);       // initialization of parameters
   Gamma=1e-4;
   ghost=(approx-1)/2;                  //radius of approx sample

/* ---------------------- initialization of arrays --------------------- */

   init_parallel();
   operate_memory(1);                     // creating arrays

   goon = ((fd=fopen(NameCPFile,"r+"))>0);
   if(fd==NULL) { //putlog("File of cp were not opened",goon);
                  Master if((fd=fopen(NameCPFile,"w+"))!=NULL) ;//putlog("File cp was successfully created",1);
                }
           else ;//putlog("File of control points opened=",(long)fd);
   if(goon)
      { do fscanf(fd,"%s\n",NameInitFile); while (!feof(fd));
        goon = strcmp(NameInitFile,"END");
      }
   if(init_data(MainFlowFile,flow)) nrerror("error of reading main flow",-1,-1);
   t_cur=0;
   count=0; enter = 0;
   precalc();
   if(goon) {if(init_data(NameInitFile,f)) nrerror("error of reading initial arrays",-1,-1);}
   fileclose(fd);

   dx[0]=2*R/N1;
   dx[1]=lfi/N2;
   dx[2]=2*R/N3;

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
 if(!goon) {
    if(rank!=0) MPI_Recv("",0,MPI_CHAR,rank-1,1,MPI_COMM_WORLD,statuses);

 fd=fileopen("node",rank);

 Master nmessage("nodes outputting has been started",t_cur,count);

 print_array2i(fd,node,0,m1,0,m3);
 print_array2d(fd,refr,0,m1,0,m3);
 print_array2d(fd,refz,0,m1,0,m3);
 fileclose(fd);

 if(rank!=size-1) MPI_Send("",0,MPI_CHAR,rank+1,1,MPI_COMM_WORLD);
             else nmessage("nodes has been dumped",t_cur,count);
            }
//--------------------------------------

//   boundary_conditions(f,nut);
   
   if(!goon)  dump(f,nut,t_cur,count);

   time_begin = MPI_Wtime();
   if(!goon) Master nmessage("work has begun",0,0);
       else Master nmessage("work continued",t_cur,count);
   Master fileopen(NameErrorFile,abs(goon));

   if(CheckStep!=0) check(f);
   if(!goon) if (OutStep!=0) printing(f,0,t_cur,count,PulsEnergy);

/*------------------------ MAIN ITERATIONS -------------------------*/
   while (t_cur < Ttot && !razlet) {
	pde(t_cur, f, df);
	dttry=dtnext;
	timestep(f, df, nut, t_cur, f1, dttry, &dtdid, &dtnext);
	nut_by_flux(f,nut,dtdid);
	t_cur+=dtdid;
	count++;
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
	    snapshot(f1,nut,t_cur,count);
        if (SnapDelta>5*dtdid && floor((t_cur-dtdid)/SnapDelta)<floor(t_cur/SnapDelta))
            snapshot(f1,nut,t_cur,count);
/*        if (floor(t_cur/ChangeParamTime)<floor((t_cur+dtdid)/ChangeParamTime))
            {
            Re = floor(t_cur/ChangeParamTime+0.5)*DeltaParam-190;
         //Re -= DeltaParam;
            Master nmessage("Rm was changed to",Rm,count);
            }*/
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

   printing(f1,dtdid,t_cur,count,PulsEnergy);
   snapshot(f,nut,t_cur,count);
   if(rank==size-1) add_control_point("END");

//   Master fileclose(NameErrorFile);

   if(t_cur>=Ttot&&!razlet) nmessage("work is succesfully done",t_cur,count);
       else nrerror("this is break of scheme",t_cur,count);
   MPI_Barrier(MPI_COMM_WORLD);
   operate_memory(-1);
   MPI_Finalize();
   nmessage("mpi_finalize is done",t_cur,count);
return 0;
}
