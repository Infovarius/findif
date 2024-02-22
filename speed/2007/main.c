#define LEVEL

#include "head.h"
#include "time.h"
#define Tunit 0.36915

int main(int argc, char** argv)
{
   double dttry, dtdid, dtnext;
   double vrho,vphi,vth, rho, r1, phi1, z1;              // for given velocity profile
   int i,j,k,l,i2,j2,k2;
   FILE *fd;
   int strl;
   double ChangeParamTime = 3, DeltaParam = 10;         // for iteration on parameters
   double dv1[7][3],dv2[7][3],dA11[7][3][3],dp1[3],dn1[3],w,dw;
   int m;

 /* Initialize MPI */
 MPI_Init(&argc,&argv);
 /*  Get rank of this process and process group size */
 MPI_Comm_rank(MPI_COMM_WORLD,&rank);
 MPI_Comm_size(MPI_COMM_WORLD,&size);

  srand(rank);
//init_tick(0,"init");
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

   Master nmessage("--------------------------------------------------------------------------",-1,-1);
   init_param(argc,argv,&dtnext);       // initialization of parameters
   Gamma=1e-4;
   ghost=(approx-1)/2;                  //radius of approx sample
   t_cur=0;
   count=0; enter = 0;

/* ---------------------- initialization of arrays --------------------- */
   goon = ((fd=fopen(NameCPFile,"r+"))>0);
   if(fd==NULL) { //putlog("File of cp were not opened",goon);
                  Master if((fd=fopen(NameCPFile,"w+"))!=NULL) ;//putlog("File cp was successfully created",1);
                }
           else ;//putlog("File of control points opened=",(long)fd);
   if(goon)
      { do fscanf(fd,"%s\n",NameInitFile); while (!feof(fd));
        goon = strcmp(NameInitFile,"END");
      }
   if(goon) {if(init_data()) nrerror("error of reading initial arrays",-1,-1);}
       else { init_parallel();  operate_memory(1);}
   fileclose(fd);

   dx[0]=2*R/N1;
   dx[1]=lfi/N2;
   dx[2]=2*R/N3;
   init_conditions();
   init_timers();

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

 print_array2d(fd,node,0,m1,0,m3);
 print_array2d(fd,refr_f,0,m1,0,m3);
 print_array2d(fd,refz_f,0,m1,0,m3);
 print_array2d(fd,refr_m,0,m1,0,m3);
 print_array2d(fd,refz_m,0,m1,0,m3);
 fileclose(fd);

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

start_tick(10,"all_dr");
 for(count=0;count<10000;count++) {
      for(i=ghost;i<mm1;i++)
   for(j=ghost;j<mm2;j++)
   for(k=ghost;k<mm3;k++)
   {
     if(isType(node[i][k],NodeMagn) && !isType(node[i][k],NodeClued))
      {
//start_tick(15,"");
      for(l=4;l<=6;l++) {
       for(m=0;m<3;m++) {
//   start_tick(17,"");
         dv1[l][m]=dr(f[l],i,j,k,m+1,0,dx[m],ghost,approx);
//        dv1[l][m] = 1.2*(f[l][i+1][j][k]-f[l][i-1][j][k])+2.3*f[l][i][j][k];
/*        dv1[l][m] = 1.2*(nut[i+3][j][k]-nut[i-3][j][k])
                               + 2.3*(nut[i+2][j][k]-nut[i-2][j][k])
                               + 3.4*(nut[i+1][j][k]-nut[i-1][j][k])
                               + 4.5*nut[i][j][k];*/
//   finish_tick(17); start_tick(18,"");
//         dv2[l][m]=dr(f[l],i,j,k,m+1,1,dx[m]*dx[m],ghost, approx);
//finish_tick(18);
	 }     //	start_tick(19,"");
/*        dA11[l][0][1] = dA11[l][1][0] = d2cross(f[l],i,j,k,2,1,ghost,approx);
        dA11[l][0][2] = dA11[l][2][0] = d2cross(f[l],i,j,k,3,1,ghost,approx);
        dA11[l][1][2] = dA11[l][2][1] = d2cross(f[l],i,j,k,3,2,ghost,approx);*/
//finish_tick(19);
        }
//finish_tick(15);
      } //  else df[4][i][j][k] = df[5][i][j][k] = df[6][i][j][k] = 0;
  } //global for

   }
finish_tick(10);
   printing(f1,dtdid,t_cur,count,PulsEnergy);
   snapshot(f,eta,t_cur,count);
   if(rank==size-1) add_control_point("END");
   print_CPU_usage();

//   Master fileclose(NameErrorFile);

   if(t_cur>=Ttot&&!razlet) nmessage("work is succesfully done",t_cur,count);
       else nrerror("this is break of scheme",t_cur,count);
   MPI_Barrier(MPI_COMM_WORLD);
   operate_memory(-1);
   MPI_Finalize();
   nmessage("mpi_finalize is done",t_cur,count);
return 0;
}
