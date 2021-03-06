//------------------------ all the outputting stuff -----------------//

#define LEVEL extern
//#include <conio.h>
#include "head.h"

double UpLimit;     //after this limit there's dump
#define PREC 5

FILE *fileopen(const char *x, int mode)  //opening of file to ff
                         /*   0-rewrite;>0-append;<0-read    */
{
FILE *ff;
char* s_mode;
if(mode>0) s_mode="a";
if(mode<0) s_mode="r";
if(mode==0) s_mode="w";
if ((ff = fopen(x,s_mode))==NULL)
	 {
		nrerror ("Can't open file !\n",t_cur,count);
		exit(-1);
	 }
return(ff);
}

void fileclose(FILE *fd)        // closing file with checking
{
if(fd==NULL) putlog("It was endeavor to close unopened file",1);
        else fclose(fd);
}

void putlog(char msg_text[],long num)
{
   FILE *log;
   char name[20];
   sprintf(name,"mess%d.log",rank);
   log=fileopen(name,num);
   time_now = (num==0)?time_begin:MPI_Wtime();
   fprintf(log,"message at t=%-7.4lf Niter=%-6d time of work=%g sec\n",
                    t_cur,count,time_now-time_begin);
   fprintf(log,"%s %d\n",msg_text,num);
   fileclose(log);
}

void read_token(FILE *inp,double *param)
{
char str[256],*pstr;
 do fgets(str,256,inp); while(strchr(str,'|')==NULL);
 if(strchr(str,'|')==NULL) return;
 if(sscanf(str,"%lf",param)==0) nrerror("Input of parameter error",-1,-1);
 pstr=strtok(str,"|");
 pstr=strtok((char *)NULL,"|");
 printf("%g\t->\t%s",*param,pstr);
}

void init_param(int argc, char** argv,double *dtnext)
{
int ver;
FILE *iop;
double d;
 if(argc<2 || (iop=fopen(argv[1],"r"))==NULL) //no ini file
    {
     nrerror("Start from no ini file!",0,0);
     Re=10.;
     l1=3.;
     l2=1.;
     l3=1.;
     parabole=0.;
     Noise=0.;
     NoiseNorm=0.;
     UpLimit=10.;
     N1=10;
     N2=10;
     N3=10;
     nvar=4;
     approx=7;                    //derivatives approximation order
     *dtnext=1e-3;
     Ns=15;
     maschtab=1e6;
     lambda = 2.0;
     max_okr = 3;
     OutStep = (CheckStep=100)/1;
     VarStep = 0;
     SnapStep = 100;
     Ttot=1.;
     }
    else {       
      if(fscanf(iop,"%d",&ver)<1 || ver>2) nrerror("parameters' file has wrong version",0,0);
      read_token(iop,&Re);
      read_token(iop,&l1);
      read_token(iop,&l2);
      read_token(iop,&l3);
      read_token(iop,&parabole);
      read_token(iop,&Noise);
      read_token(iop,&NoiseNorm);
      read_token(iop,&UpLimit);
      read_token(iop,&d);         N1 = (int)d;
      read_token(iop,&d);         N2 = (int)d;
      read_token(iop,&d);         N3 = (int)d;
      read_token(iop,&d);         nvar = (int)d;
      read_token(iop,&d);         approx = (int)d;
      read_token(iop,dtnext);
      read_token(iop,&d);         Ns = (int)d;
      read_token(iop,&maschtab);
      read_token(iop,&lambda);
      read_token(iop,&d);         max_okr = (int)d;
      read_token(iop,&d);         OutStep = (int)d;
      read_token(iop,&d);         SnapStep = (int)d;
      read_token(iop,&d);         CheckStep = (int)d;
      read_token(iop,&d);         VarStep = (int)d;
      read_token(iop,&Ttot);
      fclose(iop);
      nmessage("Parameters were extracted from file",0,0);
      }
}

void read_tilleq(FILE *ffff,char echo)
{char ch;
  if (echo=='n') while ((ch=(char)fgetc(ffff))!='=');
         else    while ((ch=(char)fgetc(ffff))!='=') printf("%c",ch);
}

int init_data()                 //returns code of error
{
 int error=0;
 int i,j,k,l,tmpr;
 float tmpd;
 char tmpc;

 FILE *inp = fileopen(NameInitFile,-1);
 read_tilleq(inp,'n');   if(fscanf(inp,"%lf",&t_cur)==0) error=1;
 read_tilleq(inp,'n');   if(fscanf(inp,"%ld",&count)==0) error=1;
 read_tilleq(inp,'n');   if(fscanf(inp,"%c%d%c%d%c%d%c",&tmpc,&pp[0],&tmpc,&pp[1],&tmpc,&pp[2],&tmpc)<7) error=1;
                         //no need unless process distribution is written
 read_tilleq(inp,'n');   if(fscanf(inp,"%d",&N1)==0) error=1;
 read_tilleq(inp,'n');   if(fscanf(inp,"%d",&N2)==0) error=1;
 read_tilleq(inp,'n');   if(fscanf(inp,"%d",&N3)==0) error=1;
 read_tilleq(inp,'n');   if(fscanf(inp,"%lf",&Re)==0) error=1;

 init_parallel();
 operate_memory(1);                     // creating arrays

 for(tmpr=0;tmpr<=rank;tmpr++)           //reading until arrays of this process
 {
 for(l=0;l<nvar;l++)                    // reading f
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
                  if(fscanf(inp,"%g",&tmpd)==0) error=2;
                  fscanf(inp,"%c",&tmpc);
                  f[l][i][j][k]=tmpd;
                  }
               fscanf(inp,"%c",&tmpc);
               }
            fscanf(inp,"%c",&tmpc);
            }
        fscanf(inp,"%c",&tmpc);
        }
 do fscanf(inp,"%c",&tmpc); while (tmpc!='{');   //reading nut
 for(i=0;i<m1;i++)
     {
     do fscanf(inp,"%c",&tmpc); while (tmpc!='{');
     for(j=0;j<m2;j++)
         {
         do fscanf(inp,"%c",&tmpc); while (tmpc!='{');
         for(k=0;k<m3;k++)
             {
             if(fscanf(inp,"%g",&tmpd)==0) error=3;
             fscanf(inp,"%c",&tmpc);
             nut[i][j][k]=tmpd;
             }
         fscanf(inp,"%c",&tmpc);
         }
     fscanf(inp,"%c",&tmpc);
     }
 fscanf(inp,"%c",&tmpc);
 }
fileclose(inp);
if(error) nrerror("Data couldn't have been read from file!!!",0,0);
     else nmessage("Data has been read from file",t_cur,count);
return(error);
}

void print_array1d(FILE *ff,double *a,int beg1,int n1)
{
int i;
fprintf(ff,"{");
for(i=beg1;i<beg1+n1;i++)
        {
        fprintf(ff,"%0.*G",PREC,a[i]);
        fprintf(ff,i<beg1+n1-1 ? "," : "}\n");
        }
}

void print_array2d(FILE *ff,double **a,int beg1,int n1,int beg2,int n2)
{
int i,j;
fprintf(ff,"{");
for(i=beg1;i<beg1+n1;i++)
    {
    fprintf(ff,"{");
    for(j=beg2;j<beg2+n2;j++)
        {
        fprintf(ff,"%0.*G",PREC,a[i][j]);
        fprintf(ff,j<beg2+n2-1 ? "," : "}");
        }
    fprintf(ff,i<beg1+n1-1 ? "," : "}\n");
    }
}

void print_array3d(FILE *ff,double ***a,
        int beg1,int n1,int beg2,int n2,int beg3,int n3)
{
int i,j,k;
fprintf(ff,"{");
for(i=beg1;i<beg1+n1;i++) {
    fprintf(ff,"{");
    for(j=beg2;j<beg2+n2;j++) {
        fprintf(ff,"{");
        for(k=beg3;k<beg3+n3;k++) {
            fprintf(ff,"%0.*G",PREC,a[i][j][k]);
            fprintf(ff,k<beg3+n3-1 ? "," : "}");
            }
        fprintf(ff,j<beg2+n2-1 ? "," : "}");
        }
    fprintf(ff,i<beg1+n1-1 ? "," : "}\n");
    }
}

void check(double ****f)   //calculate energy of pulsations of all components and know if there's crash
{
int i,j,k,l;
PulsEnergy=0;
TotalEnergy=0;
for(l=0;l<=2;l++)
   for(i=0;i<m2;i++)
       averf[l][i] = 0;
for(i=ghost;i<mm1;i++)
     for(j=ghost;j<mm2;j++)
        for(k=ghost;k<mm3;k++)
           {
           PulsEnergy+=deviation(f,i,j,k);
           for(l=0;l<=2;l++) averf[l][k-ghost] += f[l][i][j][k];
           }
for(l=0;l<=2;l++)
   for(k=ghost;k<mm3;k++)
        {
//        averf[l][k-ghost] /= n1*n2;
//        for(i=ghost;i<mm1;i++)
//            for(j=ghost;j<=ghost;j++)
               TotalEnergy += pow(averf[l][k-ghost],2.);
        }
TotalEnergy += 1.;   //if zero average field
razlet = (PulsEnergy/TotalEnergy>UpLimit);
}

void printing(double ****f1,double dtdid,double t_cur,long count,double en)
{
double temp, div, totdiv;
int i,j,k,l;
double mf[3], totmf[3], toten;     //mf[0]=max(f), mf[1]=max(df), mf[2]=max(df/f)
FILE *fv,*fnu,*fen,*fkv;

//clrscr();
time_now = MPI_Wtime();
Master printf("program is working %0.2f seconds\n",time_now-time_begin);

for(i=ghost;i<mm1;i++)
   for(j=ghost;j<mm2;j++)
      for(k=ghost;k<mm3;k++) {
           temp=0;
           for(l=0;l<3;l++)
            temp+=dr(f1[l],i,j,k,l+1,0,dx[l],ghost, approx);
           if (fabs(temp)>div) div=fabs(temp);
           }
MPI_Allreduce(&div, &totdiv, 1, MPI_DOUBLE , MPI_MAX, MPI_COMM_WORLD);
Master printf("t=%g dtdid=%g NIter=%d maxdiv=%g(local=%g)\n", t_cur, dtdid, count, totdiv, div);

   for(l=0;l<nvar;l++) {
       mf[0]=mf[1]=mf[2]=0;
       for(i=ghost;i<mm1;i++)
        for(j=ghost;j<mm2;j++)
         for(k=ghost;k<mm3;k++) {
            if (fabs(f1[l][i][j][k])>mf[0]) mf[0]=fabs(f1[l][i][j][k]);
            temp=fabs(f[l][i][j][k]-f1[l][i][j][k]);
            if (temp>mf[1]) mf[1]=temp;
            if(f1[l][i][j][k]!=0) temp/=f1[l][i][j][k];
            if (temp>mf[2]) mf[2]=temp;
          }
          MPI_Allreduce(mf, totmf, 3, MPI_DOUBLE , MPI_MAX, MPI_COMM_WORLD);
          Master printf("%d  maxf=%e(loc=%e) \tmaxdf=%e(loc=%e) \tmax(df/f)=%e(loc=%e)\n",
                          l,      totmf[0],mf[0],    totmf[1],mf[1],       totmf[2],mf[2]);
          }

         MPI_Allreduce(&en, &toten, 1, MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
         Master printf("Energy of pulsations=%g (local=%g)\n",toten,en);
         Master {fen = fileopen(NameEnergyFile,count);
                 fprintf(fen,"%8.8g  \t%e\n",t_cur,toten);
                 fclose(fen);
                }
         Master printf("number of runge-kutt calculations=%d\n",enter);

 // -------------------- average profile of velocity ---------------------
         for(i=0;i<N3;i++)    vx[i]=0;
         for(i=0;i<N3;i++) totvx[i]=0;
         for(i=ghost;i<mm1;i++)
            for(j=ghost;j<mm2;j++)
               for(k=ghost;k<mm3;k++)
                   vx[k-ghost+n[2]] += f[0][i][j][k];
         MPI_Allreduce(vx, totvx, N3, MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
         Master {for(i=0;i<N3;i++) totvx[i] /= N1*N2;
                 fv = fileopen(NameVFile,count);
                 fprintf(fv,"{%8.8f}\t",t_cur);
                 print_array1d(fv,totvx,0,N3);
                 fclose(fv);
                }
}

void dump(double ****f1,double ***nu,double t_cur,long count)
{
FILE *fd;
char *message;
int tag=1;

 if(rank!=0) MPI_Recv(message,0,MPI_CHAR,rank-1,tag,MPI_COMM_WORLD,statuses);

 fd=fileopen(NameDumpFile,rank);

 Master nmessage("dump has been started",t_cur,count);
 Master fprintf(fd,"current time = %0.10f \ncurrent iteration = %ld\n",t_cur,count);
 Master fprintf(fd,"number of processors along axes={%d,%d,%d}\n",pp[0],pp[1],pp[2]);
 Master fprintf(fd,"Number of points along x = %d\n",N1);
 Master fprintf(fd,"Number of points along y = %d\n",N2);
 Master fprintf(fd,"Number of points along z = %d\n",N3);
 Master fprintf(fd,"Reynolds number = %lf\n",Re);

 print_array3d(fd,f1[0],0,m1,0,m2,0,m3);
 print_array3d(fd,f1[1],0,m1,0,m2,0,m3);
 print_array3d(fd,f1[2],0,m1,0,m2,0,m3);
 print_array3d(fd,f1[3],0,m1,0,m2,0,m3);
 print_array3d(fd,nut,0,m1,0,m2,0,m3);
 fclose(fd);

 if(rank!=size-1) MPI_Send(message,0,MPI_CHAR,rank+1,tag,MPI_COMM_WORLD);
             else nmessage("dump is done",t_cur,count);
 add_control_point(NameDumpFile);
}

void snapshot(double ****f1,double ***nu,double t_cur,long count)
{
char str[256];
char message[10]="message";
long tag=count;
int i, j, l;
FILE *fd;

 sprintf(str,"%s_%d.snp",NameSnapFile,count);

 if(rank!=0) MPI_Recv(message,1,MPI_CHAR,rank-1,tag,MPI_COMM_WORLD,statuses);


 fd=fileopen(str,rank);
 Master nmessage("snap has been started",t_cur,count);
 Master fprintf(fd,"current time = %0.10f \ncurrent iteration = %ld\n",t_cur,count);
 Master fprintf(fd,"number of processors along axes={%d,%d,%d}\n",pp[0],pp[1],pp[2]);
 Master fprintf(fd,"Number of points along x = %d\n",N1);
 Master fprintf(fd,"Number of points along y = %d\n",N2);
 Master fprintf(fd,"Number of points along z = %d\n",N3);
 Master fprintf(fd,"Reynolds number = %lf\n",Re);

 print_array3d(fd,f1[0],0,m1,0,m2,0,m3);
 print_array3d(fd,f1[1],0,m1,0,m2,0,m3);
 print_array3d(fd,f1[2],0,m1,0,m2,0,m3);
 print_array3d(fd,f1[3],0,m1,0,m2,0,m3);
 print_array3d(fd,nut,0,m1,0,m2,0,m3);
 fclose(fd);
                  
 if(rank!=size-1) MPI_Send(message,1,MPI_CHAR,rank+1,tag,MPI_COMM_WORLD);
             else nmessage("snap is done",t_cur,count);
 add_control_point(str);            
}
