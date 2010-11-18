//------------------------ all the outputting stuff -----------------//

#define LEVEL extern
//#include <conio.h>
#include "head.h"

double UpLimit;     //after this limit there's dump

FILE *fileopen(const char *x, int mode)  //opening of file to ff
                         /*   0-rewrite;others-append   */
{
FILE *ff;
char* s_mode;
if(mode>0) s_mode="a";
if(mode<0) s_mode="r";
if(mode==0) s_mode="w";
if ((ff = fopen(x,s_mode))==NULL)
	 {
		printf ("Can't open file %s !\n",x);
		exit(-1);
	 }
return(ff);
}

void read_token(FILE *inp,double *param)
{
char str[256],*pstr;
 do fgets(str,256,inp); while(strchr(str,'|')==NULL);
 if(strchr(str,'|')==NULL) return;
 if(sscanf(str,"%lf",param)==0) nrerror("Input of parameter error",0);
 pstr=strtok(str,"|");
 pstr=strtok((char *)NULL,"|");
 printf("%g\t->\t%s",*param,pstr);
}

void init_param(int argc, char** argv,double *dtnext)
{
int error=0;
FILE *iop;
double d;
 if(argc<2 || (iop=fopen(argv[1],"r"))==NULL) //no ini file
    {
     if(argc>=2) nmessage("Start: no ini file!",0);
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
      }
}

void read_tilleq(FILE *ffff,char echo)
{char ch;
  if (echo=='n') while ((ch=(char)fgetc(ffff))!='=');
         else    while ((ch=(char)fgetc(ffff))!='=') printf("%c",ch);
}

int init_data(double *Re,double *t_cur,long *count)
                          //returns code of error
{
 int error=0;
 int i,j,k,l;
 float tmpd;
 char tmpc;
 int n1,n2,n3;

 FILE *inp = fileopen(NameInitFile,-1);
 read_tilleq(inp,'n');   if(fscanf(inp,"%lf",t_cur)==0) error=1;
 read_tilleq(inp,'n');   if(fscanf(inp,"%ld",count)==0) error=1;
 read_tilleq(inp,'n');   if(fscanf(inp,"%d",&n1)==0) error=1;
 read_tilleq(inp,'n');   if(fscanf(inp,"%d",&n2)==0) error=1;
 read_tilleq(inp,'n');   if(fscanf(inp,"%d",&n3)==0) error=1;
 read_tilleq(inp,'n');   if(fscanf(inp,"%lf",Re)==0) error=1;
// recreate_arrays(n1,n2,n3);
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
fclose(inp);
return(error);
}

void print_array1d(FILE *ff,double *a,int beg1,int n1)
{
int i;
fprintf(ff,"{");
for(i=beg1;i<beg1+n1;i++)
        {
        fprintf(ff,"%0.10f",a[i]);
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
    for(j=beg2;j<beg2+n2-1;j++)
        {
        fprintf(ff,"%0.10f,",a[i][j]);
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
            fprintf(ff,"%0.10f",a[i][j][k]);
            fprintf(ff,k<beg3+n3-1 ? "," : "}");
            }
        fprintf(ff,j<beg2+n2-1 ? "," : "}");
        }
    fprintf(ff,i<beg1+n1-1 ? "," : "}\n");
    }
}

double check(double ****f)   //give energy of pulsations
{
double PulsEnergy=0;
double **averf;
int i,j,k,l;
int n=max(n1,max(n2,n3));
TotalEnergy=0;
averf = alloc_mem_2f(3,n);
for(l=0;l<=2;l++)
   for(i=0;i<n;i++)
       averf[l][i] = 0;
for(i=ghost;i<mm1;i++)
     for(j=ghost;j<=ghost;j++)
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
free_mem_2f(averf,3,n);
return(PulsEnergy);
}

void printing(double ****f1,double dtdid,double t_cur,long count,double en)
{
double temp, div=0;
int i,j,k,l;
double mf, mda, mdr;
double *avervx, *avernu;
FILE *fv,*fnu,*fen,*fkv;

//clrscr();
boundary_conditions(f1);

for(i=ghost;i<mm1;i++)
   for(j=ghost;j<mm2;j++)
      for(k=ghost;k<mm3;k++) {
           temp=0;
           for(l=0;l<3;l++)
            temp+=dr(f1[l],i,j,k,l+1,0,dx[l],ghost, approx);
           if (fabs(temp)>div) div=fabs(temp);
           }
time_now = MPI_Wtime();
printf("program is working %0.2f seconds\n",time_now-time_begin);

printf("t=%e dtdid=%e NIter=%d %e\n", t_cur, dtdid, count, div);

   for(l=0;l<nvar;l++) {
       mf=mda=mdr=0;
       for(i=ghost;i<mm1;i++)
        for(j=ghost;j<mm2;j++)
         for(k=ghost;k<mm3;k++) {
            temp=fabs(f[l][i][j][k]-f1[l][i][j][k]);
            if (temp>mda) mda=temp;
            if(f1[l][i][j][k]!=0) temp/=f1[l][i][j][k];
            if (temp>mdr) mdr=temp;
            if (fabs(f1[l][i][j][k])>mf) mf=fabs(f1[l][i][j][k]);
          }
          printf("%d %e %e %e\n",l, mf, mda, mdr);
          }

         printf("Energy of pulsations=%g\n",en);
         fen = fileopen(NameEnergyFile,count);
         fprintf(fen,"%0.10f\n",en);
         fclose(fen);
         printf("number of runge-kutt calculations=%d",enter);

}

void dump(double ****f1,double t_cur,long count)
{
FILE *fd;
fd = fileopen(NameDumpFile,0);
nmessage("dump is done",t_cur);
fprintf(fd,"current time = %0.10f \ncurrent iteration = %ld\n",t_cur,count);
fprintf(fd,"Number of points along x = %d\n",n1);
fprintf(fd,"Number of points along y = %d\n",n2);
fprintf(fd,"Number of points along z = %d\n",n3);
fprintf(fd,"Reynolds number = %lf\n",Re);
print_array3d(fd,f1[0],0,m1,0,m2,0,m3);
print_array3d(fd,f1[1],0,m1,0,m2,0,m3);
print_array3d(fd,f1[2],0,m1,0,m2,0,m3);
print_array3d(fd,f1[3],0,m1,0,m2,0,m3);
print_array3d(fd,nut,0,m1,0,m2,0,m3);
fclose(fd);
}

void snapshot(double ****f,double t_cur,long count)
{
char str[256];
int i, j, l;
FILE *fd;

 sprintf(str,"%s_%d.snp",NameSnapFile,count);
 fd=fopen(str,"w");
nmessage("snap is done",t_cur);
fprintf(fd,"current time = %0.10f \ncurrent iteration = %ld\n",t_cur,count);
fprintf(fd,"Number of points along x = %d\n",n1);
fprintf(fd,"Number of points along y = %d\n",n2);
fprintf(fd,"Number of points along z = %d\n",n3);
fprintf(fd,"Reynolds number = %lf\n",Re);
print_array3d(fd,f1[0],0,m1,0,m2,0,m3);
print_array3d(fd,f1[1],0,m1,0,m2,0,m3);
print_array3d(fd,f1[2],0,m1,0,m2,0,m3);
print_array3d(fd,f1[3],0,m1,0,m2,0,m3);
print_array3d(fd,nut,0,m1,0,m2,0,m3);
 fclose(fd);

}