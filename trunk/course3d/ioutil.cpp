//------------------------ all the outputting stuff -----------------//

#define LEVEL extern
#include "head.h"
#include "kaskad.h"

const double UpLimit=10.;     //after this limit there's dump

FILE *fileopen(AnsiString x, int mode)  //opening of file to ff
                         /*   0-rewrite;>0-append;<0-read    */
{
FILE *ff;
char* s_mode;
if(mode>0) s_mode="a";
if(mode<0) s_mode="r";
if(mode==0) s_mode="w";
if ((ff = fopen(x.c_str(),s_mode))==NULL)
	 {
		nrerror ("Can't open file %s !\n" /*,x*/ ,t_cur);
		exit(-1);
	 }
return(ff);
}

void read_tilleq(FILE *ffff,char echo)
{char ch;
  if (echo=='n') while ((ch=char(fgetc(ffff)))!='=');
         else    while ((ch=char(fgetc(ffff)))!='=') printf("%c",ch);
}

int init_data(int *n1,int *n2,int *n3,double *Re,double ****f,double ***nut,double *t_cur,long *count)
                          //returns code of error
{
 int error=0;
 int i,j,k,l;
 float tmpd;
 char tmpc;

 FILE *inp = fileopen(NameInitFile,-1);
 read_tilleq(inp,'n');   if(fscanf(inp,"%lf",t_cur)==0) error=1;
 read_tilleq(inp,'n');   if(fscanf(inp,"%ld",count)==0) error=1;
 read_tilleq(inp,'n');   if(fscanf(inp,"%d",n1)==0) error=1;
 read_tilleq(inp,'n');   if(fscanf(inp,"%d",n2)==0) error=1;
 read_tilleq(inp,'n');   if(fscanf(inp,"%d",n3)==0) error=1;
 read_tilleq(inp,'n');   if(fscanf(inp,"%lf",Re)==0) error=1;
      MainWindow->EditRe->Text.printf("%0.10f",*Re);  //don't work
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
        fprintf(ff,"%0.10f",a[i][j]);
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
int n=max(n1,max(n2,n3));
TotalEnergy=0;
averf = alloc_mem_2f(3,n);
int i,j,k,l;
for(l=0;l<=2;l++)
   for(i=0;i<n;i++)
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
        averf[l][k-ghost] /= n1*n2;
        for(i=ghost;i<mm1;i++)
            for(j=ghost;j<mm2;j++)
               TotalEnergy += pow(averf[l][k-ghost],2.);
        }
TotalEnergy += 1.;   //if zero average field
razlet = (PulsEnergy/TotalEnergy>UpLimit/*||TotalEnergy<PulsEnergy*/);
free_mem_2f(averf,3,n);
return(PulsEnergy);
}

void printing(double ****f1,double dtdid,double t_cur,long count,double en)
{
double temp, div=0;
int i,j,k,l;
double mf, mda, mdr;
double *avervx, *avernu;
AnsiString TempScreen;
FILE *fv,*fnu,*fen,*fkv;
//fnu = fileopen(NameNuFile,count);
if(MainWindow->CheckClear->Checked)
         MainWindow->Screen->Lines->Clear();
boundary_conditions(f1);

for(i=ghost;i<mm1;i++)
   for(j=ghost;j<mm2;j++)
      for(k=ghost;k<mm3;k++) {
           temp=0;
           for(l=0;l<3;l++)
            temp+=dr(f1[l],i,j,k,l+1,0,dx[l],ghost, approx);
           if (fabs(temp)>div) div=fabs(temp);
           }
time(&time_now);
if(MainWindow->CheckScreenOutput->Checked)  {
   TempScreen.printf("program is working %ld seconds",time_now-time_begin);
   MainWindow->Screen->Lines->Append(TempScreen);
   TempScreen.printf("t=%e \tdtdid=%e \tNIter=%d", t_cur, dtdid, count);
   MainWindow->Screen->Lines->Append(TempScreen);
   TempScreen.printf("divergention=%e", div);
   }
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
          if(MainWindow->CheckScreenOutput->Checked)  {
                TempScreen.printf("%d\t%e\t%e\t%e",l, mf, mda, mdr);
                MainWindow->Screen->Lines->Append(TempScreen);
                }
          }

         if(MainWindow->CheckScreenOutput->Checked)  {
            TempScreen.printf("Energy of pulsations=%g, total energy=%g",en,TotalEnergy);
            MainWindow->Screen->Lines->Append(TempScreen);
            }
         if(MainWindow->CheckFileOutput->Checked)  {
            fen = fileopen(NameEnergyFile,count);
            fprintf(fen,"%0.10f\n",en);
            fclose(fen);
            }
         if(MainWindow->CheckScreenOutput->Checked)  {
            TempScreen.printf("number of runge-kutt calculations=%d",enter);
            MainWindow->Screen->Lines->Append(TempScreen);
            }

if(MainWindow->CheckFileOutput->Checked)  {  //if fileoutput
/*avervx = alloc_mem_1f(m3);
if(avervx == NULL)  nrerror("\nAlloc_mem: insufficient memory!\n\a",t_cur);*/

avernu = alloc_mem_1f(m3);

for(k=ghost;k<mm3;k++)
	{
	/*  avervx[k] =*/ avernu[k] = 0;
          for(i=ghost;i<mm1;i++)
             for(j=ghost;j<mm2;j++)
		 {
	  //	   avervx[k] += f1[0][i][j][k];
                   avernu[k] += nut[i][j][k];
		  }
//	   avervx[k] /= n1*n2;
	   avernu[k] /= n1*n2;
	 }
  //putting velocities to file
        fv = fileopen(NameVFile,count);
        fprintf(fv,"{%-7.5lf}",t_cur);
        print_array1d(fv,f1[0][m1/2][0],ghost,n3);
        fclose(fv);
  //putting viscosities to file
      fnu = fileopen(NameNuFile,count);
        print_array1d(fnu,avernu,ghost,n3);
        fclose(fnu);
  //for(k=0;k<m3;k++) printf("%e\n",f1[0][5][5][k]);
  //putting kaskad variables(log energy combined from them) to file
if(MainWindow->CheckNut->Checked)
  {
  fkv=fileopen(NameKaskadVarFile,count);
  fprintf(fkv,"{");
  for(i=0;i<Ns;i++)
    {
    fprintf(fkv,"{");
    for(k=0;k<n3;k++)
        {
        mf=fabs(sha[k][i]*sha[k][i]+shb[k][i]*shb[k][i]);
        if(mf>1e-300) fprintf(fkv,"%0.10lf",log(mf));
                 else fprintf(fkv,"NAN");
        fprintf(fkv,k<n3-1 ? "," : "}");
        }
    fprintf(fkv,i<Ns-1 ? "," : "}\n");
    }
  fclose(fkv);
  }//if nu_kaskad is used
 free_mem_1f(avernu,m3);
  }//if fileoutput
}

void dump(int n1,int n2,int n3,double Re,double ****f1,double ***nut,double t_cur,long count)
{
FILE *fd;                                
AnsiString tmp;
AnsiString tmp2 = NameDumpFile;
tmp.printf("%d",num_dump);
fd = fileopen(tmp2.Insert(tmp,tmp2.Pos(".")),0);
nmessage("dump is done",t_cur);
MainWindow->ChangeStatus("","—брос данных в файл є"+tmp);
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
num_dump++;
}