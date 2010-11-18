//------------------------ all the outputting stuff -----------------//

#define LEVEL extern
//#include <conio.h>
#include "head.h"

const double UpLimit=10.;     //after this limit there's dump

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
int n=max(n1,max(n2,n3));
TotalEnergy=0;
averf = alloc_mem_2f(3,n);
int i,j,k,l;
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
//fnu = fileopen(NameNuFile,count);
clrscr();
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

/*avervx = (double *)calloc(m3, sizeof(double));
if(avervx == NULL)  nrerror("\nAlloc_mem: unsuffitient memory!\n\a",t_cur);*/

avernu = (double *)calloc(m3, sizeof(double));
if(avernu == NULL)  nrerror("\nAlloc_mem: unsuffitient memory!\n\a",t_cur);

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
/*fkv=fileopen(KaskadVarFile,count);
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
fclose(fkv);*/
//putting structural functions
//struct_func(f,2,2,3);
fkv=fileopen(KaskadVarFile,count);
fprintf(fkv,"{");
for(i=0;i<kol_masht;i++)
    {
    fprintf(fkv,"{");
    for(k=0;k<n3;k++)
        {
        mf=fabs(s_func[k][i]);
        if(mf>1e-300) fprintf(fkv,"%0.10lf",log10(mf));
                 else fprintf(fkv,"NAN");
        fprintf(fkv,k<n3-1 ? "," : "}");
        }
    fprintf(fkv,i<kol_masht-1 ? "," : "}\n");
    }
fclose(fkv);
}

void dump(double ****f1,double t_cur,long count)
{
FILE *fd;
fd = fileopen(NameDumpFile,0);
nmessage("dump is done",t_cur);
fprintf(fd,"%0.10f \n  %ld\n",t_cur,count);
print_array3d(fd,f1[0],0,m1,0,m2,0,m3);
print_array3d(fd,f1[1],0,m1,0,m2,0,m3);
print_array3d(fd,f1[2],0,m1,0,m2,0,m3);
print_array3d(fd,f1[3],0,m1,0,m2,0,m3);
fclose(fd);
}

void snapshot(double ****f,long count)
{
char str[256];
int i, j, l;

 sprintf(str,"%s_%5d_%d.snp",NameSnapFile,count,rank);
 iop=fopen(str,"w");
 fprintf(iop,"%e\t%g\n",t_cur,MPI_Wtime()-time_begin);
print_array3d(iop,f1[0],0,m1,0,m2,0,m3);
print_array3d(iop,f1[1],0,m1,0,m2,0,m3);
print_array3d(iop,f1[2],0,m1,0,m2,0,m3);
print_array3d(iop,f1[3],0,m1,0,m2,0,m3);
 fclose(iop);

}