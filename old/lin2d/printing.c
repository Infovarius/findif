//------------------------ all the outputting stuff -----------------//

#define LEVEL extern
#include <conio.h>
#include "head.h"

LEVEL const char *NameNuFile = "nut.dat";
LEVEL const char *NameVFile  = "vv.dat";
LEVEL const char *NameDumpFile = "dump.dat";
LEVEL const char *NameEnergyFile = "energy.dat";

const double UpLimit=3.;     //after this limit there's dump
FILE *fileopen(const char *x, int mode)  //opening of file to ff
                         /*   0-rewrite;others-append   */
{
FILE *ff;
if(mode) ff=fopen(x,"a");
else if ((ff = fopen (x,"w+"))==NULL)
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

void print_array3d2(FILE *ff,double ***a,
        int beg1,int n1,int beg3,int n3)
{
int i,k;
fprintf(ff,"{");
for(i=beg1;i<beg1+n1;i++) {
    fprintf(ff,"{");
        for(k=beg3;k<beg3+n3;k++) {
            fprintf(ff,"%0.10f",a[i][0][k]);
            fprintf(ff,k<beg3+n3-1 ? "," : "}");
            }
    fprintf(ff,i<beg1+n1-1 ? "," : "}\n");
    }
}

void printing(double ****f1,double dtdid,double t_cur,long count)
{
static int net;
double temp, div=0;
int i,k,l;
double mf, mda, mdr, en;
double *avervx, *avernu;
FILE *fv, *fnu, *fen;
//fnu = fileopen(NameNuFile,count);
clrscr();
boundary_conditions(f1);

(count==0? net = 1 : net);

for(i=ghost;i<mm1;i++)
      for(k=ghost;k<mm3;k++) {
           temp=0;
           for(l=0;l<3;l++ for2D(l))
            temp+=dr(f1[l],i,0,k,l+1,0,dx[l],ghost, approx);
//           fprintf(fnu,"%g%15g%15g%15g%15g\n",coordin(i,0),coordin(k,2),temp,
//                dr(f1[0],i,0,k,1,0,dx[0],ghost, approx), dr(f1[2],i,0,k,3,0,dx[2],ghost, approx));
           if (fabs(temp)>div) div=fabs(temp);
           }
//           fclose(fnu);
printf("%e %e %d %e\n", t_cur, dtdid, count, div);

        for(l=0;l<nvar;l++ for2D(l)) {
          mf=mda=mdr=en=0;
          for(i=ghost;i<mm1;i++)
          for(k=ghost;k<mm3;k++) {
            if (l==0) en+=pow(f1[l][i][0][k],2);
               else en+=pow(f1[l][i][0][k],2);
            temp=fabs(f[l][i][0][k]-f1[l][i][0][k]);
            if (temp>mda) mda=temp;
            if(f1[l][i][0][k]!=0) temp/=f1[l][i][0][k];
            if (temp>mdr) mdr=temp;
            if (fabs(f1[l][i][0][k])>mf) mf=fabs(f1[l][i][0][k]);
          }
          if(l==0&&mf>UpLimit&&net) {dump(f1,t_cur,count); net = 0;};
          printf("%d %e %e %e\n",l, mf, mda, mdr);
          }
printf("Energy of pulsations=%g\n",en);
fen = fileopen(NameEnergyFile,count);
fprintf(fen,"%0.10f\n",en);
fclose(fen);
printf("number of runge-kutt calculations=%d",enter);

avervx = (double *)calloc(m3, sizeof(double));
if(avervx == NULL)  nrerror("\nAlloc_mem: insufficient memory!\n\a");

avernu = (double *)calloc(m3, sizeof(double));
if(avernu == NULL)  nrerror("\nAlloc_mem: insufficient memory!\n\a");

for(k=ghost;k<mm3;k++ for2D(k))
	{
	  avervx[k] = avernu[k] = 0;
          for(i=ghost;i<mm1;i++)
		 {
		   avervx[k] += f1[0][i][0][k];
                   avernu[k] += nut[i][0][k];
		  }
	   avervx[k] /= n1;
	   avernu[k] /= n1;
	 }

//putting velocities to file
        fv = fileopen(NameVFile,count);
        print_array1d(fv,f1[0][m1/2][0],ghost,n3);
        fclose(fv);
//putting viscosities to file
/*        fnu = fileopen(NameNuFile,count);
        print_array1d(fnu,avernu,ghost,n3);
        fclose(fnu);*/

//for(k=0;k<m3;k++) printf("%e\n",f1[0][5][0][k]);

}

void dump(double ****f1,double t_cur,long count)
{
FILE *fd;
fd = fileopen(NameDumpFile,0);
fprintf(fd,"{%0.10f ,  %ld}\n",t_cur,count);
print_array3d2(fd,f1[0],0,m1,0,m3);
print_array3d2(fd,f1[2],0,m1,0,m3);
print_array3d2(fd,f1[3],0,m1,0,m3);
fclose(fd);
}
