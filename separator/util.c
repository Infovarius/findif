//----------------------- Utilites ---------------------//
#define LEVEL extern
#include "head.h"

void nrerror(char error_text[],double t_cur,long count)
{
    FILE *err;
    nmessage(error_text,t_cur,count);
    err=fileopen(NameErrorFile,0);

    fprintf(err,"Run-time error of proc#%d at t=%-6.4lf:\n",rank,t_cur);
    fprintf(err,"%s\n",error_text);
    fprintf(err,"...now exiting to system...\n");
    fileclose(err);
    if(f) operate_memory(-1);
    add_control_point("END");
    MPI_Finalize();
    exit(1);
}

void nmessage(char msg_text[],double t_cur,long count)
{
   FILE *msg;
   msg=fileopen(NameMessageFile,1);
   time_now = (t_cur<0)?time_begin:MPI_Wtime();
   fprintf(msg,"message of proc#%d at t=%-7.4lf Niter=%-6d time of work=%g sec:\n",rank,
                    t_cur,count,time_now-time_begin);
   fprintf(msg,"%s\n",msg_text);
   fileclose(msg);
}

void add_control_point(char *name_cp)
{
FILE *cpf;
cpf = fileopen(NameCPFile,1);
if(cpf==NULL) cpf=fileopen(NameCPFile,0);
fprintf(cpf,"%s\n",name_cp);
fileclose(cpf);
}

double ****alloc_mem_4f(int mvar, int n1, int n2, int n3)
{
int i, j, k;
double  ****aaa;

   aaa = (double ****)calloc(mvar, sizeof(double ***));
   if(aaa == NULL)  nrerror("\nAlloc_mem: insufficient memory!\n",-1,-1);

   for(i = 0; i < mvar; i++) {
         aaa[i] = (double ***)calloc(n1, sizeof(double **));
         if(aaa[i] == NULL) nrerror("\nAlloc_mem: insufficient memory!\n",-1,-1);
   }

   for(i = 0; i < mvar; i++)
   for(j = 0; j < n1; j++) {
         aaa[i][j] = (double **)calloc(n2, sizeof(double *));
         if(aaa[i][j] == NULL) nrerror("\nAlloc_mem: insufficient memory!\n",-1,-1);
   }

   for(i = 0; i < mvar; i++)
   for(j = 0; j < n1; j++)
   for(k = 0; k < n2; k++) {
         aaa[i][j][k] = (double *)calloc(n3, sizeof(double));
         if(aaa[i][j][k] == NULL) nrerror("\nAlloc_mem: insufficient memory!\n",-1,-1);
      }

return(aaa);
}

void free_mem_4f(double ****aaa, int mvar, int n1, int n2, int n3)
{
int i, j, k;
   for(i = 0; i < mvar; i++)
   for(j = 0; j < n1; j++)
   for(k = 0; k < n2; k++)
	free(aaa[i][j][k]);

   for(i = 0; i < mvar; i++)
   for(j = 0; j < n1; j++)
	free(aaa[i][j]);

   for(i = 0; i < mvar; i++)
	free(aaa[i]);

   free(aaa);

   return;
}

double ***alloc_mem_3f(int mvar, int n1, int n2)
{
int i, j;
double  ***aaa;

   aaa = (double ***)calloc(mvar, sizeof(double **));
   if(aaa == NULL)  nrerror("\nAlloc_mem: insufficient memory!\n",-1,-1);

   for(i = 0; i < mvar; i++) {
         aaa[i] = (double **)calloc(n1, sizeof(double *));
         if(aaa[i] == NULL) nrerror("\nAlloc_mem: insufficient memory!\n",-1,-1);
   }

   for(i = 0; i < mvar; i++)
   for(j = 0; j < n1; j++) {
         aaa[i][j] = (double *)calloc(n2, sizeof(double));
         if(aaa[i][j] == NULL) nrerror("\nAlloc_mem: insufficient memory!\n",-1,-1);
   }

return(aaa);
}

void free_mem_3f(double ***aaa, int n1, int n2, int n3)
{
int j, k;

   for(j = 0; j < n1; j++)
   for(k = 0; k < n2; k++)
	free(aaa[j][k]);

   for(j = 0; j < n1; j++)
	free(aaa[j]);

   free(aaa);

   return;
}

double **alloc_mem_2f(int mvar, int n1)
{
int i;
double  **aa;

   aa = (double **)calloc(mvar, sizeof(double *));
   if(aa == NULL)  nrerror("\nAlloc_mem: insufficient memory!\n",-1,-1);

   for(i = 0; i < mvar; i++) {
         aa[i] = (double *)calloc(n1, sizeof(double));
         if(aa[i] == NULL) nrerror("\nAlloc_mem: insufficient memory!\n",-1,-1);
   }

return(aa);
}

void free_mem_2f(double **aa, int n1, int n2)
{
int j, k;

   for(j = 0; j < n1; j++)
	free(aa[j]);

   free(aa);

   return;
}

double *alloc_mem_1f(int n)
{
double  *a;

   a = (double *)calloc(n, sizeof(double));
   if(a == NULL)  nrerror("\nAlloc_mem: insufficient memory!\n",-1,-1);

return(a);
}

void free_mem_1f(double *a, int n)
{
int k;
   free(a);
   return;
}

void operate_memory(int dir)
{
 if(dir>0)
   {   s_func = alloc_mem_2f(n3+2,kol_masht);
       f  =alloc_mem_4f(nvar, m1, m2, m3);   //f[3]-pressure,f[0..2]-v(vector),f[4..6]-potential of induction(vector)
       f1 =alloc_mem_4f(nvar, m1, m2, m3);
       df =alloc_mem_4f(nvar, m1, m2, m3);
       df2=alloc_mem_4f(nvar, m1, m2, m3);
       df3=alloc_mem_4f(nvar, m1, m2, m3);
       df4=alloc_mem_4f(nvar, m1, m2, m3);
       df5=alloc_mem_4f(nvar, m1, m2, m3);
       B = alloc_mem_4f(3, m1, m2, m3);      // vector of magnetic induction itself
       J = alloc_mem_4f(3, m1, m2, m3);      // vector of electric current
       r_1=alloc_mem_1f(m1);                  // r^(-1)
       r_2=alloc_mem_1f(m1);                  // r^(-2)
       node=alloc_mem_2f(m1, m3);         // kind of nodes
       refr_f=alloc_mem_2f(m1, m3);         // reflection of nodes relative inner surface (for hydro)
       refz_f=alloc_mem_2f(m1, m3);
       refr_m=alloc_mem_2f(m1, m3);         // reflection of nodes relative sphere of simulation (for magn)
       refz_m=alloc_mem_2f(m1, m3);
       sinth=alloc_mem_2f(m1, m3);
       costh=alloc_mem_2f(m1, m3);
       chi=alloc_mem_2f(m1, m3);          // angle of blade tilt
       nut = alloc_mem_3f(m1, m2, m3);
       eta = alloc_mem_3f(m1, m2, m3);
       averf = alloc_mem_3f(nvar,m1,m3);
       vfi = alloc_mem_1f(N3);
       totvfi = alloc_mem_1f(N3);
       init_shell();
    } else
   {
//   free_mem_2f(s_func,n3+2,kol_masht);
       free_mem_4f(f  ,nvar, m1, m2, m3);
       free_mem_4f(f1 ,nvar, m1, m2, m3);
       free_mem_4f(df ,nvar, m1, m2, m3);
       free_mem_4f(df2,nvar, m1, m2, m3);
       free_mem_4f(df3,nvar, m1, m2, m3);
       free_mem_4f(df4,nvar, m1, m2, m3);
       free_mem_4f(df5,nvar, m1, m2, m3);
       free_mem_4f(B, 3, m1, m2, m3);
       free_mem_4f(J, 3, m1, m2, m3);
       free_mem_3f(nut, m1, m2, m3);
       free_mem_2f(refr_f,m1, m3);
       free_mem_2f(refz_f,m1, m3);
       free_mem_2f(refr_m,m1, m3);
       free_mem_2f(refz_m,m1, m3);
       free_mem_2f(sinth,m1, m3);
       free_mem_2f(costh,m1, m3);
       free_mem_2f(chi,  m1, m3);
       free_mem_2f(node,m1, m3);
       free_mem_3f(averf,nvar,m1,m3);
       free_mem_1f(vfi,N3);
       free_mem_1f(totvfi,N3);
       free_mem_1f(r_1,m1);
       free_mem_1f(r_2,m1);
       free_mem_3f(eta,m1,m2,m3);
       erase_shell();
    }
}

int isType(double nod, enum TypeNodes tip)
{
  int tmpnod;
  tmpnod = floor(nod+0.5);
  if(tip==NodeUnknown) return (!tmpnod);
                  else return (tmpnod & tip);
}

void setType(double *nod, enum TypeNodes tip)
{
  int tmpnod;
  tmpnod = floor(*nod+0.5);
  *nod = tmpnod | tip;
  return;
}

 void CopyBufferToGrid(double ****m,double ***nut,double *buffer,int x1,int y1,int z1,int x2,int y2,int z2)
 {
   int i,j,k,l,n=0;
   for(l=0;l<nvar;l++)
    for(i=x1;i<=x2;i++)
     for(j=y1;j<=y2;j++)
      for(k=z1;k<=z2;k++)
      m[l][i][j][k]=buffer[n++];
   for(i=x1;i<=x2;i++)
    for(j=y1;j<=y2;j++)
     for(k=z1;k<=z2;k++)
     nut[i][j][k]=buffer[n++];
 }

 void CopyGridToBuffer(double ****m,double ***nut,double *buffer,int x1,int y1,int z1,int x2,int y2,int z2)
 {
   int i,j,k,l,n=0;
   for(l=0;l<nvar;l++)
    for(i=x1;i<=x2;i++)
     for(j=y1;j<=y2;j++)
      for(k=z1;k<=z2;k++)
      buffer[n++]=m[l][i][j][k];
   for(i=x1;i<=x2;i++)
    for(j=y1;j<=y2;j++)
     for(k=z1;k<=z2;k++)
     buffer[n++]=nut[i][j][k];
 }

