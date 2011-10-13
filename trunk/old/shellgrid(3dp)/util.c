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
    { 
    if(!aaa[i][j][k]) {putlog("error at releasing",(long)aaa[i][j][k]);return;}
	free(aaa[i][j][k]);
	}

   for(i = 0; i < mvar; i++)
   for(j = 0; j < n1; j++)
    {
    if(!aaa[i][j]) {putlog("error at releasing",(long)aaa[i][j]);return;}
	free(aaa[i][j]);
	}

   for(i = 0; i < mvar; i++)
    {
    if(!aaa[i]) {putlog("error at releasing",(long)aaa[i]);return;}
	free(aaa[i]);
	}

   if(!aaa) {putlog("error at releasing",(long)aaa);return;}
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
    {
    if(!aaa[j][k]) {putlog("error at releasing",(long)aaa[j][k]);return;}
	free(aaa[j][k]);
	}

   for(j = 0; j < n1; j++)
    {
    if(!aaa[j]) {putlog("error at releasing",(long)aaa[j]);return;}
	free(aaa[j]);
	}

    if(!aaa) {putlog("error at releasing",(long)aaa);return;}
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
    {
    if(!aa[j]) {putlog("error at releasing",(long)aa[j]);return;}
	free(aa[j]);
    }
   if(!aa) {putlog("error at releasing",(long)aa);return;}
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
   if(!a) {putlog("error at releasing",(long)a);return;}
   free(a);
   return;
}

void operate_memory(int dir)
{
int n=max(n1,max(n2,n3));
 if(dir>0)
   {   s_func = alloc_mem_2f(n3+2,kol_masht);
       f  =alloc_mem_4f(nvar, m1, m2, m3);   //f[3]-pressure,f[0..2]-v(vector)
       f1 =alloc_mem_4f(nvar, m1, m2, m3);
       df =alloc_mem_4f(nvar, m1, m2, m3);
       df2=alloc_mem_4f(nvar, m1, m2, m3);
       df3=alloc_mem_4f(nvar, m1, m2, m3);
       df4=alloc_mem_4f(nvar, m1, m2, m3);
       df5=alloc_mem_4f(nvar, m1, m2, m3);
       nut=alloc_mem_3f(m1, m2, m3);
       averf = alloc_mem_2f(3,n);
       vx = alloc_mem_1f(N3);
       totvx = alloc_mem_1f(N3);
       avernu = alloc_mem_1f(N3);
    } else
   {
       free_mem_2f(s_func,n3+2,kol_masht);
       free_mem_4f(f  ,nvar, m1, m2, m3);
       free_mem_4f(f1 ,nvar, m1, m2, m3);
       free_mem_4f(df ,nvar, m1, m2, m3);
       free_mem_4f(df2,nvar, m1, m2, m3);
       free_mem_4f(df3,nvar, m1, m2, m3);
       free_mem_4f(df4,nvar, m1, m2, m3);
       free_mem_4f(df5,nvar, m1, m2, m3);
       free_mem_3f(nut, m1, m2, m3);
       free_mem_2f(averf,3,n);
       free_mem_1f(vx,N3);
       free_mem_1f(totvx,N3);
       free_mem_1f(avernu,N3);
    }   
}
