//----------------------- Utilites ---------------------//
#define LEVEL extern
#include "head.h"

void nrerror(char error_text[],double t_cur)
{
   FILE *err;
   nmessage(error_text,t_cur);
   err=fileopen("error.err",0);

	fprintf(err,"Run-time error at t=%-6.4lf:\n",t_cur);
	fprintf(err,"%s\n",error_text);
	fprintf(err,"...now exiting to system...\n");
        fclose(err);
	exit(1);
}

void nmessage(char msg_text[],double t_cur)
{
   FILE *msg;
   msg=fileopen("message.dat",1);
   time(&time_now);
   fprintf(msg,"message at t=%-7.4lf Niter=%-6d time of work=%ld:\n",
                    t_cur,count,time_now-time_begin);
   fprintf(msg,"%s\n",msg_text);
   fclose(msg);
}

double ****alloc_mem_4f(int mvar, int n1, int n2, int n3)
{
int i, j, k;
double  ****aaa;

   aaa = (double ****)calloc(mvar, sizeof(double ***));
   if(aaa == NULL)  nrerror("\nAlloc_mem: insufficient memory!\n\a",0);

   for(i = 0; i < mvar; i++) {
      aaa[i] = (double ***)calloc(n1, sizeof(double **));
      if(aaa[i] == NULL) nrerror("\nAlloc_mem: insufficient memory!\n\a",0);
   }

   for(i = 0; i < mvar; i++)
   for(j = 0; j < n1; j++) {
      aaa[i][j] = (double **)calloc(n2, sizeof(double *));
      if(aaa[i][j] == NULL) nrerror("\nAlloc_mem: insufficient memory!\n\a",0);
   }

	for(i = 0; i < mvar; i++)
   for(j = 0; j < n1; j++)
   for(k = 0; k < n2; k++) {
		aaa[i][j][k] = (double *)calloc(n3, sizeof(double));
		if(aaa[i][j][k] == NULL) nrerror("\nAlloc_mem: insufficient memory!\n\a",0);
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
   if(aaa == NULL)  nrerror("\nAlloc_mem: insufficient memory!\n\a",0);

   for(i = 0; i < mvar; i++) {
      aaa[i] = (double **)calloc(n1, sizeof(double *));
      if(aaa[i] == NULL) nrerror("\nAlloc_mem: insufficient memory!\n\a",0);
   }

   for(i = 0; i < mvar; i++)
   for(j = 0; j < n1; j++) {
      aaa[i][j] = (double *)calloc(n2, sizeof(double));
      if(aaa[i][j] == NULL) nrerror("\nAlloc_mem: insufficient memory!\n\a",0);
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
   if(aa == NULL)  nrerror("\nAlloc_mem: insufficient memory!\n\a",0);

   for(i = 0; i < mvar; i++) {
      aa[i] = (double *)calloc(n1, sizeof(double));
      if(aa[i] == NULL) nrerror("\nAlloc_mem: insufficient memory!\n\a",0);
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
   if(a == NULL)  nrerror("\nAlloc_mem: insufficient memory!\n\a",0);

return(a);
}

void free_mem_1f(double *a, int n)
{
int k;
   free(a);
   return;
}


