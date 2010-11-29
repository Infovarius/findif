//----------------------- Utilites ---------------------//
#define LEVEL extern
#include "head.h"

void nrerror(char error_text[])
{
	fprintf(stderr,"Run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

double ****alloc_mem_4f(int mvar, int n1, int n2, int n3)
{
int i, j, k;
double  ****aaa;

   aaa = (double ****)calloc(mvar, sizeof(double ***));
   if(aaa == NULL)  nrerror("\nAlloc_mem: insufficient memory!\n\a");

   for(i = 0; i < mvar; i++) {
      aaa[i] = (double ***)calloc(n1, sizeof(double **));
      if(aaa[i] == NULL) nrerror("\nAlloc_mem: insufficient memory!\n\a");
   }

   for(i = 0; i < mvar; i++)
   for(j = 0; j < n1; j++) {
      aaa[i][j] = (double **)calloc(n2, sizeof(double *));
      if(aaa[i][j] == NULL) nrerror("\nAlloc_mem: insufficient memory!\n\a");
   }

	for(i = 0; i < mvar; i++)
   for(j = 0; j < n1; j++)
   for(k = 0; k < n2; k++) {
		aaa[i][j][k] = (double *)calloc(n3, sizeof(double));
		if(aaa[i][j][k] == NULL) nrerror("\nAlloc_mem: insufficient memory!\n\a");
      }

return(aaa);
}

void    free_mem_4f(double ****aaa, int mvar, int n1, int n2, int n3)
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
   if(aaa == NULL)  nrerror("\nAlloc_mem: insufficient memory!\n\a");

   for(i = 0; i < mvar; i++) {
      aaa[i] = (double **)calloc(n1, sizeof(double *));
      if(aaa[i] == NULL) nrerror("\nAlloc_mem: insufficient memory!\n\a");
   }

   for(i = 0; i < mvar; i++)
   for(j = 0; j < n1; j++) {
      aaa[i][j] = (double *)calloc(n2, sizeof(double));
      if(aaa[i][j] == NULL) nrerror("\nAlloc_mem: insufficient memory!\n\a");
   }

return(aaa);
}

void    free_mem_3f(double ***aaa, int n1, int n2, int n3)
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

