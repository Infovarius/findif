#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <conio.h>


LEVEL double ****f, ****f1, ****df, ****df2, ****df3, ****df4, ****df5,**s_func,
             dx[3];

LEVEL double ***nut,
             Gamma, l1, l2, l3, p1, p2,Re;

LEVEL int    m1, m2, m3, mm1, mm2, mm3, n1, n2, n3, Ns,
             approx, ghost, nvar;

LEVEL double t_cur,Ttot;

LEVEL long count;
LEVEL long time_begin,time_now;
LEVEL long enter;
LEVEL bool razlet;
const int kol_masht=10;
#define min_d min(min(dx[0],dx[1]),dx[2])

/*-------- goy.cpp -----------*/
LEVEL double *nl,**sha,**shb;
void init_shell(void);
void sima(double *zza, double *zzb, double *xx, double *yy);
void time_step_shell(double tm);
/*------- struct_func.cpp----*/
void struct_func(double ****f,int q,double lambda,int size_okr);
/*------- util.cpp -----------*/
void nrerror(char[],double);
void nmessage(char [],double);
double ****alloc_mem_4f(int, int, int, int);
void   free_mem_4f(double ****, int, int, int, int);
double ***alloc_mem_3f(int, int, int);
void   free_mem_3f(double ***, int, int, int);
double **alloc_mem_2f(int, int);
void   free_mem_2f(double **, int, int);
double *alloc_mem_1f(int);
void   free_mem_1f(double *, int);
/*------- rkad.cpp -----------*/
void   timestep(double ****, double ****, double, double ****,
              double, double *, double *);
double rkck(double ****, double ****, double, double, double ****);
/*------- pde.cpp -----------*/
void   pde(double, double ****, double ****);
double norma(double a,double b,double c,int q);
double deviation(double ****,int,int,int);
void nut_by_flux(double ****f,double Re);
void   boundary_conditions(double ****);
void  init_conditions(double ****,double Re);
double dr(double ***, int, int, int, int, int, double, int, int);
double coordin(int, int);
/*------- printing.cpp ------*/
FILE *fileopen(const char *, int mode);
void print_array1d(FILE *ff,double *a,int beg1,int n1);
void print_array2d(FILE *ff,double **a,int beg1,int n1,int beg2,int n2);
void print_array3d(FILE *ff,double ***a,
        int beg1,int n1,int beg2,int n2,int beg3,int n3);
double check(double ****f);
void printing(double ****f1,double dtdid,double t_cur,long count,double en);
void dump(double ****f1,double t_cur,long count);


