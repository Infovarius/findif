#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include <conio.h>

#define for2D(l)  ,l==1?l++:l

LEVEL double ****f, ****f1, ****df, ****df2, ****df3, ****df4, ****df5,
             dx[3];

LEVEL double ***nut,
             Re, Gamma, l1, l2, l3, p1, p2;

LEVEL int    m1, m2, m3, mm1, mm2, mm3, n1, n2, n3,
             approx, ghost, nvar;

LEVEL long count;
LEVEL long enter;

/*------- util.c -----------*/
void   nrerror(char []);
double ****alloc_mem_4f(int, int, int, int);
void   free_mem_4f(double ****, int, int, int, int);
double ***alloc_mem_3f(int, int, int);
void   free_mem_3f(double ***, int, int, int);
/*------- rkad.c -----------*/
void   timestep(double ****, double ****, double, double ****,
              double, double *, double *);
double rkck(double ****, double ****, double, double, double ****);
/*------- pde.c -----------*/
void   pde(double, double ****, double ****);
void   boundary_conditions(double ****);
void  init_conditions(double ****);
double dr(double ***, int, int, int, int, int, double, int, int);
double coordin(int, int);
/*------- printing.c ------*/
FILE *fileopen(const char *x, int mode);
void print_array1d(FILE *ff,double *a,int beg1,int n1);
void print_array2d(FILE *ff,double **a,int beg1,int n1,int beg2,int n2);
void print_array3d(FILE *ff,double ***a,
        int beg1,int n1,int beg2,int n2,int beg3,int n3);
void printing(double ****f1,double dtdid,double t_cur,long count);
void dump(double ****f1,double t_cur,long count);
void print_array3d2(FILE *ff,double ***a,
        int beg1,int n1,int beg3,int n3);
