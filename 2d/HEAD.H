#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include <conio.h>


LEVEL double ****f, ****f1, ****df, ****df2, ****df3, ****df4, ****df5,
             dx[3];

LEVEL double ***nut,
             Re, Gamma, l1, l2, l3, p1, p2;

LEVEL int    m1, m2, m3, mm1, mm2, mm3, n1, n2, n3,
             approx, ghost, nvar;

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
void   boudary_conditions(double ****);
void  init_conditions(double ****);
double dr(double ***, int, int, int, int, int, double, int, int);
double coordin(int, int);



