#ifndef headH
#define headH
//---------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <conio.h>

#define for2D(l)  ,l==1?l++:l

LEVEL enum {ActNone=0, ActBegin=1, ActDump=2, ActPause=4, ActContinue=8, ActQuit=16}
     action;
LEVEL enum {StNone=0, StPause=1, StRun=2} RunStatus;

LEVEL char *NameMessageFile;
LEVEL char *NameErrorFile;
LEVEL char *NameInitFile;
LEVEL char *NameNuFile;
LEVEL char *NameVFile;
LEVEL char *NameDumpFile;
LEVEL char *NameEnergyFile;
LEVEL char *KaskadVarFile;

LEVEL double ****f, ****f1, ****df, ****df2, ****df3, ****df4, ****df5,**s_func,
             dx[3];

LEVEL double ***nut,
             Gamma, l1, l2, l3, p1, p2,Re;

LEVEL int    m1, m2, m3, mm1, mm2, mm3, n1, n2, n3, Ns,
             approx, ghost, nvar;

LEVEL double Noise;
LEVEL double NoiseNorm;
LEVEL double maschtab;

LEVEL double t_cur,Ttot;
LEVEL long count;
LEVEL long time_begin,time_now;

LEVEL long enter;
LEVEL bool razlet;
LEVEL double TotalEnergy, PulsEn;
LEVEL int CheckStep;  //each CheckStep divergence of program is checking
LEVEL int OutStep;    //each OutStep result is printing
LEVEL int num_dump;

//--------for GOY
#define kol_masht 5
LEVEL int max_okr;  //должно быть <=ghost, иначе окрестность выходит из фикт.обл-ти - виснет
LEVEL double lambda; //multiplier
LEVEL long num_points[kol_masht];
LEVEL double *nl,**sha,**shb;
LEVEL double *f0;                           //flux of energy from grid to shell
LEVEL double *en_flux;                      //summary of energy in cascade model

/*------- printing.cpp ---------*/
int init_data(double *Re,double *t_cur,long *count);
FILE *fileopen(char *, int mode);
void print_array1d(FILE *ff,double *a,int beg1,int n1);
void print_array2d(FILE *ff,double **a,int beg1,int n1,int beg2,int n2);
void print_array3d(FILE *ff,double ***a,
        int beg1,int n1,int beg2,int n2,int beg3,int n3);
double check(double ****f);
void printing(double ****f1,double dtdid,double t_cur,long count,double en);
void dump(int n1,int n2,int n3,double Re,double ****f,double ***nut,double t_cur,long count);
/*-------- goy.cpp -----------*/
void init_shell(void);
void sima(double *zza, double *zzb, double *xx, double *yy,int layer);
void time_step_shell(double tm,double *en_flux);
void erase_shell(void);
/*------- struct_func.cpp----*/
void struct_func(double ****f,int q);
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
/*------- kaskad.cpp --------*/
void allocate_memory();
void free_memory(double **s_func,double ****f,double ****f1,double ****df,double ****df2,
                     double ****df3,double ****df4,double ****df5,double ***nut);
void at_the_end_of_work();
void recreate_arrays(int n1new,int n2new,int n3new);
#endif
