#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include <conio.h>
#include "mpi.h"

LEVEL char *NameMessageFile;
LEVEL char *NameErrorFile;
LEVEL char *NameInitFile;
LEVEL char *NameNuFile;
LEVEL char *NameVFile;
LEVEL char NameDumpFile[80];
LEVEL char *NameSnapFile;
LEVEL char *NameEnergyFile;
LEVEL char *KaskadVarFile;
LEVEL char  fname[80], fname_stat[80];

LEVEL double ****f, ****f1, ***nut, ****df, ****df2, ****df3, ****df4, ****df5, **s_func,
             dx[3], *r_1, *r_2, **node, **refr, **refz, **sinth, **costh, **chi;
                              //1-outside;2-ghost node;3-integrating area
LEVEL double ***averf;
LEVEL double *totvfi, *vfi;
LEVEL double Gamma, lfi, R, rc, Re;

LEVEL int    m1, m2, m3, mm1, mm2, mm3, n1, n2, n3, N1, N2, N3, Ns,
             approx, ghost, nvar;

LEVEL double Noise, NoiseNorm, parabole;
LEVEL double chimax;           //maximal angle of blade tilt
LEVEL double maschtab;

LEVEL double t_cur,Ttot;
LEVEL long count;
LEVEL double time_begin,time_now;

LEVEL long enter;
LEVEL char razlet;
LEVEL double TotalEnergy, PulsEnergy;
LEVEL int CheckStep;  //each CheckStep divergence of program is checking
LEVEL int OutStep;    //each OutStep result is printing
LEVEL int SnapStep;   //each SnapStep snapshot is done
LEVEL int VarStep;    //each VarStep step in RungeKutta varies
LEVEL int num_dump;
LEVEL long numlog;

#define max(a, b)  (((a) > (b)) ? (a) : (b))
#define min(a, b)  (((a) < (b)) ? (a) : (b))
#define min_d min(min(dx[0],dx[1]),dx[2])
/*----------------parallel---------------*/
#define Master if(rank==0)
LEVEL int rank, size;
LEVEL MPI_Request SendRequest[256],RecvRequest[256];
LEVEL MPI_Status statuses[256];
LEVEL double *buf_send[6], *buf_recv[6], buf_size[3];
LEVEL int pr_neighbour[6],pp[3],n[3];
LEVEL int pr[3];
//void init_parallel();                       in file pde.c

//--------for GOY
#define kol_masht 10
LEVEL int max_okr;  //должно быть <=ghost, иначе окрестность выходит из фикт.обл-ти - виснет
LEVEL double lambda; //multiplier
LEVEL long num_points[kol_masht];
LEVEL double *nl,**sha,**shb;
LEVEL double *f0;                           //flux of energy from grid to shell
LEVEL double *en_flux;                      //summary of energy in cascade model

/*-------- goy.c  -----------*/
void init_shell(void);
void sima(double *zza, double *zzb, double *xx, double *yy);
void time_step_shell(double tm);
/*------- struct_func.c ----*/
void struct_func(double ****f,int q,double lambda,int size_okr);
/*------- util.c  -----------*/
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
/*------- rkad.c  -----------*/
void timestep(double ****f, double ****df, double t, double ****fout,
              double dttry, double *dtdid, double *dtnext);
double rkck(double ****, double ****, double, double, double ****);
/*------- pde.c  -----------*/
void   pde(double, double ****, double ****);
double norma(double a,double b,double c,int q);
double deviation(double ****,int,int,int);
void nut_by_flux(double ****f,double Re);
void   boundary_conditions(double ****);
void  init_conditions(double ****,double Re);
void init_parallel(void);
double dr(double ***, int, int, int, int, int, double, int, int);
double coordin(int, int);
double interpolation(double,double, double,double,double,double);
void CopyBufferToGrid(double ****,double ***,double *,int,int,int,int,int,int);
void CopyGridToBuffer(double ****,double ***,double *,int,int,int,int,int,int);
/*------- ioutil.c  ------*/
FILE *fileopen(const char *, int mode);
void init_param(int argc, char** argv,double *dtnext);
int init_data(double *Re,double *t_cur,long *count);
void print_array1d(FILE *ff,double *a,int beg1,int n1);
void print_array2d(FILE *ff,double **a,int beg1,int n1,int beg2,int n2);
void print_array3d(FILE *ff,double ***a,int beg1,int n1,int beg2,int n2,int beg3,int n3);
void check(double ****f);
void printing(double ****f1,double dtdid,double t_cur,long count,double en);
void dump(double ****f1,double t_cur,long count);
void snapshot(double ****f,double t_cur,long count);
void putlog(char msg_text[],long num);
