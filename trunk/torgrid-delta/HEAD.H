#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include <conio.h>
#include "mpi.h"

LEVEL char *NameMessageFile;
LEVEL char *NameErrorFile;
LEVEL char NameInitFile[200];
LEVEL char *MainFlowFile;
LEVEL char *NameVFile;
LEVEL char *NameDumpFile;
LEVEL char *NameSnapFile;
LEVEL char *NameEnergyFile;
LEVEL char *NameCPFile;
LEVEL char *NameStatFile;
LEVEL char *fname;

LEVEL double ****flow, ****dflow, ****f, ****f1, ***nut, ****df, ****df2, ****df3, ****df4, ****df5,
             dx[3], *r_1, *r_2, **refr, **refz, **sinth, **costh, **chi;
LEVEL int **node;

enum  TypeNodes {NodeUnknown=0,NodeClued=1,
                               NodeFluid=8, NodeShell=2, NodeVacuum=4,   // regions
                               /*NodeFluid=8,*/ NodeGhostFluid=16};       // for hydrodynamics
LEVEL double ***averf;
LEVEL double *totvfi, *vfi;
LEVEL double Gamma, lfi, R, rc, Re, p1;

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

LEVEL int goon;                       // continue calculation
LEVEL int CheckStep;  //each CheckStep divergence of program is checking
LEVEL int OutStep;    //each OutStep result is printing
LEVEL int SnapStep;   //each SnapStep snapshot is done
LEVEL double SnapDelta;  //each SnapDelta of simulation time snapshot is done
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
LEVEL double *buf_send[6], *buf_recv[6], buf_size[3];    //6=2*3(D)
LEVEL int pr_neighbour[6],pp[3],n[3];                    //3D
LEVEL int pr[3];                                         //3D
//void init_parallel();                       in file pde.c

//--------for GOY
#define kol_masht 10
LEVEL int max_okr;  //должно быть <=ghost, иначе окрестность выходит из фикт.обл-ти - виснет
LEVEL double lambda; //multiplier
LEVEL long num_points[kol_masht];

/*------- util.c  -----------*/
void nrerror(char *mess,double t_cur,long count);
void nmessage(char *mess,double t_cur,long count);
void add_control_point(char *);
double ****alloc_mem_4f(int, int, int, int);
void   free_mem_4f(double ****, int, int, int, int);
double ***alloc_mem_3f(int, int, int);
void   free_mem_3f(double ***, int, int, int);
double **alloc_mem_2f(int, int);
void   free_mem_2f(double **, int, int);
double *alloc_mem_1f(int);
void   free_mem_1f(double *, int);
int **alloc_mem_2i(int, int);
void   free_mem_2i(int **, int, int);
int *alloc_mem_1i(int);
void   free_mem_1i(int *, int);
void operate_memory(int dir);
int isType(int nod, enum TypeNodes tip);
void setType(int *nod, enum TypeNodes tip);
void CopyBufferToGrid(double ****,double ***,double *,int,int,int,int,int,int);
void CopyGridToBuffer(double ****,double ***,double *,int,int,int,int,int,int);
/*------- rkad.c  -----------*/
void timestep(double ****f, double ****df, double ***nu, double t, double ****fout,
              double dttry, double *dtdid, double *dtnext);
double rkck(double ****f, double ****df, double ***nut, double t, double dt, double ****fout);
/*------- pde.c  -----------*/
void   pde(double t, double ****f, double ****df);
double norma(double a,double b,double c,int q);
double deviation(double ****f,int,int,int);
void nut_by_flux(double ****f, double ***nu, double Re);
void boundary_conditions(double ****f, double ***nu);
void init_conditions(void);
void init_parallel(void);
void precalc(void);
double dr(double ***m,int ii,int jj,int kk,int dir,int or,double dx,int sh,int sm);
double coordin(int, int);
double interpolation(double,double, double,double,double,double);
/*------- ioutil.c  ------*/
FILE *fileopen(const char *, int mode);
void fileclose(FILE *fd);
void init_param(int argc, char** argv,double *dtnext);
int init_data(char *fname, double ****f);
void print_array1d(FILE *ff,double *a,int beg1,int n1);
void print_array2d(FILE *ff,double **a,int beg1,int n1,int beg2,int n2);
void print_array3d(FILE *ff,double ***a,int beg1,int n1,int beg2,int n2,int beg3,int n3);
void print_array2i(FILE *ff,int **a,int beg1,int n1,int beg2,int n2);
void check(double ****f);
void printing(double ****f1,double dtdid,double t_cur,long count,double en);
void dump(double ****f1,double ***nu,double t_cur,long count);
void snapshot(double ****f,double ***nu,double t_cur,long count);
void putlog(char msg_text[],long num);
