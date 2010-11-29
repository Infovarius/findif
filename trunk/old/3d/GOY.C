#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/*==============================================================================*/
/*                 Complex shell model.  Z(N) = A(N) + iB(N)                    */
/*                               GOY  -  epsilon                                */
/*                             - DOUBLE  PRECISION                              */
/*                                                                              */
/*==============================================================================*/
/*                                                                              */
/*  NN shells, and N0-th shell is considered as the 0-th scale                  */
/*  Forcing A: random forces act in shells N0 and N0+1 each TAUF step of time   */
/*  Forcing B: The force acts in shells number N0 and N0+1, and provides        */
/*             the permanent input flux of energy PQ                            */
/*  Forcing C: A[N0] and A[N0+1] are random, but                                */
/*             A[N0]**2 + B[N0]**2 = 1                                          */
/*             A[N0+1]**2 + B[N0+1]**2 = 2**(-2/3)                              */
/*  Forcing D: permanent force acts in shell N0 f_0 = F0*(1+i)                  */
/*==============================================================================*/


#define NN      18              /* NN - quantity of levels minus 1		*/
#define NTT     (NN+1)          /* NTT - total quantity of levels 		*/
#define N0      3               /* N0 - the scale 0 is in level	number N0       */
#define NG      20000          /* NG  - number of points of differ counts	*/
#define IG      50              /* IG - each IG-th point the counts are done    */
#define NGG     40              /* NGG  - size of massiv for time series  	*/
#define MT0     5000000          /* MT0 - start point of count                   */
#define MTMAX	(NG*IG+MT0)	/* MTMAX - total number of time's steps		*/

#define RE      1.0e+6
#define TAU	0.00005		/* TAU - step of time				*/
#define TAUF	100             /* each TAUF step the force acts (forcing A)    */
#define K0	1.0             /* K0 - viscous constant                        */
#define PQ      0.001            /* PQ - inputed flux of energy in the system   */
#define DIS     1                /* DIS=1 -ordinary dissip                      */
#define FORCE   4                /* FORCE=1 -forcing A,  FORCE=2 -forcing B     */
                                 /* FORCE=3 -forcing C                          */
#define F0      0.005
#define BEG     1                /* BEG=0 -initialization, BEG=1 -input of data */
#define CORR    1                /* CORR=1    calculation of correlations       */
#define TIMESER 0                /* TIMESER=1 output of time series             */
#define STRUCT  1                /* STRUCT=1  calculation of structure functions*/
#define JBENZI  15               /* JBENZI  the highest order of structure func.*/

#define EPS     -2.              /* EPS - the GOY-model parameter               */
#define INTNAME "goym05"
#define OUTNAME "goym05"

#define STATION 2                /* 1 - SUN,   2 - HP                           */
#define SUN     1073841824.0     /* conrand for SUN-station                     */
#define PH      16384.0          /* conrand for HP-station                      */

double	a[NTT], b[NTT], da[NTT], db[NTT], dnn[NTT],
        CK, CI, CJ;

FILE    *ff;

char    fn1[30]={"ab_"},       fn2[30]={"ab_"},
        fn3[30]={"en_time_"},  fn5[30]={"dis_time_"},
        fn6[30]={"sp_"},
        fn11[30]={"cor_uuu_"}, fn12[30]={"sin_"},
        fn13[30]={"zeta_"},    fn14[30]={"zeta_ess_"},    fn16[30]={"stru_"},
        fn17[30]={"eps_"},     fn18[30]={"epsrel_"},      fn19[30]={"eps_inf_"},
        fn20[30]={"pirel_"},   fn21[30]={"beta_"};

double	sign(double);
void	sima(double *, double *, double *, double *);
void    prepin(char *);
void    prepout(char *);

/*==============================================================================*/

void main(void)
{
  int    mt, i, j, k, nk, igg=0, iran, intrand;
  float  dum[NTT], dom[NTT], dam, fl1, fl2;

  double tem=0.0, tau1, en, om, sl, dis, dis1, teta, set, rek,
         varrand, conrand, input=0.0, output=0.0;
  double q0=(PQ*TAU), q1, fN0=0.010, aaa, bbb, ss;

  double nl[NTT],   dnl[NTT],   rk[NTT],      ro[NTT],    pirel[NTT],
         eni[NTT],  disi[NTT],  enisr[NTT],   dissr[NTT],
         tek[NTT][11],      sek[NTT][11],     tet[NTT],
         stru[JBENZI][NTT], lstru[JBENZI][NTT], zeta[JBENZI][NTT],
         eps[JBENZI][NTT],  epsrel[JBENZI][NTT];

  double x1[NTT], x2[NTT], x3[NTT], x4[NTT],
         y1[NTT], y2[NTT], y3[NTT], y4[NTT],
         d[NTT],  za[NTT], zb[NTT];

  double tim[NGG], e[NGG], disen[NGG];

  tau1 = sqrt(TAU);
  q1 = sqrt(0.5*q0);
  conrand = PH;
  if (STATION == 1) conrand = SUN;

  for (k=0; k<NTT; k++) {
    nk = k - N0;
    dnn[k] = pow(2.0,(double)nk);
    nl[k] = dnn[k]*dnn[k];
    enisr[k] = dissr[k] = 0.0;
    for (j=0; j<11; j++) tek[k][j] = sek[k][j] = 0.0;
    for (j=0; j<JBENZI; j++) {
      stru[j][k] = 0.0;
      eps[j][k] = 0.0;
    }
  }
    CK = -EPS/2.0;
    CJ = (EPS-1.0)/2.0;
    CI = 0.5;
    fl1 = (EPS-1.0)/4.0;
    fl2 = -0.5;

  for (k=0; k<NTT; k++){
    dnl[k] = K0*nl[k]/(RE);
    rk[k] = exp(-TAU*dnl[k]);
  }
  printf ("\n  EPS= %f          re= %e \n\a\a",EPS, RE);
       printf ("\n");

/*---------------------------------------     initial conditions    -----*/
/*---------------------------------------      or input of data     -----*/
  if (BEG==1){
    for (k=1; k<N0; k++){
      intrand = rand();
      varrand = ((float)intrand - conrand)/conrand;
      a[k] = 0.00001*(2.0f+varrand)*dnn[k];
      intrand = rand();
      varrand = ((float)intrand - conrand)/conrand;
      b[k] = 0.00001*(2.0f+varrand)*dnn[k];
    }
    for (k=N0; k<NTT; k++){
      intrand = rand();
      varrand = ((float)intrand - conrand)/conrand;
    printf ("\n %e",varrand);
      a[k] = 0.00001*(2.0f+varrand) / (dnn[k]*dnn[k]);
      intrand = rand();
      varrand = ((float)intrand - conrand)/conrand;
      b[k] = 0.00001*(2.0f+varrand) / (dnn[k]*dnn[k]);
    }
    printf ("\n End of initialization !\n");
  }
  else {
    prepin(fn1);
    fscanf (ff,"%e ",&dam);
    for (k=0; k<=NN; k++) fscanf (ff,"%e %e",&dum[k], &dom[k]);
    fclose(ff);
    printf ("\n End of input !\n");

    tem = dam;
    for (k=0; k<=NN; k++) {
      a[k] = (double)dum[k];
      b[k] = (double)dom[k];
    }
  }
  printf ("\n");

/*------------------------------------------------------------------------*/
  mt = 0;
  while (mt++ < MTMAX) {
    tem += TAU;

/*------------------------------------------------  Runge - Kutta  -------*/
    for (i=0; i<NTT; i++) {
      za[i] = a[i];
      zb[i] = b[i];
    }
    sima(za,zb,x1,y1);
    for (i=0; i<NTT; i++) {
      za[i] = a[i]+x1[i]/2.0f;
      zb[i] = b[i]+y1[i]/2.0f;
    }
    sima(za,zb,x2,y2);
    for (i=0; i<NTT; i++) {
      za[i] = a[i]+x2[i]/2.0f;
      zb[i] = b[i]+y2[i]/2.0f;
    }
    sima(za,zb,x3,y3);
     for (i=0; i<NTT; i++) {
      za[i] = a[i]+x3[i];
      zb[i] = b[i]+y3[i];
    }
    sima(za,zb,x4,y4);
	       		   
    for (i=0; i<NTT; i++){
      a[i] += (x1[i]+2.0f*x2[i]+2.0f*x3[i]+x4[i])/6.0f;
      b[i] += (y1[i]+2.0f*y2[i]+2.0f*y3[i]+y4[i])/6.0f;
    }
/*--------------------- input and output of energy and enstrophy  --------*/

/*----------------------- forcing ----------------------------------------*/
  if (FORCE == 1){                              /*   Forcing  A           */
    if (!(mt % TAUF)) {
      intrand = rand();
      varrand = ((float)intrand - conrand)/conrand;
      a[N0] += tau1*fN0*varrand;
      intrand = rand();
      varrand = ((float)intrand - conrand)/conrand;
      b[N0] += tau1*fN0*varrand;
      intrand = rand();
      varrand = ((float)intrand - conrand)/conrand;
      a[N0+1] += tau1*0.5*fN0*varrand;
      intrand = rand();
      varrand = ((float)intrand - conrand)/conrand;
      b[N0+1] += tau1*fN0*0.5*varrand;
    }
  }
  if (FORCE == 2){                              /*  Forcing  B            */
    aaa = a[N0];
    intrand = rand();
    varrand = ((float)intrand - conrand)/conrand;
    a[N0] += q1*varrand/a[N0];

    en = a[N0]*a[N0]-aaa*aaa;
    if (en > q0) {
      ss = sign(a[N0]);
      a[N0] = ss*sqrt(aaa*aaa+0.5*q0);
      en = a[N0]*a[N0]-aaa*aaa;
    }
      bbb = b[N0]*b[N0]+q0-en;
      ss = sign(b[N0]);
      b[N0] = ss*sqrt(bbb);

    aaa = a[N0+1];
    intrand = rand();
    varrand = ((float)intrand - conrand)/conrand;
    a[N0+1] += q1*varrand/a[N0+1];
    en = a[N0+1]*a[N0+1]-aaa*aaa;
    if (en > q0) {
      ss = sign(a[N0+1]);
      a[N0+1] = ss*sqrt(aaa*aaa+0.5*q0);
      en = a[N0+1]*a[N0+1]-aaa*aaa;
    }
      bbb = b[N0+1]*b[N0+1]+q0-en;
      ss = sign(b[N0+1]);
      b[N0+1] = ss*sqrt(bbb);
   }
  if (FORCE == 3){                              /*  Forcing  C            */
    intrand = rand();
    varrand = ((float)intrand - conrand)/conrand;
    ss = sign(a[N0]);
    a[N0] = varrand;
    b[N0] = ss*sqrt(1.0-a[N0]*a[N0]);

    intrand = rand();
    varrand = ((float)intrand - conrand)/conrand;
    ss = sign(b[N0+1]);
    b[N0+1] = varrand*0.79;
    a[N0+1] = ss*sqrt(0.63-b[N0+1]*b[N0+1]);
  }
  if (FORCE == 4){                              /*   Forcing  D          */
      a[N0] = 1;
      b[N0] = 1;
  }

/*-----------------------------  ordinary dissipation  ----------*/
  if (DIS == 1) {
    for(k = 0; k<NTT; k++) {
      a[k] *= rk[k];
      b[k] *= rk[k];
    }
  }
/*----------------------- statistics --------------------------------------*/

    if (!(mt % IG)) {                 
      en = 0.0;  
      dis = 0.;
      for (k=0; k<=NN; k++) {
        eni[k] = (a[k]*a[k]+b[k]*b[k]);		
        ro[k] = sqrt(eni[k]);
        en += eni[k];
      }
      if (mt > MT0) {
        if (CORR==1) {
          for (j=1; j<5; j++){
            for (k=j; k<NN; k++) {
              teta = a[k]*(b[k-j]*a[k+1]+a[k-j]*b[k+1])+
                     b[k]*(a[k-j]*a[k+1]-b[k-j]*b[k+1]);
              set = teta/(ro[k-j]*ro[k]*ro[k+1]);
              tek[k][j] += teta/(double)NG;
              sek[k][j] += set/(double)NG;
              if (j==1) tet[k]=teta;
           }
          }
        }
        if (STRUCT==1){
          for (k=2; k<NN; k++){
            dis1 = (fl1*tet[k-1]+fl2*tet[k])*dnn[k];
            disi[k] = dis1;
            dis += disi[k];
            dis1 *= sign(dis1);

            for (j=1; j<JBENZI; j++){
              stru[j][k] += (pow(ro[k],(double)j))/(double)NG; 
              eps[j][k]  += (pow(dis1,(double)j))/(double)NG; 
            }
          }
        }
        for (k=0; k<=NN; k++) {                   /*   average:            */
	  enisr[k] += eni[k]/(double)NG;           /*            energy     */
	  dissr[k] += disi[k]/(double)NG;           /* enstrophy dissipation */
        }
        if (TIMESER==1){
          tim[igg] = tem;
          e[igg] = en;                                /*    total  energy    */
          disen[igg] = dis;
        }
        igg ++;
      }
      if (!(mt % 10000)){
         printf ("\n");	 
         printf ("%d, t=%f, en=%e dis=%e a[4]=%e a[6]=%e a[9]=%e",mt,tem,en,dis,a[4],a[6],a[9]);
       }
    }                         
  }                                                  /* end of main loop     */
  if (STRUCT==1){
    for (k=2; k<NN; k++){
      for (j=1; j<(JBENZI-1); j++){
        lstru[j][k] = log(stru[j][k])/log(2.0);
        epsrel[j][k] = eps[j+1][k]/eps[j][k];         
        zeta[j][k] = (lstru[j][k-1]-lstru[j][k]);
      }
      pirel[k] = eps[1][k]/epsrel[JBENZI-2][k];
    }

  }
/*------------------------------- output block ------------------------------*/

  prepout(fn2);
    fprintf (ff,"%e \n",tem);
    for (k=0; k<=NN; k++) fprintf (ff,"%e %e \n",a[k], b[k]);
  fclose(ff);

  printf ("\n End of output a(n) and b(n) !\n");

                                            /* output SPECTRA    */

  prepout(fn6);   
  for (k=1; k<NN; k++) {
    enisr[k] = lstru[2][k];
    sl = (enisr[k]-enisr[k-1]);
    fprintf (ff," %d  %e  %e %e \n", k-N0,enisr[k],sl,dissr[k] );
  }
  fclose (ff);



if (TIMESER==1){
    prepout(fn3);
    for (j=0; j<NG; j++) fprintf (ff,"  %f  %f\n", tim[j], e[j]);
    fprintf (ff,"\n");
    fclose(ff);

    prepout(fn5);
    for (j=0; j<NG; j++) fprintf (ff,"  %f  %f\n", tim[j], disen[j]);
    fclose(ff);
  }
  if (CORR==1){
/*    prepout(fn11);
    for (j=1; j<5; j++){
      for (k=j; k<NN; k++) {
        rek = tek[k][j]/(sqrt(enisr[k-j])*sqrt(enisr[k])*sqrt(enisr[k+1]));
        fprintf (ff," %d   %e\n", k-N0,rek );
      }
    fprintf (ff,"\n");
    }
    fclose(ff);
*/
    prepout(fn12);
    for (j=1; j<5; j++){
      for (k=j; k<NN; k++) {
        fprintf (ff," %d   %e\n", k-N0,sek[k][j]);
      }
    fprintf (ff,"\n");
    }
    fclose(ff);
  }
  if (STRUCT==1){

    prepout(fn13);
    for (j=1; j<JBENZI; j++){
      for (k=1; k<NTT; k++) {
        fprintf (ff," %d   %e\n", (k-N0), zeta[j][k]);
      }
    fprintf (ff,"\n");
    }
    fclose(ff);

    prepout(fn14);
    for (j=1; j<(JBENZI-1); j++){
      for (k=2; k<NN; k++) {
        fprintf (ff," %d   %e\n", (k-N0), (zeta[j][k]/zeta[3][k]));
      }
    fprintf (ff,"\n");
    }
    fclose(ff);

    prepout(fn16);
    for (j=1; j<JBENZI-1; j++){
      for (k=0; k<NTT; k++) {
        fprintf (ff," %d   %e\n", (k-N0), stru[j][k]);
      }
    fprintf (ff,"\n");
    }
    fclose(ff);

    prepout(fn19);
    for (j=1; j<JBENZI-1; j++){
      for (k=2; k<NN; k++) {
        fprintf (ff," %d   %e\n", (k-N0), log(epsrel[j][k])/log(2.0));
      }
    fprintf (ff,"\n");
    }
    fclose(ff);

    prepout(fn17);
    for (j=1; j<JBENZI-1; j++){
      for (k=2; k<NN; k++) {
        fprintf (ff," %d   %e\n", k-N0,eps[j][k]);
      }
    fprintf (ff,"\n");
    }
    fclose(ff);

    prepout(fn18);
    for (j=1; j<JBENZI-1; j++){
      for (k=2; k<NTT; k++) {
        fprintf (ff," %d   %e\n", k-N0,epsrel[j][k]);
      }
    fprintf (ff,"\n");
    }
    fclose(ff);

    prepout(fn20);
    for (k=2; k<NN; k++) {
      fprintf (ff," %e   %e\n", stru[3][k],pirel[k]);
    }
    fclose(ff);

    prepout(fn21);
    for (k=8; k<25; k++) {
      for (j=2; j<10; j++){

        fprintf (ff," %e   %e\n", (epsrel[j-1][k]/epsrel[JBENZI-2][k]),
                                  (epsrel[j][k]/epsrel[JBENZI-2][k]));
      }
    fprintf (ff,"\n");
    }
    fclose(ff);
  }
  printf ("\n End of GOY !\n\n\a\a");
}

/*=======================================================================*/

void sima(double *zza, double *zzb, double *x, double *y)
{
  int	i, j, k;
  double	ti, tj, tk;

  for(i =0; i<NTT; i++) {
    x[i] = 0.0;
    y[i] = 0.0;
  }
  for (k=1; k<NN; k++) {
    i = k - 1;
    j = k + 1;
    tk = dnn[k]*CK;
    ti = dnn[k]*CI;
    tj = dnn[k]*CJ;


    y[i] += ti*(zza[j]*zza[k]-zzb[j]*zzb[k]);
    y[j] += tj*(zza[i]*zza[k]-zzb[i]*zzb[k]);
    y[k] += tk*(zza[j]*zza[i]-zzb[j]*zzb[i]);

    x[i] += ti*(zza[j]*zzb[k]+zzb[j]*zza[k]);
    x[j] += tj*(zza[i]*zzb[k]+zzb[i]*zza[k]);
    x[k] += tk*(zza[j]*zzb[i]+zzb[j]*zza[i]);
  }

  for(i = 0; i<NTT; i++) {
    x[i] = TAU*x[i];
    y[i] = TAU*y[i];
  }
}
/*=======================================================================*/

double	sign(double x)
{  
  double ss;
  ss=1.0;
  if (x < 0.0) ss = -1.0;
  return ss;
}  
/*=======================================================================*/
void prepin(char *x)
{
  strcat(x,INTNAME);
  if((ff = fopen (x,"r+"))==NULL){
    printf ("Can't open file %s !\n\a\a",x);
    exit(-1);
  }
}
/*=======================================================================*/
void prepout(char *x)
{
  strcat(x,OUTNAME);
  if((ff = fopen (x,"w+"))==NULL){
    printf ("Can't open file %s !\n\a\a",x);
    exit(-1);
  }
}
/*=======================================================================*/





