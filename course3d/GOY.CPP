//----------------------- Calculation of Shell model ----------//
#define LEVEL extern
#include "head.h"


int N0=2;
double dtsh=1e-5;
double **shza,**shzb,**shx1,**shy1,**shx2,**shy2,**shx3,**shy3,**shx4,**shy4;
double EPS,CK,CJ,CI;
double *rk,*dnn;

void init_shell(void)
{
int i,j,k,nk;

sha=alloc_mem_2f(n3,Ns+1);
shb=alloc_mem_2f(n3,Ns+1);
shza=alloc_mem_2f(n3,Ns+1);
shzb=alloc_mem_2f(n3,Ns+1);
shx1=alloc_mem_2f(n3,Ns+1);
shy1=alloc_mem_2f(n3,Ns+1);
shx2=alloc_mem_2f(n3,Ns+1);
shy2=alloc_mem_2f(n3,Ns+1);
shx3=alloc_mem_2f(n3,Ns+1);
shy3=alloc_mem_2f(n3,Ns+1);
shx4=alloc_mem_2f(n3,Ns+1);
shy4=alloc_mem_2f(n3,Ns+1);
dnn=alloc_mem_1f(Ns+1);
rk=alloc_mem_1f(Ns+1);
nl =alloc_mem_1f(Ns+1);

  for (k=0; k<=Ns; k++) {
    nk = k - N0;
    dnn[k] = pow(2.0,(double)nk);
    nl[k] = dnn[k]*dnn[k];
    rk[k] = exp(-dtsh*nl[k]/Re);
  }
 for (j=0; j<n3; j++)
    for (i=0; i<=Ns; i++) {
         sha[j][i] = (i==0||i==1)? 1:0;
         shb[j][i] = (i==0||i==1)? 1.5:0;
    }
    EPS = -2.;
    CK = -EPS/2.0;
    CJ = (EPS-1.0)/2.0;
    CI = 0.5;

}

void sima(double *zza, double *zzb, double *xx, double *yy)
{
  int  i, j, k;
  double  ti, tj, tk;

  for(i =0; i<=Ns; i++) {
    xx[i] = 0.0;
    yy[i] = 0.0;
  }
  for (k=1; k<Ns; k++) {
    i = k - 1;
    j = k + 1;
    tk = dnn[k]*CK;
    ti = dnn[k]*CI;
    tj = dnn[k]*CJ;

    yy[i] += ti*(zza[j]*zza[k]-zzb[j]*zzb[k]);
    yy[j] += tj*(zza[i]*zza[k]-zzb[i]*zzb[k]);
    yy[k] += tk*(zza[j]*zza[i]-zzb[j]*zzb[i]);

    xx[i] += ti*(zza[j]*zzb[k]+zzb[j]*zza[k]);
    xx[j] += tj*(zza[i]*zzb[k]+zzb[i]*zza[k]);
    xx[k] += tk*(zza[j]*zzb[i]+zzb[j]*zza[i]);
  }

  for(i = 0; i<=Ns; i++) {
    xx[i] = dtsh*xx[i];
    yy[i] = dtsh*yy[i];
  }
}

void time_step_shell(double tm)
{
int i,j,k;
double tm_cur, en;

for (tm_cur=0;tm_cur<tm;tm_cur+=dtsh)
    enter++;
 for (j=0; j<n3; j++) {
    for (i=0,en=0; i<=Ns; i++)
      en+=sha[j][i]*sha[j][i]+shb[j][i]*shb[j][i];

    for (i=0; i<=Ns; i++) {
      shza[j][i] = sha[j][i];
      shzb[j][i] = shb[j][i];
    }
    sima(shza[j],shzb[j],shx1[j],shy1[j]);
    for (i=0; i<=Ns; i++) {
      shza[j][i] = sha[j][i]+shx1[j][i]/2.0f;
      shzb[j][i] = shb[j][i]+shy1[j][i]/2.0f;
    }
    sima(shza[j],shzb[j],shx2[j],shy2[j]);
    for (i=0; i<=Ns; i++) {
      shza[j][i] = sha[j][i]+shx2[j][i]/2.0f;
      shzb[j][i] = shb[j][i]+shy2[j][i]/2.0f;
    }
    sima(shza[j],shzb[j],shx3[j],shy3[j]);
     for (i=0; i<=Ns; i++) {
      shza[j][i] = sha[j][i]+shx3[j][i];
      shzb[j][i] = shb[j][i]+shy3[j][i];
    }
    sima(shza[j],shzb[j],shx4[j],shy4[j]);

    for (i=0; i<=Ns; i++){
      sha[j][i] += (shx1[j][i]+2.0f*shx2[j][i]+2.0f*shx3[j][i]+shx4[j][i])/6.0f;
      shb[j][i] += (shy1[j][i]+2.0f*shy2[j][i]+2.0f*shy3[j][i]+shy4[j][i])/6.0f;
    }

    for(k = 0; k<=Ns; k++) {
      sha[j][k] *= rk[k];
      shb[j][k] *= rk[k];
    }
  }
}


