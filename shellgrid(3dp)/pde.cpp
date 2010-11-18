//----------------------- Calculation of PDE right part  ----------//
#define LEVEL extern
#include "head.h"

inline double norma(double a,double b,double c,int order)
{
if(order==2) return ((a)*(a)+(b)*(b)+(c)*(c));
   else  return  pow(((a)*(a)+(b)*(b)+(c)*(c)),order/2.);
}

void pde(double t, double ****f, double ****df)
{
   int i,j,k,l,m;
   double dv1[3][3],dv2[3][3],dp1[3],dp2[3],dn1[3];

   boundary_conditions(f);

   for(i=ghost;i<mm1;i++)
   for(j=ghost;j<mm2;j++)
   for(k=ghost;k<mm3;k++) {
      for(l=0;l<3;l++)
      for(m=0;m<3;m++) {
         dv1[l][m]=dr(f[l],i,j,k,m+1,0,dx[m],ghost, approx);
         dv2[l][m]=dr(f[l],i,j,k,m+1,1,dx[m]*dx[m],ghost, approx);
      }
      for(m=0;m<3;m++) {
         dp1[m]=dr(f[3],i,j,k,m+1,0,dx[m],ghost, approx);
         dn1[m]=dr(nut,i,j,k,m+1,0,dx[m],ghost, approx);
      }

      for(l=0;l<3;l++)
       df[l][i][j][k]= (dv2[l][0] + dv2[l][1] + dv2[l][2])*nut[i][j][k]
                     - dp1[l] + (dn1[0]-f[0][i][j][k])*dv1[l][0]
                              + (dn1[1]-f[1][i][j][k])*dv1[l][1]
                              + (dn1[2]-f[2][i][j][k])*dv1[l][2];
      df[3][i][j][k]= (-(dv1[0][0] + dv1[1][1] + dv1[2][2]))/Gamma;
   }

   return;
}

double deviation(double ****f,int i,int j,int k)
{
const size_okr=min(1,ghost);
double flux;
int kol=0,l;
   for(l=1;l<=size_okr;l++)
            {
            flux = 0;
            if(i>0)
               {flux += norma(f[0][i][j][k]-f[0][i-l][j][k],f[1][i][j][k]-f[1][i-l][j][k],f[2][i][j][k]-f[2][i-l][j][k],2);
            	kol++;
               }
            if(i<m1)
               {flux += norma(f[0][i][j][k]-f[0][i+l][j][k],f[1][i][j][k]-f[1][i+l][j][k],f[2][i][j][k]-f[2][i+l][j][k],2);
            	kol++;
               }
            if(j>0)
               {flux += norma(f[0][i][j][k]-f[0][i][j-l][k],f[1][i][j][k]-f[1][i][j-l][k],f[2][i][j][k]-f[2][i][j-l][k],2);
            	kol++;
               }
            if(j<m2)
               {flux += norma(f[0][i][j][k]-f[0][i][j+l][k],f[1][i][j][k]-f[1][i][j+l][k],f[2][i][j][k]-f[2][i][j+l][k],2);
            	kol++;
               }
            if(k>0)
               {flux += norma(f[0][i][j][k]-f[0][i][j][k-l],f[1][i][j][k]-f[1][i][j][k-l],f[2][i][j][k]-f[2][i][j][k-l],2);
                kol++;
               }
            if(k<m3)
               {flux += norma(f[0][i][j][k]-f[0][i][j][k+l],f[1][i][j][k]-f[1][i][j][k+l],f[2][i][j][k]-f[2][i][j][k+l],2);
                kol++;
               }
           };
   flux /= kol;
return(flux);
}

void nut_by_flux(double ****f,double dt) //calculating nu_turbulent by velocity fluctuations
{
double maschtab = 100000;
int i,j,k,l;
double koef;
/*struct_func(f,2,2,3);
for(i=0;i<n3;i++)
    {
    koef=sqrt(s_func[i][0]/(pow(sha[i][1],2.)+pow(shb[i][1],2.)));
    sha[i][1] *= koef;
    shb[i][1] *= koef;
    koef=sqrt(s_func[i][1]/(pow(sha[i][0],2.)+pow(shb[i][0],2.)));
    sha[i][0] *= koef;
    shb[i][0] *= koef;
    }*/
/*clrscr();
for(j=0;j<n3;j++)
    {
    printf("%lf  %lf",s_func[j][0],s_func[j][1]);
    double en;
    for (i=0,en=0; i<=Ns; i++)
      en+=sha[j][i]*sha[j][i]+shb[j][i]*shb[j][i];
    printf("   totEn=%lf\n",en);
    }  */
//time_step_shell(dt);
for(k=0;k<n3;k++)
   {
/*   double tmp = maschtab*pow(
           (nl[2]*s_func[k][0] + nl[1]*s_func[k][1] + nl[0]*s_func[k][2])*pow(dx[2],4),
                         1./3);*/
   double tmp=0;
   for(i=ghost;i<mm1;i++)
       for(j=ghost;j<mm2;j++)
            nut[i][j][k+ghost] = (1. + tmp)/Re;
   }
}

void  boundary_conditions(double ****f)
{
   int i, j, k, l, g, req_numS=0, req_numR=0;

   // streamwise surfaces
      //periodic for velocities and gradient-periodic for pressure
 if (pr_neighbour[0]!=-1) {
   CopyGridToBuffer(f,buf_send[0],ghost,ghost,ghost,2*ghost-1,mm2-1,mm3-1);
   MPI_Isend(buf_send[0],buf_size[0],MPI_DOUBLE,pr_neighbour[0],tag,MPI_COMM_WORLD,&SendRequest[req_numS++]);
   MPI_Irecv(buf_recv[0],buf_size[0],MPI_DOUBLE,pr_neighbour[0],tag,MPI_COMM_WORLD,&RecvRequest[req_numR++]);
   } else {
    for(l=0;l<nvar;l++)
      for(j=ghost;j<mm2;j++)
       for(k=ghost;k<mm3;k++)
        for(g=0;g<ghost;g++)
           f[l][g][j][k] = f[l][n1+g][j][k] - ((l==3)?(p2-p1):0);
    CopyGridToBuffer(f,buf_recv[0],0,ghost,ghost,ghost-1,mm2-1,mm3-1);
   }

 if (pr_neighbour[1]!=-1) {
   CopyGridToBuffer(f,buf_send[1],n1,ghost,ghost,mm1-1,mm2-1,mm3-1);
   MPI_Isend(buf_send[1],buf_size[0],MPI_DOUBLE,pr_neighbour[1],tag,MPI_COMM_WORLD,&SendRequest[req_numS++]);
   MPI_Irecv(buf_recv[1],buf_size[0],MPI_DOUBLE,pr_neighbour[1],tag,MPI_COMM_WORLD,&RecvRequest[req_numR++]);
   } else {
    for(l=0;l<nvar;l++)
      for(j=ghost;j<mm2;j++)
       for(k=ghost;k<mm3;k++)
        for(g=0;g<ghost;g++)
          f[l][mm1+g][j][k] = f[l][ghost+g][j][k] + ((l==3)?(p2-p1):0);
    CopyGridToBuffer(f,buf_recv[1],mm1,ghost,ghost,m1-1,mm2-1,mm3-1);
   }

   // spanwise surfaces
       //periodic conditions for velocities and pressure
 if (pr_neighbour[2]!=-1) {
   CopyGridToBuffer(f,buf_send[2],ghost,ghost,ghost,mm1-1,2*ghost-1,mm3-1);
   MPI_Isend(buf_send[2],buf_size[1],MPI_DOUBLE,pr_neighbour[2],tag,MPI_COMM_WORLD,&SendRequest[req_numS++]);
   MPI_Irecv(buf_recv[2],buf_size[1],MPI_DOUBLE,pr_neighbour[2],tag,MPI_COMM_WORLD,&RecvRequest[req_numR++]);
   } else {
    for(l=0;l<nvar;l++)
     for(i=ghost;i<mm1;i++)
      for(k=ghost;k<mm3;k++)
       for(g=0;g<ghost;g++)
         f[l][i][g][k] = f[l][i][n2+g][k];
    CopyGridToBuffer(f,buf_recv[2],ghost,0,ghost,mm1-1,ghost-1,mm3-1);
   }

 if (pr_neighbour[3]!=-1) {
   CopyGridToBuffer(f,buf_send[3],ghost,n2,ghost,mm1-1,mm2-1,mm3-1);
   MPI_Isend(buf_send[3],buf_size[1],MPI_DOUBLE,pr_neighbour[3],tag,MPI_COMM_WORLD,&SendRequest[req_numS++]);
   MPI_Irecv(buf_recv[3],buf_size[1],MPI_DOUBLE,pr_neighbour[3],tag,MPI_COMM_WORLD,&RecvRequest[req_numR++]);
   } else {
    for(l=0;l<nvar;l++)
     for(i=ghost;i<mm1;i++)
      for(k=ghost;k<mm3;k++)
       for(g=0;g<ghost;g++)
        f[l][i][mm2+g][k] = f[l][i][ghost+g][k];
    CopyGridToBuffer(f,buf_recv[3],ghost,mm2,ghost,mm1-1,m2-1,mm3-1);
   }

   // on horizontal surfaces
      //sticking for velocities and free conditions for pressure
 if (pr_neighbour[4]!=-1) {
   CopyGridToBuffer(f,buf_send[4],ghost,ghost,ghost,mm1-1,mm2-1,2*ghost-1);
   MPI_Isend(buf_send[4],buf_size[2],MPI_DOUBLE,pr_neighbour[4],tag,MPI_COMM_WORLD,&SendRequest[req_numS++]);
   MPI_Irecv(buf_recv[4],buf_size[2],MPI_DOUBLE,pr_neighbour[4],tag,MPI_COMM_WORLD,&RecvRequest[req_numR++]);
   } else {
    for(l=0;l<nvar;l++)
     for(i=ghost;i<mm1;i++)
      for(j=ghost;j<mm2;j++)
       for(g=0;g<ghost;g++)
         f[l][i][j][g] = f[l][i][j][2*ghost-1-g]*((l==3)?1:-1);
    CopyGridToBuffer(f,buf_recv[4],ghost,ghost,0,mm1-1,mm2-1,ghost-1);
   }

 if (pr_neighbour[5]!=-1) {
   CopyGridToBuffer(f,buf_send[5],ghost,ghost,n3,mm1-1,mm2-1,mm3-1);
   MPI_Isend(buf_send[5],buf_size[2],MPI_DOUBLE,pr_neighbour[5],tag,MPI_COMM_WORLD,&SendRequest[req_numS++]);
   MPI_Irecv(buf_recv[5],buf_size[2],MPI_DOUBLE,pr_neighbour[5],tag,MPI_COMM_WORLD,&RecvRequest[req_numR++]);
   } else {
    for(l=0;l<nvar;l++)
     for(i=ghost;i<mm1;i++)
      for(j=ghost;j<mm2;j++)
       for(g=0;g<ghost;g++)
         f[l][i][j][mm3+g] = f[l][i][j][mm3-1-g]*((l==3)?1:-1);
    CopyGridToBuffer(f,buf_recv[5],ghost,ghost,mm3,mm1-1,mm2-1,m3-1);
   }

  MPI_Waitall(req_numR,RecvRequest,statuses);

    CopyBufferToGrid(f,buf_recv[0],0,ghost,ghost,ghost-1,mm2-1,mm3-1);
    CopyBufferToGrid(f,buf_recv[1],mm1,ghost,ghost,m1-1,mm2-1,mm3-1);
    CopyBufferToGrid(f,buf_recv[2],ghost,0,ghost,mm1-1,ghost-1,mm3-1);
    CopyBufferToGrid(f,buf_recv[3],ghost,mm2,ghost,mm1-1,m2-1,mm3-1);
    CopyBufferToGrid(f,buf_recv[4],ghost,ghost,0,mm1-1,mm2-1,ghost-1);
    CopyBufferToGrid(f,buf_recv[5],ghost,ghost,mm3,mm1-1,mm2-1,m3-1);

   return;
}

void  init_conditions(double ****f,double Re)
{
   int i,j,k,l;
   double Noise=0., Noise1=0;
//   double k1,k2,k3;

// --------------- initial conditions
//   k1=2*M_PI/l1;  k3=M_PI/l3;

   for(i=0;i<m1;i++)
   for(j=0;j<m2;j++)
   for(k=0;k<m3;k++) {
        f[0][i][j][k]=(0.5+Noise*((double)rand()-RAND_MAX/2)/RAND_MAX)*
                       coordin(k,2)*(l3-coordin(k,2))*4/l3/l3;
        f[1][i][j][k]=Noise1*cos(2*M_PI*coordin(j,1)/l2)*cos(2*M_PI*coordin(k,2)/l3)
                      + Noise*((double)rand()-RAND_MAX/2)/RAND_MAX*
                       coordin(k,2)*(l3-coordin(k,2))*4/l3/l3;
        f[2][i][j][k]=Noise1*sin(2*M_PI*coordin(j,1)/l2)*sin(2*M_PI*coordin(k,2)/l3)
                      + Noise*((double)rand()-RAND_MAX/2)/RAND_MAX*
                       coordin(k,2)*(l3-coordin(k,2))*4/l3/l3;
        f[3][i][j][k]=p1+(i-0.5)*(p2-p1)/n1;
        nut[i][j][k]=(
        (0.39+14.8*exp(-2.13*pow(2*coordin(k,2)-l3,2)))*0.1*0
                    +1)/Re;
   }
//   struct_func(f,2,2,3);
}

void init_parallel()
{
   int divisors[100],kd=0,nd1,nd2;
   int i,j;
   int kp1,kp2,kp3;
   int pr[3];
   double vtime,mintime;
   double k1=1,                   // time of calculation per node
          k2=1,                   // time of sending per node
          k3=1;                   // latent time per surface

//--------------meshing-----------------
  for(i=1;i<=size;i++)
    if(size%i==0) divisors[kd++] = i;
  mintime=-1;
  for(i=0;i<kd;i++)
    for(j=0;j<kd;j++)
       if(size%(divisors[i]*divisors[j])==0)
         {
         kp1=divisors[i]; kp2=divisors[j]; kp3=size/kp1/kp2;
         vtime = k1*ceil(double(N1)/kp1)*ceil(double(N2)/kp2)*ceil(double(N3)/kp3)+
                 k2*((kp1-1)*N2*N3+N1*(kp2-1)*N3+N1*N2*(kp3-1))+
                 k3*((kp1-1)*kp2*kp3+kp1*(kp2-1)*kp3+kp1*kp2*(kp3-1));
         if(mintime==-1 || vtime<mintime) { mintime=vtime; nd1=i; nd2=j; }
         }
  pp[0]=divisors[nd1]; pp[1]=divisors[nd2]; pp[2]=size/pp[0]/pp[1];                  // number of procs along axes
  n1 = ceil(double(N1)/pp[0]);         if(pr[0]>=N1-pp[0]*(n1-1)) n1--;              // dimensions of subregion
  n2 = ceil(double(N2)/pp[1]);         if(pr[0]>=N1-pp[0]*(n1-1)) n1--;
  n3 = ceil(double(N3)/pp[2]);         if(pr[2]>=N3-pp[2]*(n3-1)) n3--;          
  pr[0] = rank%pp[0]; pr[1] = (rank/pp[0])%pp[1]; pr[2] = (rank/pp[0]/pp[1])%pp[2];  // coordinates of subregion

   m1 = n1+2*ghost;
   m2 = n2+2*ghost;
   m3 = n3+2*ghost;
   mm1 = ghost+n1;
   mm2 = ghost+n2;
   mm3 = ghost+n3;

  pr_neighbour[0] = (pr[0]>0       ? rank-1           :-1);                           // neighbours of subregion
  pr_neighbour[1] = (pr[0]<pp[0]-1 ? rank+1           :-1);
  pr_neighbour[2] = (pr[1]>0       ? rank-pp[0]       :-1);
  pr_neighbour[3] = (pr[1]<pp[1]-1 ? rank+pp[0]       :-1);
  pr_neighbour[4] = (pr[2]>0       ? rank-pp[0]*pp[1] :-1);
  pr_neighbour[5] = (pr[2]<pp[2]-1 ? rank+pp[0]*pp[1] :-1);

 buf_size[0]=n2*n3*nvar*ghost;
 buf_size[1]=n1*n3*nvar*ghost;
 buf_size[2]=n1*n2*nvar*ghost;

 for(i=0;i<nvar-1;i++)
  for(j=0;j<=1;j++) {
   buf_send[j+2*i] = alloc_mem_1f(buf_size[i]);
   buf_recv[j+2*i] = alloc_mem_1f(buf_size[i]);
  }

   iop=fopen(fname_stat,"w");
   fprintf(iop,"%d\n",rank);
   fprintf(iop,"%d\t%d\t%d\n",pr[0],pr[1],pr[2]);
   fprintf(iop,"%d\t%d\t%d\n",pp[0],pp[1],pp[2]);
   for(i=0;i<6;i++)
     fprintf(iop,"%d ",pr_neighbour[i]);
   fprintf(iop,"\n");
   fclose(iop);

}

static double kf3[2][3][3]={{{-3./2.0, 2.0, -1./2.0}, {-1./2.0, 0.0, 1./2.0}, {1./2.0, -2.0, 3./2.0}},
                     {{1.0, -2.0, 1.0}, {1.0, -2.0, 1.0}, {1.0, -2.0, 1.0}}};
static double kf5[2][5][5]={{{-25./12.0, 4.0, -3.0, 4./3.0, -1./4.0}, {-1./4.0, -5./6.0, 3./2.0, -1./2.0, 1./12.0},
							{1./12.0, -2./3.0, 0.0, 2./3.0, -1./12.0}, {-1./12.0, 1./2.0, -3./2.0, 5./6.0, 1./4.0},
							{1./4.0, -4./3.0, 3.0, -4.0, 25./12.0}}, {{35./12.0, -26./3.0, 19./2.0, -14./3.0, 11./12.0},
                     {11./12.0, -5./3.0, 1./2.0, 1./3.0, -1./12.0}, {-1./12.0, 4./3.0, -5./2.0, 4./3.0, -1./12.0},
                     {-1./12.0, 1./3.0, 1./2.0, -5./3.0, 11./12.0}, {11./12.0, -14./3.0, 19./2.0, -26./3.0, 35./12.0}}};
static double kf7[2][7][7]={{{-49./20.0, 6.0, -15./2.0, 20./3.0, -15./4.0, 6./5.0, -1./6.0},	{-1./6.0, -77./60.0, 5./2.0, -5./3.0, 5./6.0,-1./4.0, 1./30.0},
                     {1./30.0, -2./5.0, -7./12.0, 4./3.0, -1./2.0, 2./15.0, -1./60.0}, {-1./60.0, 3./20.0, -3./4.0, 0.0, 3./4.0, -3./20.0, 1./60.0},
                     {1./60.0, -2./15.0, 1./2.0, -4./3.0, 7./12.0, 2./5.0, -1./30.0}, {-1./30.0, 1./4.0, -5./6.0, 5./3.0, -5./2.0, 77./60.0, 1./6.0},
							{1./6.0, -6./5.0, 15./4.0, -20./3.0, 15./2.0, -6.0, 49./20.0}}, {{203./45.0, -87./5.0, 117./4.0, -254./9.0, 33./2.0, -27./5.0, 137./180.0},
							{137./180.0, -49./60.0, -17./12.0, 47./18.0, -19./12.0, 31./60.0, -13./180.0}, {-13./180.0, 19./15.0, -7./3.0, 10./9.0, 1./12.0, -1./15.0, 1./90.0},
                     {1./90.0, -3./20.0, 3./2.0, -49./18.0, 3./2.0, -3./20.0, 1./90.0}, {1./90.0, -1./15.0, 1./12.0, 10./9.0, -7./3.0, 19./15.0, -13./180.0},
                     {-13./180.0, 31./60.0, -19./12.0, 47./18.0, -17./12.0, -49./60.0, 137./180.0}, {137./180.0, -27./5.0, 33./2.0, -254./9.0, 117./4.0, -87./5.0, 203./45.0}}};

double dr(double ***m, int ii, int jj, int kk, int dr, int or, double dx, int sh,  int sm)
/*        matrix     , point                 , direct, order , differ   , shift , sample */
/*                                           , 1,2,3 ,  0,1    dx,dx^2  , 0-left , 3,5,7 */
{
double tmp=0.0;
int i;

switch (sm*dr) {
	case 3 : for(i=0; i<sm; i++) tmp += m[ii+i-sh][jj][kk]*kf3[or][sh][i]; break;
	case 6 : for(i=0; i<sm; i++) tmp += m[ii][jj+i-sh][kk]*kf3[or][sh][i]; break;
	case 9 : for(i=0; i<sm; i++) tmp += m[ii][jj][kk+i-sh]*kf3[or][sh][i]; break;
	case 5 : for(i=0; i<sm; i++) tmp += m[ii+i-sh][jj][kk]*kf5[or][sh][i]; break;
	case 10: for(i=0; i<sm; i++) tmp += m[ii][jj+i-sh][kk]*kf5[or][sh][i]; break;
	case 15: for(i=0; i<sm; i++) tmp += m[ii][jj][kk+i-sh]*kf5[or][sh][i]; break;
	case 7 : for(i=0; i<sm; i++) tmp += m[ii+i-sh][jj][kk]*kf7[or][sh][i]; break;
	case 14: for(i=0; i<sm; i++) tmp += m[ii][jj+i-sh][kk]*kf7[or][sh][i]; break;
	case 21: for(i=0; i<sm; i++) tmp += m[ii][jj][kk+i-sh]*kf7[or][sh][i]; break;
	default :
    	nrerror("\nNO SUCH SAMPLE for derivative Bye ...",0);
	}
return(tmp/dx);
}

double coordin(int i, int dir)
                      //0-x,1-y,2-z
{
 return dx[dir]*((double)(i-ghost)+0.5);
}

 void CopyBufferToGrid(double ****m,double *buffer,int x1,int y1,int z1,int x2,int y2,int z2)
 {
   int i,j,k,l,n=0;
   for(l=0;l<nvar;l++)
    for(i=x1;i<=x2;i++)
     for(j=y1;j<=y2;j++)
      for(k=z1;k<=z2;k++)
      m[l][i][j][k]=buffer[n++];
 }

 void CopyGridToBuffer(double ****m,double *buffer,int x1,int y1,int z1,int x2,int y2,int z2)
 {
   int i,j,k,l,n=0;
   for(l=0;l<nvar;l++)
    for(i=x1;i<=x2;i++)
     for(j=y1;j<=y2;j++)
      for(k=z1;k<=z2;k++)
      buffer[n++]=m[l][i][j][k];
 }


