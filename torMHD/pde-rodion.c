//----------------------- Calculation of PDE right part  ----------//
#define LEVEL extern
#include "head.h"

#define Cos cos
#define Sin sin



void pde(double t, double ****f, double ****df)
{
   int i,j,k,l,m,n,lg;
   double dv1[7][3],dv2[7][3],dp1[3],tmp[7],tmp1[7];


   proc_communication(f);



   for(i=ghost;i<mm1;i++)
   for(j=ghost;j<mm2;j++)
   for(k=ghost;k<mm3;k++)
   if(cells[i][j][k]==0)
   {
	   for(n=-1;n<=2;n+=2)
		   for(lg=0,l=n;n*l<=ghost;l+=n) 
			if(cells[i+l][j][k]==1) {
			if(lg==0) lg=l;
			for(m=0;m<4;m++)
				   f[m+3][i+l][j][k]=(m==0?1:-1)*f[m+3][i-(l-2*lg+n)][j][k];
		   }

	   for(n=-1;n<=2;n+=2)
		   for(lg=0,l=n;n*l<=ghost;l+=n) 
			if(cells[i][j+l][k]==1) {
			if(lg==0) lg=l;
			for(m=0;m<4;m++)
				   f[m+3][i][j+l][k]=(m==0?1:-1)*f[m+3][i][j-(l-2*lg+n)][k];
		   }

	   for(n=-1;n<=2;n+=2)
		   for(lg=0,l=n;n*l<=ghost;l+=n) 
			if(cells[i][j][k+l]==1) {
			if(lg==0) lg=l;
			for(m=0;m<4;m++)
				   f[m+3][i][j][k+l]=(m==0?1:-1)*f[m+3][i][j][k-(l-2*lg+n)];
		   }


       for(m=0;m<3;m++)
         dp1[m]=dr(f[3],i,j,k,m+1,0,dx[m],ghost, approx);

      for(l=1;l<=3;l++) {
       for(m=0;m<3;m++) {
         dv1[l][m]=dr(f[l+3],i,j,k,m+1,0,dx[m],ghost, approx);
         dv2[l][m]=dr(f[l+3],i,j,k,m+1,1,dx[m]*dx[m],ghost, approx);
         }
      }

      df[4][i][j][k]=(dv2[1][0]+dv2[1][1]+dv2[1][2])/Re
		     -dp1[0]    +2*W*f[5][i][j][k]
		     -f[4][i][j][k]*dv1[1][0]
		     -f[5][i][j][k]*dv1[1][1]
		     -f[6][i][j][k]*dv1[1][2]
		     ;
      df[5][i][j][k]=(dv2[2][0]+dv2[2][1]+dv2[2][2])/Re
		     -dp1[1]    -2*W*f[4][i][j][k]
		     -f[4][i][j][k]*dv1[2][0]
		     -f[5][i][j][k]*dv1[2][1]
		     -f[6][i][j][k]*dv1[2][2]
		     ;
      df[6][i][j][k]=(dv2[3][0]+dv2[3][1]+dv2[3][2])/Re
		     -dp1[2]
		     -f[4][i][j][k]*dv1[3][0]
		     -f[5][i][j][k]*dv1[3][1]
		     -f[6][i][j][k]*dv1[3][2];
      df[3][i][j][k]= -(dv1[1][0]+dv1[2][1]+dv1[3][2])/Gamma;

      for(l=1;l<=3;l++)								//forcing
       df[l+3][i][j][k]+=R0*V[l-1][i][j][k];
   } else
   for(l=0;l<=3;l++)
      df[l+3][i][j][k]=0;

    calculate_curl(f,B);

   for(i=ghost;i<mm1;i++)
   for(j=ghost;j<mm2;j++)
   for(k=ghost;k<mm3;k++)
   {
	if(cells[i][j][k]==1)
     for(l=0;l<=3;l++)
      f[l+3][i][j][k]=0;

	for(l=0;l<=2;l++)
       for(m=0;m<3;m++) {
         dv2[l][m]=dr(f[l],i,j,k,m+1,1,dx[m]*dx[m],ghost, approx);
       }

       df[0][i][j][k]=(f[5][i][j][k]*B[2][i][j][k]-f[6][i][j][k]*B[1][i][j][k])
                     +(dv2[0][0]+dv2[0][1]+dv2[0][2])/Rm;
       df[1][i][j][k]=(f[6][i][j][k]*B[0][i][j][k]-f[4][i][j][k]*B[2][i][j][k])
                     +(dv2[1][0]+dv2[1][1]+dv2[1][2])/Rm;
       df[2][i][j][k]=(f[4][i][j][k]*B[1][i][j][k]-f[5][i][j][k]*B[0][i][j][k])
                     +(dv2[2][0]+dv2[2][1]+dv2[2][2])/Rm;

	  df[4][i][j][k]-=
		     -(dv2[2][0]+dv2[2][1]+dv2[2][2])*B[1][i][j][k]
		     +(dv2[1][0]+dv2[1][1]+dv2[1][2])*B[2][i][j][k]
		     ;
      df[5][i][j][k]-=
		     +(dv2[2][0]+dv2[2][1]+dv2[2][2])*B[0][i][j][k]
		     -(dv2[0][0]+dv2[0][1]+dv2[0][2])*B[2][i][j][k]
		     ;
      df[6][i][j][k]-=
		     -(dv2[1][0]+dv2[1][1]+dv2[1][2])*B[0][i][j][k]
		     +(dv2[0][0]+dv2[0][1]+dv2[0][2])*B[1][i][j][k]
		     ;

   }

   
   pdecount++;
   return;
}

void normal(double ****f)
{
double tmp, Bmax1, Bmax;
int i,j,k,l;

Bmax=0;

   for(i=ghost;i<mm1;i++)
   for(j=ghost;j<mm2;j++)
   for(k=ghost;k<mm3;k++)
        {
         tmp=sqrt(pow(f[0][i][j][k],2.)+pow(f[1][i][j][k],2.)+pow(f[2][i][j][k],2.));
         Bmax=(tmp>Bmax)?tmp:Bmax;
        }

 putLog("Start normal");

MPI_Allreduce(&Bmax, &Bmax1, 1, MPI_DOUBLE , MPI_MAX, MPI_COMM_WORLD);

 putLog("finis normal");

   for(i=ghost;i<mm1;i++)
   for(j=ghost;j<mm2;j++)
   for(k=ghost;k<mm3;k++)
        {
         f[0][i][j][k]/=Bmax1;
         f[1][i][j][k]/=Bmax1;
         f[2][i][j][k]/=Bmax1;
        }
}

void proc_communication(double ****f)
{
   int i,tag=10;
   int  req_numS,req_numR;
   int reslen;

  timeE0=MPI_Wtime();

  req_numS=req_numR=0;

  // exchanging in 1-direction
 if(pr_neighbour[2]>-1)
  if(pr_neighbour[2]==rank) CopyGridToBuffer(f,buf_recv[2],0,n2,0,m1-1,mm2-1,m3-1);
         else { CopyGridToBuffer(f,buf_send[2],0,ghost,0,m1-1,2*ghost-1,m3-1);
                MPI_Isend(buf_send[2],buf_size[1],MPI_DOUBLE,pr_neighbour[2],tag+2,MPI_COMM_WORLD,&SendRequest[req_numS++]);
                MPI_Irecv(buf_recv[2],buf_size[1],MPI_DOUBLE,pr_neighbour[2],tag+3,MPI_COMM_WORLD,&RecvRequest[req_numR++]);
                }
 if(pr_neighbour[3]>-1)
  if(pr_neighbour[3]==rank) CopyGridToBuffer(f,buf_recv[3],0,ghost,0,m1-1,2*ghost-1,m3-1);
         else { CopyGridToBuffer(f,buf_send[3],0,n2,0,m1-1,mm2-1,m3-1);
                MPI_Isend(buf_send[3],buf_size[1],MPI_DOUBLE,pr_neighbour[3],tag+3,MPI_COMM_WORLD,&SendRequest[req_numS++]);
                MPI_Irecv(buf_recv[3],buf_size[1],MPI_DOUBLE,pr_neighbour[3],tag+2,MPI_COMM_WORLD,&RecvRequest[req_numR++]);
                }


   MPI_Waitall(req_numR,RecvRequest,statuses);

  for(i=0;i<req_numR;i++)
   if(statuses[i].MPI_ERROR) {putLog("bc:error during transfer R1=");
                               MPI_Error_string(statuses[i].MPI_ERROR,msg,&reslen);
                               msg[reslen++] = ','; msg[reslen]= 0;
                               putLog(msg);
                               }

  if(pr_neighbour[2]>-1) CopyBufferToGrid(f,buf_recv[2],0,0,0,m1-1,ghost-1,m3-1);
  if(pr_neighbour[3]>-1) CopyBufferToGrid(f,buf_recv[3],0,mm2,0,m1-1,m2-1,m3-1);

   MPI_Waitall(req_numS,SendRequest,statuses);

  for(i=0;i<req_numS;i++)
   if(statuses[i].MPI_ERROR) {putLog("bc:error during transfer S1=");
                               MPI_Error_string(statuses[i].MPI_ERROR,msg,&reslen);
                               msg[reslen++] = ','; msg[reslen]= 0;
                               putLog(msg);
                               }

  req_numS=req_numR=0;

// exchanging in 0-direction
 if(pr_neighbour[0]>-1)
  if(pr_neighbour[0]==rank) CopyGridToBuffer(f,buf_recv[0],n1,0,0,mm1-1,m2-1,m3-1);
         else { CopyGridToBuffer(f,buf_send[0],ghost,0,0,2*ghost-1,m2-1,m3-1);
                MPI_Isend(buf_send[0],buf_size[0],MPI_DOUBLE,pr_neighbour[0],tag,MPI_COMM_WORLD,&SendRequest[req_numS++]);
                MPI_Irecv(buf_recv[0],buf_size[0],MPI_DOUBLE,pr_neighbour[0],tag+1,MPI_COMM_WORLD,&RecvRequest[req_numR++]);
                }
 if(pr_neighbour[1]>-1)
  if(pr_neighbour[1]==rank) CopyGridToBuffer(f,buf_recv[1],ghost,0,0,2*ghost-1,m2-1,m3-1);
         else { CopyGridToBuffer(f,buf_send[1],n1,0,0,mm1-1,m2-1,m3-1);
                MPI_Isend(buf_send[1],buf_size[0],MPI_DOUBLE,pr_neighbour[1],tag+1,MPI_COMM_WORLD,&SendRequest[req_numS++]);
                MPI_Irecv(buf_recv[1],buf_size[0],MPI_DOUBLE,pr_neighbour[1],tag,MPI_COMM_WORLD,&RecvRequest[req_numR++]);
                }

   MPI_Waitall(req_numR,RecvRequest,statuses);

  for(i=0;i<req_numR;i++)
   if(statuses[i].MPI_ERROR) {putLog("bc:error during transfer R0=");
                               MPI_Error_string(statuses[i].MPI_ERROR,msg,&reslen);
                               msg[reslen++] = ','; msg[reslen]= 0;
                               putLog(msg);
                               }

  if(pr_neighbour[0]>-1) CopyBufferToGrid(f,buf_recv[0],0,0,0,ghost-1,m2-1,m3-1);
  if(pr_neighbour[1]>-1) CopyBufferToGrid(f,buf_recv[1],mm1,0,0,m1-1,m2-1,m3-1);

  MPI_Waitall(req_numS,SendRequest,statuses);

  for(i=0;i<req_numS;i++)
   if(statuses[i].MPI_ERROR) {putLog("bc:error during transfer S0=");
                               MPI_Error_string(statuses[i].MPI_ERROR,msg,&reslen);
                               msg[reslen++] = ','; msg[reslen]= 0;
                               putLog(msg);
                               }

  req_numS=req_numR=0;

// exchanging in z-direction
 if(pr_neighbour[4]>-1)
  if(pr_neighbour[4]==rank) CopyGridToBuffer(f,buf_recv[4],0,0,n3,m1-1,m2-1,mm3-1);
         else { CopyGridToBuffer(f,buf_send[4],0,0,ghost,m1-1,m2-1,2*ghost-1);
                MPI_Isend(buf_send[4],buf_size[2],MPI_DOUBLE,pr_neighbour[4],tag+4,MPI_COMM_WORLD,&SendRequest[req_numS++]);
		MPI_Irecv(buf_recv[4],buf_size[2],MPI_DOUBLE,pr_neighbour[4],tag+5,MPI_COMM_WORLD,&RecvRequest[req_numR++]);
                }
 if(pr_neighbour[5]>-1)
  if(pr_neighbour[5]==rank) CopyGridToBuffer(f,buf_recv[5],0,0,ghost,m1-1,m2-1,2*ghost-1);
        else { CopyGridToBuffer(f,buf_send[5],0,0,n3,m1-1,m2-1,mm3-1);
               MPI_Isend(buf_send[5],buf_size[2],MPI_DOUBLE,pr_neighbour[5],tag+5,MPI_COMM_WORLD,&SendRequest[req_numS++]);
               MPI_Irecv(buf_recv[5],buf_size[2],MPI_DOUBLE,pr_neighbour[5],tag+4,MPI_COMM_WORLD,&RecvRequest[req_numR++]);
               }

   MPI_Waitall(req_numR,RecvRequest,statuses);

  for(i=0;i<req_numR;i++)
   if(statuses[i].MPI_ERROR) {putLog("bc:error during transfer R2=");
                               MPI_Error_string(statuses[i].MPI_ERROR,msg,&reslen);
                               msg[reslen++] = ','; msg[reslen]= 0;
                               putLog(msg);
                               }

  if(pr_neighbour[4]>-1) CopyBufferToGrid(f,buf_recv[4],0,0,0,m1-1,m2-1,ghost-1);
  if(pr_neighbour[5]>-1) CopyBufferToGrid(f,buf_recv[5],0,0,mm3,m1-1,m2-1,m3-1);

   MPI_Waitall(req_numS,SendRequest,statuses);

  for(i=0;i<req_numS;i++)
   if(statuses[i].MPI_ERROR) {putLog("bc:error during transfer S2=");
                               MPI_Error_string(statuses[i].MPI_ERROR,msg,&reslen);
                               msg[reslen++] = ','; msg[reslen]= 0;
                               putLog(msg);
                               }

  timeE1+=(MPI_Wtime()-timeE0);

  return;
}


void  init_conditions()
{
   int i,j,k,n;
   double r, rho, r1, z1, x1, y1, u1, v1, x0, y0, z0, tmp, tmpx,tmpy,tmpz;
   double rms;

   for(i=0;i<m1;i++)
   for(j=0;j<m2;j++)
   for(k=0;k<m3;k++)
     {
       tmpx=coordin(i,0);
       tmpy=coordin(j,1);
       tmpz=coordin(k,2);
					tmp=1.e-4;
                     f[0][i][j][k]=tmp*((double)rand()-RAND_MAX/2)/RAND_MAX;
                     f[1][i][j][k]=tmp*((double)rand()-RAND_MAX/2)/RAND_MAX;
                     f[2][i][j][k]=tmp*((double)rand()-RAND_MAX/2)/RAND_MAX;


                     f[3][i][j][k]=0;
                     f[4][i][j][k]=0;
                     f[5][i][j][k]=0;
					 f[6][i][j][k]=0;

       r1=sqrt(pow(tmpx,2.)+pow(tmpy,2.));
       r=sqrt(pow(r1-R0,2.)+pow(tmpz,2.));
       u1=atan2(tmpy,tmpx+1.e-10);
       v1=atan2(tmpz,r1-R0+1.e-10);
       if(r<=a0) {
        tmp=(fabs(u1+M_PI)<M_PI/15 || fabs(u1-M_PI)<M_PI/15 || fabs(u1)<M_PI/15 || fabs(u1-M_PI/2)<M_PI/15 || fabs(u1+M_PI/2)<M_PI/15?1:0);
		V[0][i][j][k]=-Hksi*r1*sin(u1)/R0+tmp*cos(u1)*r*sin(v1)/(1+a0/R0*cos(v1));
        V[1][i][j][k]= Hksi*r1*cos(u1)/R0+tmp*sin(u1)*r*sin(v1)/(1+a0/R0*cos(v1));
        V[2][i][j][k]=-tmp*r*cos(v1)/(1+a0/R0*cos(v1));
        cells[i][j][k]=0;

       } else {
        V[0][i][j][k]=V[1][i][j][k]=V[2][i][j][k]=0;

        cells[i][j][k]=1;
      }
      }
}

void init_parallel(int argc, char** argv)
{
   FILE *iop;
   int divisors[100],kd=0,nd1,nd2;
   int i,j;
   int kp1,kp2,kp3;
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
         vtime = k1*ceil((double)N1/kp1)*ceil((double)N2/kp2)*ceil((double)N3/kp3)+
                 k2*((kp1-1)*N2*N3+N1*(kp2-1)*N3+N1*N2*(kp3-1))+
                 k3*((kp1-1)*kp2*kp3+kp1*(kp2-1)*kp3+kp1*kp2*(kp3-1));
         if((mintime<0 || vtime<mintime) && ceil((double)N1/kp1)==floor((double)N1/kp1) && ceil((double)N2/kp2)==floor((double)N2/kp2) && ceil((double)N3/kp3)==floor((double)N3/kp3)) { mintime=vtime; nd1=i; nd2=j; }
         }
 if(argc < 5) {
  pp[0]=divisors[nd1]; pp[1]=divisors[nd2]; pp[2]=size/pp[0]/pp[1];                 // number of procs along axes
 } else {
  pp[0]=atoi(argv[2]);
  pp[1]=atoi(argv[3]);
  pp[2]=atoi(argv[4]);
  if (pp[0]*pp[1]*pp[2]!=size) putLog("------------- Wrong division of mesh: size is not p1*p2*p3 --------------");
  if (N1%pp[0]!=0 || N2%pp[1]!=0 || N3%pp[2]!=0) putLog("------------- Wrong division of mesh: n is not multiple p --------------");
 }
  pr[0] = rank%pp[0]; pr[1] = (rank/pp[0])%pp[1]; pr[2] = (rank/pp[0]/pp[1])%pp[2];  // coordinates of current subregion

/* dimensions of subregion:         global indicies of subregion origin          if there's nonequal subregions*/
  n1 = floor((double)N1/pp[0]);     n[0] = n1*pr[0] + min(pr[0],N1-pp[0]*n1);    if(pr[0]<N1-pp[0]*n1) n1++;              // dimensions of subregion
  n2 = floor((double)N2/pp[1]);     n[1] = n2*pr[1] + min(pr[1],N2-pp[1]*n2);    if(pr[1]<N2-pp[1]*n2) n2++;
  n3 = floor((double)N3/pp[2]);     n[2] = n3*pr[2] + min(pr[2],N3-pp[2]*n3);    if(pr[2]<N3-pp[2]*n3) n3++;
   if(n1<ghost || n2<ghost || n3<ghost) putLog("Too small mesh or incorrect number of processes");

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


  if(pr[0]==0)             pr_neighbour[0] = rank+pp[0]-1;     //periodic conditions
  if(pr[0]==pp[0]-1)       pr_neighbour[1] = rank-(pp[0]-1);     //periodic conditions
  if(pr[1]==0)             pr_neighbour[2] = rank+pp[0]*(pp[1]-1);     //periodic conditions
  if(pr[1]==pp[1]-1)       pr_neighbour[3] = rank-pp[0]*(pp[1]-1);     //periodic conditions
  if(pr[2]==0)             pr_neighbour[4] = rank+pp[0]*pp[1]*(pp[2]-1);     //periodic conditions
  if(pr[2]==pp[2]-1)       pr_neighbour[5] = rank-pp[0]*pp[1]*(pp[2]-1);     //periodic conditions


 buf_size[0]=m2*m3*(nvar)*ghost;
 buf_size[1]=m1*m3*(nvar)*ghost;
 buf_size[2]=m1*m2*(nvar)*ghost;

 for(i=0;i<3;i++)           //3-D
  for(j=0;j<=1;j++) {
   buf_send[j+2*i] = alloc_mem_1f(buf_size[i]);
   buf_recv[j+2*i] = alloc_mem_1f(buf_size[i]);
  }

   printf("%d\t",rank);
   printf("%d\t%d\t%d\t",pr[0],pr[1],pr[2]);
   printf("%d\t%d\t%d\t",pp[0],pp[1],pp[2]);
   printf("%d\t%d\t%d\t",n[0],n[1],n[2]);
   printf("%d\t%d\t%d\t",n1,n2,n3);
   for(i=0;i<6;i++)
     printf("%d\t",pr_neighbour[i]);
   printf("\n");

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

double dr(double ***m, int ii, int jj, int kk, int dir, int or, double dx, int sh,  int sm)
/*        matrix     , point                 , direct, order , differ   , shift , sample */
/*                                           , 1,2,3 ,  0,1    dx,dx^2  , 0-left , 3,5,7 */
{
double tmp=0.0;
int i;
//start_tick(9);
if(or==0)
switch (sm) {
     case 7 :  switch (dir) {
                  case 1 : tmp = kf7[or][sh][6]*(m[ii+3][jj][kk]-m[ii-3][jj][kk])
                               + kf7[or][sh][5]*(m[ii+2][jj][kk]-m[ii-2][jj][kk])
                               + kf7[or][sh][4]*(m[ii+1][jj][kk]-m[ii-1][jj][kk]); break;
                  case 2 : tmp = kf7[or][sh][6]*(m[ii][jj+3][kk]-m[ii][jj-3][kk])
                               + kf7[or][sh][5]*(m[ii][jj+2][kk]-m[ii][jj-2][kk])
                               + kf7[or][sh][4]*(m[ii][jj+1][kk]-m[ii][jj-1][kk]); break;
                  case 3 : tmp = kf7[or][sh][6]*(m[ii][jj][kk+3]-m[ii][jj][kk-3])
                               + kf7[or][sh][5]*(m[ii][jj][kk+2]-m[ii][jj][kk-2])
                               + kf7[or][sh][4]*(m[ii][jj][kk+1]-m[ii][jj][kk-1]); break;
                  }; break;
     case 3 :  switch (dir) {
                  case 1 : tmp = kf3[or][sh][2]*(m[ii+1][jj][kk]-m[ii-1][jj][kk]); break;
                  case 2 : tmp = kf3[or][sh][2]*(m[ii][jj+1][kk]-m[ii][jj-1][kk]); break;
                  case 3 : tmp = kf3[or][sh][2]*(m[ii][jj][kk+1]-m[ii][jj][kk-1]); break;
                  }; break;
     case 5 :  switch (dir) {
                  case 1 : tmp = kf5[or][sh][4]*(m[ii+2][jj][kk]-m[ii-2][jj][kk])
                               + kf5[or][sh][3]*(m[ii+1][jj][kk]-m[ii-1][jj][kk]); break;
                  case 2 : tmp = kf5[or][sh][4]*(m[ii][jj+2][kk]-m[ii][jj-2][kk])
                               + kf5[or][sh][3]*(m[ii][jj+1][kk]-m[ii][jj-1][kk]); break;
                  case 3 : tmp = kf5[or][sh][4]*(m[ii][jj][kk+2]-m[ii][jj][kk-2])
                               + kf5[or][sh][3]*(m[ii][jj][kk+1]-m[ii][jj][kk-1]); break;
                  }; break;
    }
if(or==1)
switch (sm) {
     case 7 :  switch (dir) {
                  case 1 : tmp = kf7[or][sh][6]*(m[ii+3][jj][kk]+m[ii-3][jj][kk])
                               + kf7[or][sh][5]*(m[ii+2][jj][kk]+m[ii-2][jj][kk])
                               + kf7[or][sh][4]*(m[ii+1][jj][kk]+m[ii-1][jj][kk])
                               + kf7[or][sh][3]*m[ii][jj][kk]; break;
                  case 2 : tmp = kf7[or][sh][6]*(m[ii][jj+3][kk]+m[ii][jj-3][kk])
                               + kf7[or][sh][5]*(m[ii][jj+2][kk]+m[ii][jj-2][kk])
                               + kf7[or][sh][4]*(m[ii][jj+1][kk]+m[ii][jj-1][kk])
                               + kf7[or][sh][3]*m[ii][jj][kk]; break;
                  case 3 : tmp = kf7[or][sh][6]*(m[ii][jj][kk+3]+m[ii][jj][kk-3])
                               + kf7[or][sh][5]*(m[ii][jj][kk+2]+m[ii][jj][kk-2])
                               + kf7[or][sh][4]*(m[ii][jj][kk+1]+m[ii][jj][kk-1])
                               + kf7[or][sh][3]*m[ii][jj][kk]; break;
                  }; break;
     case 3 :  switch (dir) {
                  case 1 : tmp = kf3[or][sh][2]*(m[ii+1][jj][kk]+m[ii-1][jj][kk])+ kf3[or][sh][1]*m[ii][jj][kk]; break;
                  case 2 : tmp = kf3[or][sh][2]*(m[ii][jj+1][kk]+m[ii][jj-1][kk])+ kf3[or][sh][1]*m[ii][jj][kk]; break;
                  case 3 : tmp = kf3[or][sh][2]*(m[ii][jj][kk+1]+m[ii][jj][kk-1])+ kf3[or][sh][1]*m[ii][jj][kk]; break;
                  }; break;
     case 5 :  switch (dir) {
                  case 1 : tmp = kf5[or][sh][4]*(m[ii+2][jj][kk]+m[ii-2][jj][kk])
                               + kf5[or][sh][3]*(m[ii+1][jj][kk]+m[ii-1][jj][kk])
                               + kf5[or][sh][2]*m[ii][jj][kk]; break;
                  case 2 : tmp = kf5[or][sh][4]*(m[ii][jj+2][kk]+m[ii][jj-2][kk])
                               + kf5[or][sh][3]*(m[ii][jj+1][kk]+m[ii][jj-1][kk])
                               + kf5[or][sh][2]*m[ii][jj][kk]; break;
                  case 3 : tmp = kf5[or][sh][4]*(m[ii][jj][kk+2]+m[ii][jj][kk-2])
                               + kf5[or][sh][3]*(m[ii][jj][kk+1]+m[ii][jj][kk-1])
                               + kf5[or][sh][2]*m[ii][jj][kk]; break;
                  }; break;
    }
/*switch (sm*dir) {
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
    	nrerror("\nNO SUCH SAMPLE for derivative. Bye ...",0,0);
	} */
return(tmp/dx);
}

double coordin(int i, int dir)
		      //0-r,1-phi,2-z
{
 switch (dir)
 {
   case 0:  return dx[dir]*(i-ghost+0.5+n[dir])-LX;
   case 1:  return dx[dir]*(i-ghost+0.5+n[dir])-LY;
   case 2:  return dx[dir]*(i-ghost+0.5+n[dir])-LZ;
 }
    return(0);
}

void calculate_curl(double ****f,double ****b)
{
   int i,j,k,l,m,n;
   double dA[7][3];

   for(i=ghost;i<mm1;i++)
   for(j=ghost;j<mm2;j++)
   for(k=ghost;k<mm3;k++)
   {
      for(l=0;l<=2;l++) {
       for(m=0;m<3;m++)
         dA[l][m]=dr(f[l],i,j,k,m+1,0,dx[m],ghost, approx);
       }
      b[0][i][j][k] = dA[2][1] - dA[1][2];
      b[1][i][j][k] = dA[0][2] - dA[2][0];
      b[2][i][j][k] = dA[1][0] - dA[0][1];
    }
}

void  communication_test()
{
   int i,j,k,n;
   double tmp, tmpx,tmpy,tmpz;
   double rms;

   for(i=0;i<m1;i++)
   for(j=0;j<m2;j++)
   for(k=0;k<m3;k++) {
       tmpx=coordin(i,0);
       tmpy=coordin(j,1);
       tmpz=coordin(k,2);

	   for(n=0;n<nvar;n++) {
		   df[n][i][j][k]=0;
		   f[n][i][j][k]=(n+1)*cos(2*M_PI*tmpx/LX+n)*cos(3*M_PI*tmpy/LY+4-n)*cos(4*M_PI*tmpz/LZ+1+n);
	   }
     }
   for(i=ghost;i<mm1;i++)
   for(j=ghost;j<mm2;j++)
   for(k=ghost;k<mm3;k++) {

	   for(n=0;n<nvar;n++) {
		   df[n][i][j][k]=f[n][i][j][k];
	   }
     }

  if(pr[0]==0)             pr_neighbour[0] = rank+pp[0]-1;     //periodic conditions
  if(pr[0]==pp[0]-1)       pr_neighbour[1] = rank-(pp[0]-1);     //periodic conditions
  if(pr[1]==0)             pr_neighbour[2] = rank+pp[0]*(pp[1]-1);     //periodic conditions
  if(pr[1]==pp[1]-1)       pr_neighbour[3] = rank-pp[0]*(pp[1]-1);     //periodic conditions
  if(pr[2]==0)             pr_neighbour[4] = rank+pp[0]*pp[1]*(pp[2]-1);     //periodic conditions
  if(pr[2]==pp[2]-1)       pr_neighbour[5] = rank-pp[0]*pp[1]*(pp[2]-1);     //periodic conditions

   timeE0=MPI_Wtime();
   proc_communication(df);
   timeE0=MPI_Wtime()-timeE0;

   rms=0;
   for(i=0;i<m1;i++)
   for(j=0;j<m2;j++)
   for(k=0;k<m3;k++) {
    tmp=0;
	for(n=0;n<nvar;n++)
	 tmp+=fabs(df[n][i][j][k]-f[n][i][j][k]);
    if(tmp>rms) 
		rms=tmp;
   }

  MPI_Allreduce(&rms, &tmp, 1, MPI_DOUBLE , MPI_MAX, MPI_COMM_WORLD);
   if(tmp<1.e-10)
    putLog("----------- Communication test OK --------------");
   else
    putLog("----------- Communication test failed! --------------");
  sprintf(msg,"----------- Communication time = %f ----------",timeE0);
  putLog(msg);

  pr_neighbour[0] = (pr[0]>0       ? rank-1           :-1);
  pr_neighbour[1] = (pr[0]<pp[0]-1 ? rank+1           :-1);
  pr_neighbour[2] = (pr[1]>0       ? rank-pp[0]       :-1);
  pr_neighbour[3] = (pr[1]<pp[1]-1 ? rank+pp[0]       :-1);
  pr_neighbour[4] = (pr[2]>0       ? rank-pp[0]*pp[1] :-1);
  pr_neighbour[5] = (pr[2]<pp[2]-1 ? rank+pp[0]*pp[1] :-1);

}
