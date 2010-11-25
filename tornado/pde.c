//----------------------- Calculation of PDE right part  ----------//
#define LEVEL extern
#include "head.h"
#define Re 1

double norma(double a,double b,double c,int order)
{
if(order==2) return ((a)*(a)+(b)*(b)+(c)*(c));
   else  return  pow(((a)*(a)+(b)*(b)+(c)*(c)),order/2.);
}

void pde(double t, double ****f, double ****df)
{
   int i,j,k,l,m;
   double dv1[4][3],dv2[4][3],dp1[3],dn1[3],w,dw;

   boundary_conditions(f,nut);

   for(i=0;i<m1;i++)
   for(j=ghost;j<mm2;j++)
   for(k=0;k<m3;k++)
     if(isType(node[i][k],NodeFluid) && !isType(node[i][k],NodeClued))
     {
      for(m=0;m<3;m++) {
         dp1[m]=dr(f[0],i,j,k,m+1,0,dx[m],ghost, approx);
         dn1[m]=dr(nut,i,j,k,m+1,0,dx[m],ghost, approx);
      }
      for(l=1;l<=3;l++) {
       for(m=0;m<3;m++) {
	 dv1[l][m]=dr(f[l],i,j,k,m+1,0,dx[m],ghost, approx);
	 dv2[l][m]=dr(f[l],i,j,k,m+1,1,dx[m]*dx[m],ghost, approx);
	 }
      }

      df[1][i][j][k]=nut[i][j][k]*(dv2[1][0]+dv2[1][1]+dv2[1][2]+r_1[i]*dv1[1][0]
				   -f[1][i][j][k]*r_2[i]-2*dv1[2][1]*r_1[i])
		     -dp1[0]+Ta*Ta/r_1[i]       +2*Ta*f[2][i][j][k]                   //forces of inertion
//                     +(j+n[2]<=ghost+2 ? coordin(k,2)*f[2][i][j][k] : 0)                //helical force
		     + (dn1[0]-f[1][i][j][k])*dv1[1][0]
		     + (dn1[1]-f[2][i][j][k])*dv1[1][1]
		     + (dn1[2]-f[3][i][j][k])*dv1[1][2]
		     - (dn1[1]-f[2][i][j][k])*r_1[i]*f[2][i][j][k]
		     ;
      df[2][i][j][k]=nut[i][j][k]*(dv2[2][0]+dv2[2][1]+dv2[2][2]+r_1[i]*dv1[2][0]
				   -f[2][i][j][k]*r_2[i]+2*dv1[1][1]*r_1[i])
		     -dp1[1]/*-dw/r_1[i]*/ -2*Ta*f[1][i][j][k]                   //forces of inertion
//                     -(j+n[2]<=ghost+2 ? (coordin(i,0)-rc)*f[2][i][j][k] :0)            //helical force
		     + (dn1[0]-f[1][i][j][k])*dv1[2][0]
		     + (dn1[1]-f[2][i][j][k])*dv1[2][1]
		     + (dn1[2]-f[3][i][j][k])*dv1[2][2]
		     + (dn1[1]-f[1][i][j][k])*r_1[i]*f[2][i][j][k]
		     ;
      df[3][i][j][k]=nut[i][j][k]*(dv2[3][0]+dv2[3][1]+dv2[3][2]+r_1[i]*dv1[3][0])
		     +(i==ghost? Gr : 0)
		     -dp1[2]
		     + (dn1[0]-f[1][i][j][k])*dv1[3][0]
		     + (dn1[1]-f[2][i][j][k])*dv1[3][1]
		     + (dn1[2]-f[3][i][j][k])*dv1[3][2];
      df[0][i][j][k]= -(dv1[1][0]+dv1[2][1]+dv1[3][2]+f[1][i][j][k]*r_1[i])/Gamma;

  } //global for
   return;
}

double deviation(double ****f,int i,int j,int k)
{
const int size_okr=min(1,ghost);
double flux = 0;
int kol=0,l;
   for(l=1;l<=size_okr;l++)
            {
            if(i-l>0)
	       {flux += norma(f[1][i][j][k]-f[1][i-l][j][k],f[2][i][j][k]-f[2][i-l][j][k],f[3][i][j][k]-f[3][i-l][j][k],2);
            	kol++;
               }
            if(i+l<m1)
               {flux += norma(f[1][i][j][k]-f[1][i+l][j][k],f[2][i][j][k]-f[2][i+l][j][k],f[3][i][j][k]-f[3][i+l][j][k],2);
		kol++;
               }
            if(j-l>0)
               {flux += norma(f[1][i][j][k]-f[1][i][j-l][k],f[2][i][j][k]-f[2][i][j-l][k],f[3][i][j][k]-f[3][i][j-l][k],2);
            	kol++;
               }
            if(j+l<m2)
               {flux += norma(f[1][i][j][k]-f[1][i][j+l][k],f[2][i][j][k]-f[2][i][j+l][k],f[3][i][j][k]-f[3][i][j+l][k],2);
            	kol++;
               }
	    if(k-l>0)
               {flux += norma(f[1][i][j][k]-f[1][i][j][k-l],f[2][i][j][k]-f[2][i][j][k-l],f[3][i][j][k]-f[3][i][j][k-l],2);
                kol++;
               }
            if(k+l<m3)
               {flux += norma(f[1][i][j][k]-f[1][i][j][k+l],f[2][i][j][k]-f[2][i][j][k+l],f[3][i][j][k]-f[3][i][j][k+l],2);
		kol++;
               }
           };
   flux /= kol;
return(flux);
}

void nut_by_flux(double ****f, double ***nut, double dt) //calculating nu_turbulent by velocity fluctuations
{
int i,j,k,l;
double koef,r1,z1,rho,tmp;
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
for(k=0;k<m3;k++)
   {
/*   double tmp = maschtab*pow(
           (nl[2]*s_func[k][0] + nl[1]*s_func[k][1] + nl[0]*s_func[k][2])*pow(dx[2],4),
                         1./3);*/
/*   for(i=0;i<m1;i++)
     {
     r1 = coordin(i,0);   z1 = coordin(k,2);
     rho = 1 - sqrt(r1*r1 + z1*z1)/R;
     if(rho<0) tmp = 0;
	  else tmp = sqrt(2*Re)*(0.32*rho*(1-rho)+0.013*(2*rho-1)*(2*rho-1)*rho);   // Шлихтинг
	  else tmp = sqrt(2*Re)*0.276467*rho*(2.08 - 2.8*rho + rho*rho)*(0.64 - 1.2*rho + rho*rho);   // Рейнольдс
       for(j=ghost;j<mm2;j++)
	if(isType(node[i][k],NodeFluid) && !isType(node[i][k],NodeClued))
	  nut[i][j][k] = (1. + maschtab*tmp)/Re;
     }*/

   }
}

void  boundary_conditions(double ****f, double ***nut)
{
   int i, j, k, l, req_numS=0, req_numR=0;
   int r1,r2,z1,z2;
   int /*flag,cnt,*/z[4],tag=10;
   double vrho,vphi,vth;
   z[1]=z[2]=z[3]=-1; z[0]=1;  //  влияет на вид гран.условий (-1:жесткие, 1:свободные)

  /*============================ divertor =================================*/
  if(t_cur<00)                  // divertors are off till t sec
  if(n[2]==0)
  for(i=0;i<m1;i++)
    for(j=ghost;j<=ghost;j++)
      for(k=0;k<m3;k++)
         {
         vrho = f[3][i][j][k]*costh[i][k]+f[1][i][j][k]*sinth[i][k];
         vth  = -f[3][i][j][k]*sinth[i][k]+f[1][i][j][k]*costh[i][k];
	 vphi = sqrt(pow(f[2][i][j][k],2)+vth*vth);           //sqrt(vfi*vfi+vth*vth)
         f[1][i][j][k] = vrho*sinth[i][k]+vphi*costh[i][k]*sin(chi[i][k]);
	 f[2][i][j][k] = vphi*cos(chi[i][k]);
	 f[3][i][j][k] = vrho*costh[i][k]-vphi*sinth[i][k]*sin(chi[i][k]);
	 }

  /*-------------------------------- exchanging of ghosts -------------------------------------*/
// exchanging in phi-direction - periodical directions first
 if(pr_neighbour[2]>-1)
  if(pr_neighbour[2]==rank) CopyGridToBuffer(f,nut,buf_recv[2],0,n2,0,m1-1,mm2-1,m3-1);
         else { CopyGridToBuffer(f,nut,buf_send[2],0,ghost,0,m1-1,2*ghost-1,m3-1);
		MPI_Isend(buf_send[2],buf_size[1],MPI_DOUBLE,pr_neighbour[2],tag+2,MPI_COMM_WORLD,&SendRequest[req_numS++]);
		MPI_Irecv(buf_recv[2],buf_size[1],MPI_DOUBLE,pr_neighbour[2],tag+3,MPI_COMM_WORLD,&RecvRequest[req_numR++]);
		}
 if(pr_neighbour[3]>-1)
  if(pr_neighbour[3]==rank) CopyGridToBuffer(f,nut,buf_recv[3],0,ghost,0,m1-1,2*ghost-1,m3-1);
	 else { CopyGridToBuffer(f,nut,buf_send[3],0,n2,0,m1-1,mm2-1,m3-1);
		MPI_Isend(buf_send[3],buf_size[1],MPI_DOUBLE,pr_neighbour[3],tag+3,MPI_COMM_WORLD,&SendRequest[req_numS++]);
		MPI_Irecv(buf_recv[3],buf_size[1],MPI_DOUBLE,pr_neighbour[3],tag+2,MPI_COMM_WORLD,&RecvRequest[req_numR++]);
		}

    MPI_Waitall(req_numR,RecvRequest,statuses);
   if(statuses[0].MPI_ERROR) putlog("bc:error during transfer=",numlog++);

  if(pr_neighbour[2]>-1) CopyBufferToGrid(f,nut,buf_recv[2],0,0,0,m1-1,ghost-1,m3-1);
  if(pr_neighbour[3]>-1) CopyBufferToGrid(f,nut,buf_recv[3],0,mm2,0,m1-1,m2-1,m3-1);

// exchanging in r-direction
 if(pr_neighbour[0]>-1)
  if(pr_neighbour[0]==rank) CopyGridToBuffer(f,nut,buf_recv[0],n1,0,0,mm1-1,m2-1,m3-1);
         else { CopyGridToBuffer(f,nut,buf_send[0],ghost,0,0,2*ghost-1,m2-1,m3-1);
		MPI_Isend(buf_send[0],buf_size[0],MPI_DOUBLE,pr_neighbour[0],tag,MPI_COMM_WORLD,&SendRequest[req_numS++]);
		MPI_Irecv(buf_recv[0],buf_size[0],MPI_DOUBLE,pr_neighbour[0],tag+1,MPI_COMM_WORLD,&RecvRequest[req_numR++]);
		}
 if(pr_neighbour[1]>-1)
  if(pr_neighbour[1]==rank) CopyGridToBuffer(f,nut,buf_recv[1],ghost,0,0,2*ghost-1,m2-1,m3-1);
         else { CopyGridToBuffer(f,nut,buf_send[1],n1,0,0,mm1-1,m2-1,m3-1);
		MPI_Isend(buf_send[1],buf_size[0],MPI_DOUBLE,pr_neighbour[1],tag+1,MPI_COMM_WORLD,&SendRequest[req_numS++]);
		MPI_Irecv(buf_recv[1],buf_size[0],MPI_DOUBLE,pr_neighbour[1],tag,MPI_COMM_WORLD,&RecvRequest[req_numR++]);
		}

    MPI_Waitall(req_numR,RecvRequest,statuses);
  if(statuses[0].MPI_ERROR) putlog("bc:error during transfer=",numlog++);
  if(pr_neighbour[0]>-1) CopyBufferToGrid(f,nut,buf_recv[0],0,0,0,ghost-1,m2-1,m3-1);
  if(n[0]==0)      // at axis
	for(i=0;i<ghost;i++)
	for(j=0;j<m2;j++)
	for(k=0;k<m3;k++)
		f[2][i][j][k] *= -1;      // continuity of velocity at axis
  if(pr_neighbour[1]>-1) CopyBufferToGrid(f,nut,buf_recv[1],mm1,0,0,m1-1,m2-1,m3-1);

// exchanging in z-direction
 if(pr_neighbour[4]>-1)
  if(pr_neighbour[4]==rank) CopyGridToBuffer(f,nut,buf_recv[4],0,0,mm3,m1-1,m2-1,m3-1);
	 else { CopyGridToBuffer(f,nut,buf_send[4],0,0,ghost,m1-1,m2-1,2*ghost-1);
		MPI_Isend(buf_send[4],buf_size[2],MPI_DOUBLE,pr_neighbour[4],tag+4,MPI_COMM_WORLD,&SendRequest[req_numS++]);
		MPI_Irecv(buf_recv[4],buf_size[2],MPI_DOUBLE,pr_neighbour[4],tag+5,MPI_COMM_WORLD,&RecvRequest[req_numR++]);
		}
 if(pr_neighbour[5]>-1)
  if(pr_neighbour[5]==rank) CopyGridToBuffer(f,nut,buf_recv[5],0,0,0,m1-1,m2-1,ghost-1);
        else { CopyGridToBuffer(f,nut,buf_send[5],0,0,n3,m1-1,m2-1,mm3-1);
	       MPI_Isend(buf_send[5],buf_size[2],MPI_DOUBLE,pr_neighbour[5],tag+5,MPI_COMM_WORLD,&SendRequest[req_numS++]);
	       MPI_Irecv(buf_recv[5],buf_size[2],MPI_DOUBLE,pr_neighbour[5],tag+4,MPI_COMM_WORLD,&RecvRequest[req_numR++]);
	       }

    MPI_Waitall(req_numR,RecvRequest,statuses);
  if(statuses[0].MPI_ERROR) putlog("bc:error during transfer=",numlog++);
										  

  if(pr_neighbour[4]>-1) CopyBufferToGrid(f,nut,buf_recv[4],0,0,0,m1-1,m2-1,ghost-1);
  if(pr_neighbour[5]>-1) CopyBufferToGrid(f,nut,buf_recv[5],0,0,mm3,m1-1,m2-1,m3-1);

//    MPI_Barrier(MPI_COMM_WORLD);
//    MPI_Startall(req_numR,RecvRequest);
//    MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,statuses);
//    if(flag==0) putlog("bc:error during transfer=",numlog++);
//    MPI_Testall(req_numR,RecvRequest,&flag,statuses);
//    MPI_Get_count(statuses,MPI_DOUBLE,&cnt);
//    MPI_Waitall(req_numR,RecvRequest,statuses);

/*----------------------- filling of ghost nodes ------------------------------*/

  for(i=0;i<m1;i++)
     for(j=0;j<m2;j++)
	for(k=0;k<m3;k++)
	  if(isType(node[i][k],NodeGhostFluid) && !isType(node[i][k],NodeClued)) {
		r2=(r1=floor(refr[i][k]))+1;
		z2=(z1=floor(refz[i][k]))+1;
		for(l=0;l<nvar;l++)
		   f[l][i][j][k] = ( (refr[i][k]-r2)*(f[l][r1][j][z1]*(refz[i][k]-z2)-f[l][r1][j][z2]*(refz[i][k]-z1))
				   + (refr[i][k]-r1)*(f[l][r2][j][z2]*(refz[i][k]-z1)-f[l][r2][j][z1]*(refz[i][k]-z2))
				    ) * z[l];
//                   nut[i][j][k] = (refr[i][k]-r2)*(nut[r1][j][z1]*(refz[i][k]-z2)-nut[r1][j][z2]*(refz[i][k]-z1))
//                                + (refr[i][k]-r1)*(nut[r2][j][z2]*(refz[i][k]-z1)-nut[r2][j][z1]*(refz[i][k]-z2));
		   nut[i][j][k] = 1./Re;
		}

  return;
}

void  init_conditions()
{
   int i,j,k,l;
   double r;
//   double k1,k2,k3;

// -------- filling of nodes' types +reflections rel circle + nut ---------------
   for(i=0;i<m1;i++)
   for(k=0;k<m3;k++)
       {
       node[i][k] = NodeUnknown;
       if(pr_neighbour[0]!=-1 && i<ghost ||
           pr_neighbour[1]!=-1 && i>=mm1  ||
           pr_neighbour[4]!=-1 && k<ghost ||
           pr_neighbour[5]!=-1 && k>=mm3) setType(&node[i][k],NodeClued);
       }
   for(i=0;i<m1;i++)
   for(k=0;k<m3;k++)
       {
     //regions
        refr[i][k] = i;
	refz[i][k] = k;
        if(i+n[0]<ghost)  {setType(&node[i][k],NodeGhostFluid);
			   refr[i][k] = 2*ghost-i-1;}
	 else if(i+n[0]>=N1+ghost)
			  {setType(&node[i][k],NodeGhostFluid);
			   refr[i][k] = 2*mm1-i-1;}
	if(k+n[2]<ghost)  {setType(&node[i][k],NodeGhostFluid);
			   refz[i][k] = 2*ghost-k-1;}
	 else if(k+n[2]>=N3+ghost)
			  {setType(&node[i][k],NodeGhostFluid);
                           refz[i][k] = 2*mm3-k-1;}
        if(!isType(node[i][k],NodeGhostFluid)) setType(&node[i][k],NodeFluid);
       }

   for(i=0;i<m1;i++) { r_1[i] = 1./coordin(i,0); r_2[i] = r_1[i]*r_1[i]; }

// --------------- initial conditions -----------------------------------------
//   k1=2*M_PI/lfi;  k3=M_PI/l3;

if(!goon) {
   for(i=0;i<m1;i++)
   for(j=0;j<m2;j++)
   for(k=0;k<m3;k++)
      if(isType(node[i][k],NodeFluid)) {
	f[0][i][j][k]=0;
	f[1][i][j][k]=NoiseNorm*cos(2*M_PI*coordin(j,1)/R)*sin(2*M_PI*coordin(k,2)/R)
		      + Noise*((double)rand()-RAND_MAX/2)/RAND_MAX*
		       (R*R-pow(coordin(i,0),2) - pow(coordin(k,2),2))/R/R;
	f[2][i][j][k]=(parabole+Noise*((double)rand()-RAND_MAX/2)/RAND_MAX)*
		       (R*R-pow(coordin(i,0),2) - pow(coordin(k,2),2))/R/R;
	f[3][i][j][k]=NoiseNorm*sin(2*M_PI*coordin(j,1)/R)*sin(2*M_PI*coordin(k,2)/R)
		      + Noise*((double)rand()-RAND_MAX/2)/RAND_MAX*
		       (R*R-pow(coordin(i,0),2) - pow(coordin(k,2),2))/R/R;
	}
   for(i=0;i<m1;i++)
   for(j=0;j<m2;j++)
   for(k=0;k<m3;k++)
	nut[i][j][k]=1./Re;
//   struct_func(f,2,2,3);
   nmessage("Arrays were filled with initial values - calculation from beginning",-1,-1);
   } else nmessage("Arrays were filled with initial values - calculation is continuing",t_cur,count);
}

void init_parallel()
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
         if(mintime<0 || vtime<mintime) { mintime=vtime; nd1=i; nd2=j; }
         }
if(!goon) {                       //reading sizes from file when continuing
  pp[0]=divisors[nd1]; pp[1]=divisors[nd2]; pp[2]=size/pp[0]/pp[1];                 // number of procs along axes
          }
  pr[0] = rank%pp[0]; pr[1] = (rank/pp[0])%pp[1]; pr[2] = (rank/pp[0]/pp[1])%pp[2];  // coordinates of current subregion

/* dimensions of subregion:         global indicies of subregion origin          if there's nonequal subregions*/
  n1 = floor((double)N1/pp[0]);     n[0] = n1*pr[0] + min(pr[0],N1-pp[0]*n1);    if(pr[0]<N1-pp[0]*n1) n1++;              // dimensions of subregion
  n2 = floor((double)N2/pp[1]);     n[1] = n2*pr[1] + min(pr[1],N2-pp[1]*n2);    if(pr[1]<N2-pp[1]*n2) n2++;
  n3 = floor((double)N3/pp[2]);     n[2] = n3*pr[2] + min(pr[2],N3-pp[2]*n3);    if(pr[2]<N3-pp[2]*n3) n3++;
   if(n1<ghost || n2<ghost || n3<ghost) nrerror("Too small mesh or incorrect number of processes",0,0);

   m1 = n1+2*ghost;
   m2 = n2+2*ghost;
   m3 = n3+2*ghost;
   mm1 = ghost+n1;
   mm2 = ghost+n2;
   mm3 = ghost+n3;
                                                                                      // neighbours of subregion
  pr_neighbour[0] = (pr[0]>0       ? rank-1  :
  			(pr[1]<pp[1]/2 ? rank+pp[0]*pp[1]/2 : rank-pp[0]*pp[1]/2 )    // at axis
                    );
  pr_neighbour[1] = (pr[0]<pp[0]-1 ? rank+1           :-1);
  pr_neighbour[2] = (pr[1]>0       ? rank-pp[0]       :-1);
  if(pr[1]==0)       pr_neighbour[2] = rank+pp[0]*(pp[1]-1);     //periodic conditions
  pr_neighbour[3] = (pr[1]<pp[1]-1 ? rank+pp[0]       :-1);
  if(pr[1]==pp[1]-1) pr_neighbour[3] = rank-pp[0]*(pp[1]-1);     //periodic conditions
  pr_neighbour[4] = (pr[2]>0       ? rank-pp[0]*pp[1] :-1);
  pr_neighbour[5] = (pr[2]<pp[2]-1 ? rank+pp[0]*pp[1] :-1);

 buf_size[0]=m2*m3*(nvar+1)*ghost;
 buf_size[1]=m1*m3*(nvar+1)*ghost;
 buf_size[2]=m1*m2*(nvar+1)*ghost;

 for(i=0;i<3;i++)           //3-D
  for(j=0;j<=1;j++) {
   buf_send[j+2*i] = alloc_mem_1f(buf_size[i]);
   buf_recv[j+2*i] = alloc_mem_1f(buf_size[i]);
  }

   iop=fopen(NameStatFile,"w");
   fprintf(iop,"%d\n",rank);
   fprintf(iop,"%d\t%d\t%d\n",pr[0],pr[1],pr[2]);
   fprintf(iop,"%d\t%d\t%d\n",pp[0],pp[1],pp[2]);
   fprintf(iop,"%d\t%d\t%d\n",n[0],n[1],n[2]);
   fprintf(iop,"%d\t%d\t%d\n",n1,n2,n3);
   for(i=0;i<6;i++)
     fprintf(iop,"%d ",pr_neighbour[i]);
   fprintf(iop,"\n");
   fileclose(iop);

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

switch (sm*dir) {
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
	}
return(tmp/dx);
}

double coordin(int i, int dir)
                      //0-r,1-phi,2-z
{
 switch (dir)
 {
   case 0:  return dx[dir]*(i-ghost+0.5+n[dir]);
   case 1:  return dx[dir]*(i-ghost+0.5+n[dir]);
   case 2:  return dx[dir]*(i-ghost+0.5+n[dir]);
 }
    return(0);
}
