//----------------------- Calculation of PDE right part  ----------//
#define LEVEL extern
#include "head.h"
#define eta0 1
#define EPS 1e-10

double norma(double a,double b,double c,int order)
{
if(order==2) return ((a)*(a)+(b)*(b)+(c)*(c));
   else  return  pow(((a)*(a)+(b)*(b)+(c)*(c)),order/2.);
}

void omega(double t, double *w, double *dw)
{
double a=1,
       omega0=1,
       t_begin=0.1,
       t_end=omega0/a+t_begin;
if(t<t_end && t>t_begin) *dw=-a;
                    else *dw=0;
if(t<t_begin) *w=omega0;
 else if(t<t_end) *w=omega0-a*(t-t_begin);
 else *w=0;
}

void pde(double t, double ****f, double ****df)
{
   int i,j,k,l,m;
   double M;
   double dv1[7][3],dv2[7][3],dA11[7][3][3],dp1[3],dn1[3],w,dw;
   double mu[3];

   boundary_conditions(f,nut);

   for(i=0;i<3;i++) mu[i] = 2*dx[i]*dx[i];
//   omega(t,&w,&dw);
   for(i=0;i<m1;i++)
   for(j=ghost;j<mm2;j++)
   for(k=0;k<m3;k++)
     if(isType(node[i][k],NodeFluid) && !isType(node[i][k],NodeClued))
     {
//      if(j<=ghost && n[1]==0) continue;
      M = 1./(rc*coordin(i,0)+1);
      memset(dp1,0,3*sizeof(double));      memset(dn1,0,3*sizeof(double));
      memset(dv1,0,21*sizeof(double));     memset(dv2,0,21*sizeof(double));
      memset(dA11,0,63*sizeof(double));
      for(m=0;m<3;m++) {
         dp1[m]=dr(f[0],i,j,k,m+1,0,dx[m],(approx-1)/2, approx);
         dn1[m]=dr(nut,i,j,k,m+1,0,dx[m], (approx-1)/2, approx);
      }
      dp1[1] *= M; dn1[1] *= M;
      for(l=1;l<=3;l++) {
       for(m=0;m<3;m++) {
         dv1[l][m]=dr(f[l],i,j,k,m+1,0,dx[m],(approx-1)/2, approx);
         dv2[l][m]=dr(f[l],i,j,k,m+1,1,dx[m]*dx[m],(approx-1)/2, approx);
	 }
         dv1[l][1] *= M;     dv2[l][1]  *= M*M;
        dA11[l][0][1] = dA11[l][1][0] = d2cross(f[l],i,j,k,2,1,(approx-1)/2,approx);
        dA11[l][0][2] = dA11[l][2][0] = d2cross(f[l],i,j,k,3,1,(approx-1)/2,approx);
        dA11[l][1][2] = dA11[l][2][1] = d2cross(f[l],i,j,k,3,2,(approx-1)/2,approx);
       }

      df[1][i][j][k]=nut[i][j][k]*(dv2[1][0]+dv2[1][1]+dv2[1][2]+r_1[i]*dv1[1][0]
				   -f[1][i][j][k]*r_2[i]-2*dv1[2][1]*r_1[i])
		     -dp1[0]/Gamma/f[0][i][j][k]//+w*w/r_1[i]       +2*w*f[2][i][j][k]                   //forces of inertion
                -(dv1[1][0]<0 ? mu[0]*dv1[1][0]*(dp1[0]*dv1[1][0]+2*f[0][i][j][k]*dv2[1][0]) : 0)
              +(j+n[1]<=ghost+2 && t_cur>-1 ? coordin(k,2)*f[2][i][j][k] : 0)                //helical force
		     + (dn1[0]-f[1][i][j][k])*dv1[1][0]
		     + (dn1[1]-f[2][i][j][k])*dv1[1][1]
		     + (dn1[2]-f[3][i][j][k])*dv1[1][2]
		+nut[i][j][k]/3*(dv2[1][0]+dv1[1][0]*r_1[i]-f[1][i][j][k]*r_2[i]-dv1[2][1]*r_1[i]+dA11[2][0][1]+dA11[3][0][2])
/*                     - dvv(f[1],f[1],1,i,j,k,(approx-1)/2,approx)
                     - dvv(f[2],f[1],2,i,j,k,(approx-1)/2,approx)
                     - dvv(f[3],f[1],3,i,j,k,(approx-1)/2,approx)*/
                     - f[1][i][j][k]*f[1][i][j][k]*r_1[i]
		     - (dn1[1]-f[2][i][j][k])*r_1[i]*f[2][i][j][k]
		     ;
      df[2][i][j][k]=nut[i][j][k]*(dv2[2][0]+dv2[2][1]+dv2[2][2]+r_1[i]*dv1[2][0]
                                   -f[2][i][j][k]*r_2[i]+2*dv1[1][1]*r_1[i])
			 + (t_cur>-1 ? p1*pow(rc*coordin(i,0)+1,-1):0) - dp1[1]/Gamma/f[0][i][j][k]
                -(dv1[2][1]<0 ? mu[1]*dv1[2][1]*(dp1[1]*dv1[2][1]+2*f[0][i][j][k]*dv2[2][1]) : 0)
					//-dw/r_1[i] -2*w*f[1][i][j][k]                   //forces of inertion
		     + (dn1[0]-f[1][i][j][k])*dv1[2][0]
		     + (dn1[1]-f[2][i][j][k])*dv1[2][1]
		     + (dn1[2]-f[3][i][j][k])*dv1[2][2]
		+nut[i][j][k]/3*(dv1[1][1]*r_1[i]+dA11[1][0][1]+dv2[2][1]+dA11[3][1][2])
/*                     - dvv(f[1],f[2],1,i,j,k,(approx-1)/2,approx)
                     - dvv(f[2],f[2],2,i,j,k,(approx-1)/2,approx)
                     - dvv(f[3],f[2],3,i,j,k,(approx-1)/2,approx)*/
                     - f[1][i][j][k]*f[2][i][j][k]*r_1[i]
		     + (dn1[1]-f[1][i][j][k])*r_1[i]*f[2][i][j][k]
		     ;
      df[3][i][j][k]=nut[i][j][k]*(dv2[3][0]+dv2[3][1]+dv2[3][2]+r_1[i]*dv1[3][0])
		     -dp1[2]/Gamma/f[0][i][j][k]
                -(dv1[3][2]<0 ? mu[2]*dv1[3][2]*(dp1[2]*dv1[3][2]+2*f[0][i][j][k]*dv2[3][2]) : 0)
              -(j+n[1]<=ghost+2 && t_cur>-1 ? (coordin(i,0)-rc)*f[2][i][j][k] :0)            //helical force
		     + (dn1[0]-f[1][i][j][k])*dv1[3][0]
		     + (dn1[1]-f[2][i][j][k])*dv1[3][1]
		     + (dn1[2]-f[3][i][j][k])*dv1[3][2]
      		+nut[i][j][k]/3*(dv1[1][2]*r_1[i]+dA11[1][0][2]+dA11[2][1][2]+dv2[3][2])
/*                     - dvv(f[1],f[2],1,i,j,k,(approx-1)/2,approx)
                     - dvv(f[2],f[2],2,i,j,k,(approx-1)/2,approx)
                     - dvv(f[3],f[2],3,i,j,k,(approx-1)/2,approx)*/
                     - f[1][i][j][k]*f[3][i][j][k]*r_1[i]
                     ;
//      df[0][i][j][k]= -f[0][i][j][k]*(dv1[1][0]+dv1[2][1]+dv1[3][2]+f[1][i][j][k]*r_1[i])/Gamma;
      df[0][i][j][k]= -(dvv(f[0],f[1],1,i,j,k,(approx-1)/2,approx)
      		       +dvv(f[0],f[2],2,i,j,k,(approx-1)/2,approx)
                       +dvv(f[0],f[3],3,i,j,k,(approx-1)/2,approx)
                       +f[0][i][j][k]*f[1][i][j][k]*r_1[i]);
/*      df[1][i][j][k] += Gamma*df[0][i][j][k]*f[1][i][j][k];
      df[2][i][j][k] += Gamma*df[0][i][j][k]*f[2][i][j][k];
      df[3][i][j][k] += Gamma*df[0][i][j][k]*f[3][i][j][k];*/
//      df[0][i][j][k] = df[1][i][j][k] = df[2][i][j][k] = df[3][i][j][k] = 0;
      }   else df[0][i][j][k] = df[1][i][j][k] = df[2][i][j][k] = df[3][i][j][k] = 0;
   for(i=0;i<m1;i++)
   for(k=0;k<m3;k++)
      if(isType(node[i][k],NodeGhostFluid))
      for(j=0;j<m2;j++)
        for(l=1;l<=3;l++)
                f[l][i][j][k] = 0;
   for(i=0;i<m1;i++)
   for(j=ghost;j<mm2;j++)
   for(k=0;k<m3;k++)
     if(isType(node[i][k],NodeMagn) && !isType(node[i][k],NodeClued))
      {
      for(l=4;l<=6;l++) {
       for(m=0;m<3;m++) {
         dv1[l][m]=dr(f[l],i,j,k,m+1,0,dx[m],(approx-1)/2, approx);
         dv2[l][m]=dr(f[l],i,j,k,m+1,1,dx[m]*dx[m],(approx-1)/2, approx);
	 }
        dA11[l][0][1] = dA11[l][1][0] = d2cross(f[l],i,j,k,2,1,(approx-1)/2,approx);
        dA11[l][0][2] = dA11[l][2][0] = d2cross(f[l],i,j,k,3,1,(approx-1)/2,approx);
        dA11[l][1][2] = dA11[l][2][1] = d2cross(f[l],i,j,k,3,2,(approx-1)/2,approx);
        }

       df[4][i][j][k]=Rm*(f[2][i][j][k]*(dv1[5][0]-dv1[4][1]+f[5][i][j][k]*r_1[i])-f[3][i][j][k]*(dv1[4][2]-dv1[6][0]))
//                      +Rm*f[2][i][j][k]	 // for induced field
                     +eta[i][j][k]*(dv2[4][0]+dv2[4][1]+dv2[4][2]+r_1[i]*dv1[4][0]-f[4][i][j][k]*r_2[i]-2*dv1[5][1]*r_1[i])
                     +(eta0-eta[i][j][k])*(dv2[4][0]+dv1[4][0]*r_1[i]-f[4][i][j][k]*r_2[i]-dv1[5][1]*r_1[i]+dA11[5][0][1]+dA11[6][0][2]);
       df[5][i][j][k]=Rm*(f[3][i][j][k]*(dv1[6][1]-dv1[5][2])-f[1][i][j][k]*(dv1[5][0]-dv1[4][1]+f[5][i][j][k]*r_1[i]))
//                     -Rm*f[1][i][j][k]	// for induced field
                     +eta[i][j][k]*(dv2[5][0]+dv2[5][1]+dv2[5][2]+r_1[i]*dv1[5][0]-f[5][i][j][k]*r_2[i]+2*dv1[4][1]*r_1[i])
                     +(eta0-eta[i][j][k])*(dv1[4][1]*r_1[i]+dA11[4][0][1]+dv2[5][1]+dA11[6][1][2]);
       df[6][i][j][k]=Rm*(f[1][i][j][k]*(dv1[4][2]-dv1[6][0])-f[2][i][j][k]*(dv1[6][1]-dv1[5][2]))
                     +eta[i][j][k]*(dv2[6][0]+dv2[6][1]+dv2[6][2]+r_1[i]*dv1[6][0])
                     +(eta0-eta[i][j][k])*(dv1[4][2]*r_1[i]+dA11[4][0][2]+dA11[5][1][2]+dv2[6][2]);
/*       df[4][i][j][k]=eta[i][j][k]*(dv2[4][0]+dv2[4][1]+dv2[4][2]+r_1[i]*dv1[4][0]-f[4][i][j][k]*r_2[i]-2*dv1[5][1]*r_1[i]);
       df[5][i][j][k]=f[3][i][j][k]*(dv1[6][1]-dv1[5][2])-f[1][i][j][k]*(dv1[5][0]-dv1[4][1]+f[5][i][j][k]*r_1[i])
                     +eta[i][j][k]*(dv2[5][0]+dv2[5][1]+dv2[5][2]+r_1[i]*dv1[5][0]-f[5][i][j][k]*r_2[i]+2*dv1[4][1]*r_1[i])
                     +(1./Rm-eta[i][j][k])*(dv1[4][1]*r_1[i]+dA11[4][0][1]+dv2[5][1]+dA11[6][1][2]);
       df[6][i][j][k]=0;*/
      }   else df[4][i][j][k] = df[5][i][j][k] = df[6][i][j][k] = 0;
   return;
}

void fill_velocity(double t, double ****f)
{
double r1, phi1, z1, rho, vrho, vth, vphi;
int pp[3];
int error=0, tmpr;
int i,j,k,l;
float tmpf; double tmpd; long tmpi;  char tmpc;
FILE *inp;
/*       for(i=0;i<m1;i++)
         for(j=0;j<m2;j++)
            for(k=0;k<m3;k++)
            if(isType(node[i][k],NodeFluid))
            {
            r1 = coordin(i,0);  phi1 = coordin(j,1);  z1 = coordin(k,2);
            rho=sqrt(r1*r1 + z1*z1);
            vrho = 0;
            vth  = vtheta_given(0,rho,Rfl,phi1)/(1+r1*rc);
/*            vphi = vfi_given(t_cur*Tunit,rho,Rfl);      // nonstationary time=t_cur*Tunit, max magnitude=0.0990906
            vth = 0;                        *
            vphi = vfi_given(0,rho,Rfl);
//            vth = rho*vphi/Rfl;
            f[1][i][j][k] = vrho*sinth[i][k]+vth*costh[i][k];
            f[2][i][j][k] = vphi;
            f[3][i][j][k] = vrho*costh[i][k]-vth*sinth[i][k];
            }    else f[1][i][j][k] = f[2][i][j][k] = f[3][i][j][k] = 0;
*/
 inp = fileopen("velocity.snp",-1);

 read_tilleq(inp,'=','y');   if(fscanf(inp,"%lf",&tmpd)==0) error=1; //	current time
 read_tilleq(inp,'=','y');   if(fscanf(inp,"%ld",&tmpi)==0) error=1; // current iteration
 read_tilleq(inp,'=','y');   if(fscanf(inp,"%c%d%c%d%c%d%c",&tmpc,&pp[0],&tmpc,&pp[1],&tmpc,&pp[2],&tmpc)<7) error=1; // number of processors along axes
 if(pp[0]*pp[1]*pp[2]!=1) nrerror("Can't read velocity data from many processors.",-1,-1);
 read_tilleq(inp,'=','y');   if(fscanf(inp,"%ld",&tmpi)==0) error=1; // Number of points along x
 read_tilleq(inp,'=','y');   if(fscanf(inp,"%ld",&tmpi)==0) error=1; // Number of points along y
 read_tilleq(inp,'=','y');   if(fscanf(inp,"%ld",&tmpi)==0) error=1; // Number of points along z
 read_tilleq(inp,'=','y');   if(fscanf(inp,"%lf",&tmpd)==0) error=1; //Reynolds number 
 read_tilleq(inp,0x0A,'y');
    if(error) nrerror("Error in reading velocity",t_cur,tmpi);
 for(tmpr=0;tmpr<=0*rank;tmpr++)           //reading until arrays of this process (if rank=0, one processor snap)
 {
 for(l=0;l<=3;l++)                    // reading f
    for(i=0;i<N1+2*ghost;i++)
    for(j=0;j<N2+2*ghost;j++)
    for(k=0;k<N3+2*ghost;k++)
       {
       if(fread(&tmpf,sizeof(float),1,inp)<1) nrerror("Error in reading velocity",0.,(k+100*(j+100*(i+l*100))));
	   if(l==0) continue;
       if(i>=n[0] && i<n[0]+m1 &&
          j>=n[1] && j<n[1]+m2 &&
          k>=n[2] && k<n[2]+m3)
                {
                if(isType(node[i-n[0]][k-n[2]],NodeFluid)) f[l][i-n[0]][j-n[1]][k-n[2]] = tmpf;
                        else f[l][i-n[0]][j-n[1]][k-n[2]] = 0;
                }
       }
 }
 fileclose(inp);    
}

double deviation(double ****f,int i,int j,int k)
{
const int size_okr=min(max_okr,ghost);
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

void nut_by_flux(double ****f,double dt) //calculating nu_turbulent by velocity fluctuations
{
//int i,j,k,l;
//double koef;
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
//for(k=0;k<m3;k++)
   {
/*   double tmp = maschtab*pow(
	   (nl[2]*s_func[k][0] + nl[1]*s_func[k][1] + nl[0]*s_func[k][2])*pow(dx[2],4),
                         1./3);*/
/*   double tmp=0;
   for(i=0;i<m1;i++)
     {
     r1 = coordin(i,0);   z1 = coordin(k,2);
     rho = 1 - sqrt(r1*r1 + z1*z1)/R;
     if(rho<0) tmp = 0;
	  else tmp = sqrt(2*Re)*(0.32*rho*(1-rho)+0.013*(2*rho-1)*(2*rho-1)*rho);   // Шлихтинг
//          else tmp = sqrt(2*Re)*0.276467*rho*(2.08 - 2.8*rho + rho*rho)*(0.64 - 1.2*rho + rho*rho);   // Рейнольдс
       for(j=ghost;j<mm2;j++)
          if(isType(node[i][k],NodeFluid) && !isType(node[i][k],NodeClued))
	  nut[i][j][k] = (1. + maschtab*tmp)/Re;  
	 }
	  */
   }
}

void ghost_filling(double ****f, double ***nut)
{
   int i,j,k,l;
   int i1,i2,k1,k2;
   int r1,r2,z1,z2;
   double rfict_1,rin_1,rrel;
   double ftemp[7];
   int /*flag,cnt,*/z[7],znorm,ztau;
   double vn, An;
   z[1]=z[2]=z[3]=-1; z[0]=1;  //  влияет на вид гран.условий (-1:жесткие, 1:свободные)
   znorm = -1; ztau = -1;

/*----------------------- filling of ghost nodes ------------------------------*/

  for(i=0;i<m1;i++)
     for(j=0;j<m2;j++)
	for(k=0;k<m3;k++)
        {
          if(isType(node[i][k],NodeGhostFluid) && !isType(node[i][k],NodeClued)) {
                r2=(r1=floor(refr_f[i][k]))+1;
                z2=(z1=floor(refz_f[i][k]))+1;
        if(r1<0 || r2>m1-1) printf("refr doesn't fit: (%d,%d)",i,k);
        if(z1<0 || z2>m3-1) printf("refz doesn't fit: (%d,%d)",i,k);
        if(r1<0 || r2>m1-1) nrerror("refr doesn't fit",t_cur,1000*i+k);
        if(z1<0 || z2>m3-1) nrerror("refz doesn't fit",t_cur,1000*i+k);
                for(l=0;l<4;l++)
		 {
		 f[l][i][j][k] = 0;
                 ftemp[l] = (refr_f[i][k]-r2)*(f[l][r1][j][z1]*(refz_f[i][k]-z2)-f[l][r1][j][z2]*(refz_f[i][k]-z1))
                                   + (refr_f[i][k]-r1)*(f[l][r2][j][z2]*(refz_f[i][k]-z1)-f[l][r2][j][z1]*(refz_f[i][k]-z2))
                                    ;
/*                 if(i==r1 && k==z1) f[l][i][j][k] /= (1-(refr[i][k]-r2)*(refz[i][k]-z2)*z[l]);
                 if(i==r1 && k==z2) f[l][i][j][k] /= (1+(refr[i][k]-r2)*(refz[i][k]-z1)*z[l]);
                 if(i==r2 && k==z1) f[l][i][j][k] /= (1-(refr[i][k]-r1)*(refz[i][k]-z1)*z[l]);
                 if(i==r2 && k==z2) f[l][i][j][k] /= (1+(refr[i][k]-r1)*(refz[i][k]-z2)*z[l]);
     если есть внутренние фиктивные, то этот отрывок не должен влиять
     */
     		  }

//                   nut[i][j][k] = (refr[i][k]-r2)*(nut[r1][j][z1]*(refz[i][k]-z2)-nut[r1][j][z2]*(refz[i][k]-z1))
//                                + (refr[i][k]-r1)*(nut[r2][j][z2]*(refz[i][k]-z1)-nut[r2][j][z1]*(refz[i][k]-z2));
                vn = ( ftemp[1]*(refr_f[i][k]-i)+ftemp[3]*(refz_f[i][k]-k) )/
                     ( (refr_f[i][k]-i)*(refr_f[i][k]-i) + (refz_f[i][k]-k)*(refz_f[i][k]-k) );
/*                for(l=0;l<nvar;l++)
                   f[l][i][j][k] = z[l]*ftemp[l];*/
//                rfict_1 = r_1[i]; rrel = 2*R*rfict_1 - 1;
                f[0][i][j][k] = z[0]*ftemp[0];
                f[1][i][j][k] = ztau*ftemp[1] + (znorm-ztau)*vn*(refr_f[i][k]-i);
                f[2][i][j][k] = ztau*ftemp[2]
//                             * (rrel+1-rfict_1*(i-i1)*dx[0]) / (rrel+1+rfict_1*(i-i1)*dx[0])
                              ;
                f[3][i][j][k] = ztau*ftemp[3] + (znorm-ztau)*vn*(refz_f[i][k]-k);
                nut[i][j][k] = 1./Re;
                }
          if(isType(node[i][k],NodeGhostMagn)) {
                i1=floor(refr_m[i][k]+0.5);
                rin_1 = r_1[i1]; rfict_1 = r_1[i]; rrel = (rin_1==0)? 1 : rfict_1/rin_1;
                k1=floor(refz_m[i][k]+0.5);
                An = ( f[4][i1][j][k1]*(i1-i)+f[6][i1][j][k1]*(k1-k) )/
                     ( (i1-i)*(i1-i) + (k1-k)*(k1-k) );
/*                for(l=4;l<=6;l++)
                   f[l][i][j][k] = z[l]*f[l][i1][j][k1];*/
                f[4][i][j][k] = ztau*f[4][i1][j][k1] + (znorm-ztau)*An*(i1-i);
//                f[5][i][j][k] = ztau*f[5][i1][j][k1];
                f[5][i][j][k] = ztau*f[5][i1][j][k1] *
                              (rrel+1-rfict_1*(i-i1)*dx[0]) / (rrel+1+rfict_1*(i-i1)*dx[0]) ;
/*                switch (i1-i)
                {
                  case -1: case -3: case -5:
                           f[5][i][j][k] = f[5][i1][j][k1] *
                              (rin+rfict-(i-i1)*dx[0]) / (rin+rfict+(i-i1)*dx[0]); //simple condition but doesn't work
                  break;
		  case 1:  f[5][i][j][k] = f[5][i1][j][k1] * rfict*(3*rin-rfict) / (rin*(3*rfict-rin));
                  break;
                  case 3: f[5][i][j][k] = f[5][i1][j][k1] * (1. + (rin*rin-rfict*rfict)/6./rfict/rin);
                           i2 = (2*i1+i)/3;
                           f[5][i][j][k]+= f[5][i2][j][k1] * 3.*(rin*rin-rfict*rfict)*(3*rfict+rin)/2./rfict/(5*rfict+rin)/(rfict+2*rin);
                  break;
                  case 5: f[5][i][j][k] = f[5][i1][j][k1] * (1. + (rin*rin-rfict*rfict)/10./rfict/rin);
                           i2 = (4*i1+i)/5;
                           f[5][i][j][k]+= f[5][i2][j][k1] * (rin*rin-rfict*rfict)/rfict/(rfict+4*rin);
                           i2 = (3*i1+2*i)/5;
                           f[5][i][j][k]+= f[5][i2][j][k1] * (rin*rin-rfict*rfict)*(13*rfict+7*rin)/2./rfict/(2*rfict+3*rin)/(7*rfict+3*rin);
                  break;
                }*/
                f[6][i][j][k] = ztau*f[6][i1][j][k1] + (znorm-ztau)*An*(k1-k);
                if((i1-i)*(k1-k)!=0) {
                                     f[4][i][j][k] = -f[4][i1][j][k1];
                                     f[6][i][j][k] = -f[6][i1][j][k1];
                                     }
                }
        }
}

void interprocessor_communication(double ****f, double***nut)
{
   int req_numS=0, req_numR=0,tag=10;
   int reslen;
   char msg_err[100];       //for putlog+mpi_error

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

   MPI_Waitall(req_numS,SendRequest,statuses);
   if(req_numS && statuses[0].MPI_ERROR) {putlog("bc:error during send=",numlog++);
                               MPI_Error_string(statuses[0].MPI_ERROR,msg_err,&reslen);
                               msg_err[reslen++] = ','; msg_err[reslen]= 0;
                               putlog(msg_err,numlog++);
                               }    else numlog++;
   req_numS = 0;
   MPI_Waitall(req_numR,RecvRequest,statuses);
   if(req_numR && statuses[0].MPI_ERROR) {putlog("bc:error during receive=",numlog++);
                               MPI_Error_string(statuses[0].MPI_ERROR,msg_err,&reslen);
                               msg_err[reslen++] = ','; msg_err[reslen]= 0;
                               putlog(msg_err,numlog++);
                               }    else numlog++;
   req_numR = 0;
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

   MPI_Waitall(req_numS,SendRequest,statuses);
   if(req_numS && statuses[0].MPI_ERROR) {putlog("bc:error during send=",numlog++);
                               MPI_Error_string(statuses[0].MPI_ERROR,msg_err,&reslen);
                               msg_err[reslen++] = ','; msg_err[reslen]= 0;
                               putlog(msg_err,numlog++);
                               }    else numlog++;
   req_numS = 0;
  MPI_Waitall(req_numR,RecvRequest,statuses);
   if(req_numR && statuses[0].MPI_ERROR) {putlog("bc:error during receive=",numlog++);
                               MPI_Error_string(statuses[0].MPI_ERROR,msg_err,&reslen);
                               msg_err[reslen++] = ','; msg_err[reslen]= 0;
                               putlog(msg_err,numlog++);
                               }    else numlog++;
   req_numR = 0;
  if(pr_neighbour[0]>-1) CopyBufferToGrid(f,nut,buf_recv[0],0,0,0,ghost-1,m2-1,m3-1);
  if(pr_neighbour[1]>-1) CopyBufferToGrid(f,nut,buf_recv[1],mm1,0,0,m1-1,m2-1,m3-1);

// exchanging in z-direction
 if(pr_neighbour[4]>-1)
  if(pr_neighbour[4]==rank) CopyGridToBuffer(f,nut,buf_recv[4],0,0,n3,m1-1,m2-1,mm3-1);
         else { CopyGridToBuffer(f,nut,buf_send[4],0,0,ghost,m1-1,m2-1,2*ghost-1);
                MPI_Isend(buf_send[4],buf_size[2],MPI_DOUBLE,pr_neighbour[4],tag+4,MPI_COMM_WORLD,&SendRequest[req_numS++]);
				MPI_Irecv(buf_recv[4],buf_size[2],MPI_DOUBLE,pr_neighbour[4],tag+5,MPI_COMM_WORLD,&RecvRequest[req_numR++]);
                }
 if(pr_neighbour[5]>-1)
  if(pr_neighbour[5]==rank) CopyGridToBuffer(f,nut,buf_recv[5],0,0,ghost,m1-1,m2-1,2*ghost-1);
        else { CopyGridToBuffer(f,nut,buf_send[5],0,0,n3,m1-1,m2-1,mm3-1);
               MPI_Isend(buf_send[5],buf_size[2],MPI_DOUBLE,pr_neighbour[5],tag+5,MPI_COMM_WORLD,&SendRequest[req_numS++]);
               MPI_Irecv(buf_recv[5],buf_size[2],MPI_DOUBLE,pr_neighbour[5],tag+4,MPI_COMM_WORLD,&RecvRequest[req_numR++]);
               }

   MPI_Waitall(req_numS,SendRequest,statuses);
   if(req_numS && statuses[0].MPI_ERROR) {putlog("bc:error during send=",numlog++);
                               MPI_Error_string(statuses[0].MPI_ERROR,msg_err,&reslen);
                               msg_err[reslen++] = ','; msg_err[reslen]= 0;
                               putlog(msg_err,numlog++);
                               }    else numlog++;
   req_numS = 0;
  MPI_Waitall(req_numR,RecvRequest,statuses);
   if(req_numR && statuses[0].MPI_ERROR) {putlog("bc:error during receive=",numlog++);
                               MPI_Error_string(statuses[0].MPI_ERROR,msg_err,&reslen);
                               msg_err[reslen++] = ','; msg_err[reslen]= 0;
                               putlog(msg_err,numlog++);
                               }    else numlog++;
   req_numR = 0;
  if(pr_neighbour[4]>-1) CopyBufferToGrid(f,nut,buf_recv[4],0,0,0,m1-1,m2-1,ghost-1);
  if(pr_neighbour[5]>-1) CopyBufferToGrid(f,nut,buf_recv[5],0,0,mm3,m1-1,m2-1,m3-1);

//    MPI_Barrier(MPI_COMM_WORLD);
//    MPI_Startall(req_numR,RecvRequest);
//    MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,statuses);
//    if(flag==0) putlog("bc:error during transfer=",numlog++);
//    MPI_Testall(req_numR,RecvRequest,&flag,statuses);
//    MPI_Get_count(statuses,MPI_DOUBLE,&cnt);
//    MPI_Waitall(req_numR,RecvRequest,statuses);
}

void  boundary_conditions(double ****f, double ***nut)
{
   int i, j, k, l;
   double vrho,vphi,vth;

  /*============================ divertor =================================*/
  if(t_cur<-0.1)                  // divertors are off till t sec
  if(n[1]==0)
  for(i=0;i<m1;i++)
    for(j=ghost;j<=ghost;j++)
      for(k=0;k<m3;k++)
      if(isType(node[i][k],NodeFluid))
         {
	 vrho = f[3][i][j][k]*costh[i][k]+f[1][i][j][k]*sinth[i][k];
         vth  = -f[3][i][j][k]*sinth[i][k]+f[1][i][j][k]*costh[i][k];
         vphi = sqrt(pow(f[2][i][j][k],2)+vth*vth);           //sqrt(vfi*vfi+vth*vth)
         f[1][i][j][k] = vrho*sinth[i][k]+vphi*costh[i][k]*sin(chi[i][k]);
         f[2][i][j][k] = vphi*cos(chi[i][k]);
         f[3][i][j][k] = vrho*costh[i][k]-vphi*sinth[i][k]*sin(chi[i][k]);
/*         f[1][i][j][k] = f[3][i][j][k] = 0;
         f[2][i][j][k] = (R*R-pow(coordin(i,0),2) - pow(coordin(k,2),2))/R/R;*/
         }

   timeE0 = MPI_Wtime();
   interprocessor_communication(f,nut);
   timeE1+=(MPI_Wtime()-timeE0);

   ghost_filling(f,nut);
   if(pp[0]>2 || pp[2]>2)
        {
        timeE0 = MPI_Wtime();
        interprocessor_communication(f,nut);
        timeE1+=(MPI_Wtime()-timeE0);
        }
  return;
}

void  init_conditions()
{
   int i,j,k,l,ll;
   double r, rho, r1, z1;
   int ghost1=(approx-1)/2;
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
        r1 = coordin(i,0);   z1 = coordin(k,2);
        rho=sqrt(r1*r1 + z1*z1);
     //regions
        if(rho>Rsh) setType(&node[i][k],NodeVacuum);
          else if(rho>Rfl) setType(&node[i][k],NodeShell);
          else             setType(&node[i][k],NodeFluid);
     //for hydrodynamics
        if(isType(node[i][k],NodeFluid))
            { if(isType(node[i][k],NodeGhostFluid)) node[i][k] -= NodeGhostFluid;
              for(l=-ghost1;l<=ghost1;l++)
              for(ll=-ghost1;ll<=ghost1;ll++)
                   if(i+l>=0&&i+l<m1 && k+ll>=0&&k+ll<m3)
                        if(!isType(node[i+l][k+ll],NodeFluid))
                          setType(&node[i+l][k+ll],NodeGhostFluid);
                        else;
                   else if(!isType(node[i][k],NodeClued)) nmessage("out of boundary",i,k);
            }

        refr_f[i][k] = r1*(2*Rfl/rho-1);      //physical coordinates
        refz_f[i][k] = z1*(2*Rfl/rho-1);
        refr_f[i][k] = (refr_f[i][k]+R)/dx[0]-0.5-n[0]+ghost;   // simulation indices
        refz_f[i][k] = (refz_f[i][k]+R)/dx[2]-0.5-n[2]+ghost;
        if(fabs(refr_f[i][k]-i)<1 && fabs(refz_f[i][k]-k)<1 && !isType(node[i][k],NodeFluid))
                 { setType(&node[i][k],NodeFluid);
                   if(isType(node[i][k],NodeGhostFluid)) node[i][k] -= NodeGhostFluid;
                   for(l=-ghost1;l<=ghost1;l++)
                     { if(i+l>=0&&i+l<m1) if(!isType(node[i+l][k],NodeFluid))
                                              setType(&node[i+l][k],NodeGhostFluid);
                       if(k+l>=0&&k+l<m3) if(!isType(node[i][k+l],NodeFluid))
                                              setType(&node[i][k+l],NodeGhostFluid);
                     }
                  }
     //for magnetism
        refr_m[i][k] = i;
        refz_m[i][k] = k;
        if(i+n[0]<ghost)  {setType(&node[i][k],NodeGhostMagn);
                           refr_m[i][k] = 2*ghost-i-1;}
         else if(i+n[0]>=N1+ghost)
                          {setType(&node[i][k],NodeGhostMagn);
                           refr_m[i][k] = 2*mm1-i-1;}
        if(k+n[2]<ghost)  {setType(&node[i][k],NodeGhostMagn);
                           refz_m[i][k] = 2*ghost-k-1;}
         else if(k+n[2]>=N3+ghost)
                          {setType(&node[i][k],NodeGhostMagn);
                           refz_m[i][k] = 2*mm3-k-1;}
        if(!isType(node[i][k],NodeGhostMagn)) setType(&node[i][k],NodeMagn);
//                                         else node[i][k] -= NodeGhostMagn;
                                          //deleting fictive cells-for constant outside field
     //for divertor's blade
	sinth[i][k]=r1/rho;
    costh[i][k]=z1/rho;
    chi[i][k]  = chimax*M_PI/180.*rho/Rfl;
       }

   for(i=0;i<m1;i++) { r_1[i] = rc/(rc*coordin(i,0)+1); r_2[i] = r_1[i]*r_1[i]; }
   for(i=0;i<m1;i++)
   for(j=0;j<m2;j++)
   for(k=0;k<m3;k++)
       {
       if(isType(node[i][k],NodeFluid)) eta[i][j][k]=1.;
       if(isType(node[i][k],NodeShell)) eta[i][j][k]=etash;
       if(isType(node[i][k],NodeVacuum)) eta[i][j][k]=etavac;
       for(l=0;l<3;l++) B[l][i][j][k] = 0;
       }

// --------------- initial conditions -----------------------------------------
//   k1=2*M_PI/lfi;  k3=M_PI/l3;

if(!goon) {
   for(i=0;i<m1;i++)
   for(j=0;j<m2;j++)
   for(k=0;k<m3;k++)
     {
      if(isType(node[i][k],NodeFluid)) {
        f[0][i][j][k]=1;
        f[1][i][j][k]=0*Noise*((double)rand()-RAND_MAX/2)/RAND_MAX*
                       (Rfl*Rfl-pow(coordin(i,0),2) - pow(coordin(k,2),2))*4/Rfl/Rfl;
        f[2][i][j][k]=NoiseNorm*cos(2*M_PI*coordin(j,1)/R)*sin(2*M_PI*coordin(k,2)/Rfl)
		      + (parabole+0*Noise*((double)rand()-RAND_MAX/2)/RAND_MAX)*
                       (Rfl*Rfl-pow(coordin(i,0),2) - pow(coordin(k,2),2))*4/Rfl/Rfl;
        f[3][i][j][k]=NoiseNorm*sin(2*M_PI*coordin(j,1)/R)*sin(2*M_PI*coordin(k,2)/Rfl)
                      + 0*Noise*((double)rand()-RAND_MAX/2)/RAND_MAX*
                       (Rfl*Rfl-pow(coordin(i,0),2) - pow(coordin(k,2),2))*4/Rfl/Rfl;
//        f[1][i][j][k]=f[2][i][j][k]=f[3][i][j][k]=0;
        nut[i][j][k]=(
//        (0.39+14.8*exp(-2.13*pow(2*coordin(k,2)-l3,2)))*0.1*0
                    +1.)/Re;
                                        }
            else f[0][i][j][k]=f[1][i][j][k]=f[2][i][j][k]=f[3][i][j][k] = 0;
      if(isType(node[i][k],NodeFluid)/*||isType(node[i][k],NodeShell)*/ )
                 {   }
//                 f[4][i][j][k]=f[5][i][j][k]=f[6][i][j][k]=0;
                 /*f[4][i][j][k]=coordin(i,0)*cos(coordin(j,1))*sin(coordin(j,1))*/;
                 /*f[5][i][j][k]=coordin(i,0)*cos(coordin(j,1))*cos(coordin(j,1))*/;
//                 f[6][i][j][k]=0;
                  if(isType(node[i][k],NodeMagn))
//                  if(isType(node[i][k],NodeFluid))
                     {
                     f[4][i][j][k]=0*coordin(j,1)+Noise*((double)rand()-RAND_MAX/2)/RAND_MAX + NoiseNorm*cos(2*M_PI*coordin(j,1)/R)*sin(2*M_PI*coordin(k,2)/Rfl);
                     f[5][i][j][k]=Noise*((double)rand()-RAND_MAX/2)/RAND_MAX + NoiseNorm*cos(2*M_PI*coordin(j,1)/R)*sin(2*M_PI*coordin(k,2)/Rfl);
                     f[6][i][j][k]=Noise*((double)rand()-RAND_MAX/2)/RAND_MAX + NoiseNorm*cos(2*M_PI*coordin(j,1)/R)*sin(2*M_PI*coordin(k,2)/Rfl);
                     }
                   else
                     {
                     f1[4][i][j][k]=f[4][i][j][k]=0*coordin(j,1);
                     f1[5][i][j][k]=f[5][i][j][k]=0;
                     f1[6][i][j][k]=f[6][i][j][k]=0;
                     }
      }
//   struct_func(f,2,2,3);
   fill_velocity(0.3, f);
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
//         if((mintime<0 || vtime<mintime) && ceil((double)N1/kp1)==floor((double)N1/kp1) && ceil((double)N2/kp2)==floor((double)N2/kp2) && ceil((double)N3/kp3)==floor((double)N3/kp3)) { mintime=vtime; nd1=i; nd2=j; }
         }
if(!goon) {                       //reading sizes from file when continuing
  pp[0]=divisors[nd1]; pp[1]=divisors[nd2]; pp[2]=size/pp[0]/pp[1];                 // number of procs along axes
  if(pp[0]>2) pp[0]=2;
  if(pp[2]>2) pp[2]=2;
  pp[1]=size/pp[0]/pp[2];
          }
  pr[0] = rank%pp[0]; pr[1] = (rank/pp[0])%pp[1]; pr[2] = (rank/pp[0]/pp[1])%pp[2];  // coordinates of current subregion

/* dimensions of subregion:         global indicies of subregion origin          if there's nonequal subregions*/
  n1 = floor((double)N1/pp[0]);     n[0] = n1*pr[0] + min(pr[0],N1-pp[0]*n1);    if(pr[0]<N1-pp[0]*n1) n1++;              // dimensions of subregion
  n2 = floor((double)N2/pp[1]);     n[1] = n2*pr[1] + min(pr[1],N2-pp[1]*n2);    if(pr[1]<N2-pp[1]*n2) n2++;
  n3 = floor((double)N3/pp[2]);     n[2] = n3*pr[2] + min(pr[2],N3-pp[2]*n3);    if(pr[2]<N3-pp[2]*n3) n3++;
if(rank==size-1)   {
   iop=fopen(NameStatFile,"w");
   fprintf(iop,"%d\n",rank);
   fprintf(iop,"%d\t%d\t%d\n",pr[0],pr[1],pr[2]);
   fprintf(iop,"%d\t%d\t%d\n",pp[0],pp[1],pp[2]);
   fprintf(iop,"%d\t%d\t%d\n",n[0],n[1],n[2]);
   fprintf(iop,"%d\t%d\t%d\n",n1,n2,n3);
       }
   if(n1<ghost || n2<ghost || n3<ghost) nrerror("Too small mesh or incorrect number of processes",0,0);

   m1 = n1+2*ghost;
   m2 = n2+2*ghost;
   m3 = n3+2*ghost;
   mm1 = ghost+n1;
   mm2 = ghost+n2;
   mm3 = ghost+n3;

  pr_neighbour[0] = (pr[0]>0       ? rank-1           :-1);                           // neighbours of subregion
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
if(rank==size-1) {
  for(i=0;i<6;i++)

   fprintf(iop,"%d ",pr_neighbour[i]);
   fprintf(iop,"\n");
   fileclose(iop);
   }
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

double dvv(double ***f1, double ***f2, int dir, int ii, int jj, int kk, int sh, int sm)
/* find d(f1 f2)/dx in dir direction*/
{
double tmp;
switch (sm) {
     case 7 :  switch (dir) {
                  case 1 : tmp = kf7[0][sh][6]*(f1[ii+3][jj][kk]*f2[ii+3][jj][kk]-f1[ii-3][jj][kk]*f2[ii-3][jj][kk])
                               + kf7[0][sh][5]*(f1[ii+2][jj][kk]*f2[ii+2][jj][kk]-f1[ii-2][jj][kk]*f2[ii-2][jj][kk])
                               + kf7[0][sh][4]*(f1[ii+1][jj][kk]*f2[ii+1][jj][kk]-f1[ii-1][jj][kk]*f2[ii-1][jj][kk]); break;
                  case 2 : tmp = kf7[0][sh][6]*(f1[ii][jj+3][kk]*f2[ii][jj+3][kk]-f1[ii][jj-3][kk]*f2[ii][jj-3][kk])
                               + kf7[0][sh][5]*(f1[ii][jj+2][kk]*f2[ii][jj+2][kk]-f1[ii][jj-2][kk]*f2[ii][jj-2][kk])
                               + kf7[0][sh][4]*(f1[ii][jj+1][kk]*f2[ii][jj+1][kk]-f1[ii][jj-1][kk]*f2[ii][jj-1][kk]); break;
                  case 3 : tmp = kf7[0][sh][6]*(f1[ii][jj][kk+3]*f2[ii][jj][kk+3]-f1[ii][jj][kk-3]*f2[ii][jj][kk-3])
                               + kf7[0][sh][5]*(f1[ii][jj][kk+2]*f2[ii][jj][kk+2]-f1[ii][jj][kk-2]*f2[ii][jj][kk-2])
                               + kf7[0][sh][4]*(f1[ii][jj][kk+1]*f2[ii][jj][kk+1]-f1[ii][jj][kk-1]*f2[ii][jj][kk-1]); break;
                  }; break;
     case 3 :  switch (dir) {
                  case 1 : tmp = kf3[0][sh][2]*(f1[ii+1][jj][kk]*f2[ii+1][jj][kk]-f1[ii-1][jj][kk]*f2[ii-1][jj][kk]); break;
                  case 2 : tmp = kf3[0][sh][2]*(f1[ii][jj+1][kk]*f2[ii][jj+1][kk]-f1[ii][jj-1][kk]*f2[ii][jj-1][kk]); break;
                  case 3 : tmp = kf3[0][sh][2]*(f1[ii][jj][kk+1]*f2[ii][jj][kk+1]-f1[ii][jj][kk-1]*f2[ii][jj][kk-1]); break;
                  }; break;
     case 5 :  switch (dir) {
                  case 1 : tmp = kf5[0][sh][4]*(f1[ii+2][jj][kk]*f2[ii+2][jj][kk]-f1[ii-2][jj][kk]*f2[ii-2][jj][kk])
                               + kf5[0][sh][3]*(f1[ii+1][jj][kk]*f2[ii+1][jj][kk]-f1[ii-1][jj][kk]*f2[ii-1][jj][kk]); break;
                  case 2 : tmp = kf5[0][sh][4]*(f1[ii][jj+2][kk]*f2[ii][jj+2][kk]-f1[ii][jj-2][kk]*f2[ii][jj-2][kk])
                               + kf5[0][sh][3]*(f1[ii][jj+1][kk]*f2[ii][jj+1][kk]-f1[ii][jj-1][kk]*f2[ii][jj-1][kk]); break;
                  case 3 : tmp = kf5[0][sh][4]*(f1[ii][jj][kk+2]*f2[ii][jj][kk+2]-f1[ii][jj][kk-2]*f2[ii][jj][kk-2])
                               + kf5[0][sh][3]*(f1[ii][jj][kk+1]*f2[ii][jj][kk+1]-f1[ii][jj][kk-1]*f2[ii][jj][kk-1]); break;
                  }; break;
    }
  return tmp;
}

double dr(double ***m, int ii, int jj, int kk, int dir, int or, double dx, int sh,  int sm)
/*        matrix     , point                 , direct, order , differ   , shift , sample */
/*                                           , 1,2,3 ,  0,1    dx,dx^2  , 0-left , 3,5,7 */
{
double tmp=0.0;

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

double d2cross(double ***m, int ii, int jj, int kk, int dir1,int dir2,  int sh, int sm)
/*             matrix     , point                 , directions of dervs,shift , sample */
/*                                                , 1, 2, 3           , 0-left, 3,5,7 */
{         //order ==0 (first),dx=dx[dir1]*dx[dir2]
double tmp=0.0;
int i1,i2;
int dirr=6-dir1-dir2;
switch (sm*dirr) {
	case 3 : for(i1=0; i1<sm; i1++)
		     for(i2=0; i2<sm; i2++)
		       tmp += m[ii][jj+i1-sh][kk+i2-sh]*kf3[0][sh][i1]*kf3[0][sh][i2]; break;
	case 6 : for(i1=0; i1<sm; i1++)
		     for(i2=0; i2<sm; i2++)
		       tmp += m[ii+i1-sh][jj][kk+i2-sh]*kf3[0][sh][i1]*kf3[0][sh][i2]; break;
	case 9 : for(i1=0; i1<sm; i1++)
		     for(i2=0; i2<sm; i2++)
		       tmp += m[ii+i1-sh][jj+i2-sh][kk]*kf3[0][sh][i1]*kf3[0][sh][i2]; break;
	case 5 : for(i1=0; i1<sm; i1++)
		     for(i2=0; i2<sm; i2++)
		       tmp += m[ii][jj+i1-sh][kk+i2-sh]*kf5[0][sh][i1]*kf5[0][sh][i2]; break;
	case 10: for(i1=0; i1<sm; i1++)
		     for(i2=0; i2<sm; i2++)
		       tmp += m[ii+i1-sh][jj][kk+i2-sh]*kf5[0][sh][i1]*kf5[0][sh][i2]; break;
	case 15: for(i1=0; i1<sm; i1++)
		     for(i2=0; i2<sm; i2++)
		       tmp += m[ii+i1-sh][jj+i2-sh][kk]*kf5[0][sh][i1]*kf5[0][sh][i2]; break;
	case 7 : for(i1=0; i1<sm; i1++)
		     for(i2=0; i2<sm; i2++)
		       tmp += m[ii][jj+i1-sh][kk+i2-sh]*kf7[0][sh][i1]*kf7[0][sh][i2]; break;
	case 14: for(i1=0; i1<sm; i1++)
		     for(i2=0; i2<sm; i2++)
		       tmp += m[ii+i1-sh][jj][kk+i2-sh]*kf7[0][sh][i1]*kf7[0][sh][i2]; break;
	case 21: for(i1=0; i1<sm; i1++)
		     for(i2=0; i2<sm; i2++)
		       tmp += m[ii+i1-sh][jj+i2-sh][kk]*kf7[0][sh][i1]*kf7[0][sh][i2]; break;
	default :
	nrerror("\nNO SUCH SAMPLE for derivative. Bye ...",0,0);
	}
return(tmp/dx[dir1-1]/dx[dir2-1]);
}

double coordin(int i, int dir)
		      //0-r,1-phi,2-z
{
 switch (dir)
 {
   case 0:  return dx[dir]*(i-ghost+0.5+n[dir])-R;
   case 1:  return dx[dir]*(i-ghost+0.5+n[dir]);
   case 2:  return dx[dir]*(i-ghost+0.5+n[dir])-R;
 }
    return(0);
}

void calculate_curl(double ****f,double ****b,enum TypeNodes tip)
{
   int i,j,k,l,m;
   double dA[7][3];
   for(l=0;l<=2;l++) for(m=0;m<3;m++) dA[l][m] = 0;
   for(i=ghost;i<mm1;i++)
   for(j=ghost;j<mm2;j++)
   for(k=ghost;k<mm3;k++)
   if(isType(node[i][k],tip) && !isType(node[i][k],NodeClued))
   {
      for(l=0;l<=2;l++) {
       for(m=0;m<3;m++)
         dA[l][m]=dr(f[l],i,j,k,m+1,0,dx[m],ghost, approx);
       }
      b[0][i][j][k] = dA[2][1] - dA[1][2];
      b[1][i][j][k] = dA[0][2] - dA[2][0];
      b[2][i][j][k] = dA[1][0] - dA[0][1] + f[1][i][j][k]*r_1[i];
    }
}
