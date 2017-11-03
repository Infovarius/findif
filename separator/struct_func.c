#define LEVEL extern

#include "head.h"

/*extern double lx,ly,lz,
		 vx[2][Nx+2][Ny+2][Nz+2],
		 vy[2][Nx+2][Ny+2][Nz+2],
		 vz[2][Nx+2][Ny+2][Nz+2];*/
FILE *fsf;

void struct_func(double ****f,int q,double lambda,int size_okr)
       //structural function of order q and multiplier lambda getting from size_okr^3 cube
       // are putting in array s_func
{
int i,j,k,l,m,n,in,jn,dist;
double d,rx,ry,rz,logl;
logl=log(lambda);
for(k=ghost;k<mm3;k+=1)
	{
        for(i=0;i<=kol_masht-1;i++) s_func[k-ghost][i] = num_points[i] = 0;
        for(i=ghost;i<mm1;i++)
        for(j=ghost;j<mm2;j++)
                for(l=-size_okr;l<=size_okr;l++)
                for(m=-size_okr;m<=size_okr;m++)
		for(n=-min(size_okr,k-ghost);n<=min(size_okr,mm3-1-k);n++)
			 {
			 if ((l==0)&&(m==0)&&(n==0)) continue;
                         in = i+l; jn = j+m;
                         if (i+l<ghost) in += mm1-1-ghost;
                         if (j+m<ghost) jn += mm2-1-ghost;
                         if (i+l>=mm1) in -= mm1-1-ghost;
                         if (j+m>=mm2) jn -= mm2-1-ghost;
			 rx = abs(l)*dx[0];
                         ry = abs(m)*dx[1];
                         rz = abs(n)*dx[2];
                         d = sqrt(rx*rx+ry*ry+rz*rz);
                         d = log(d/min_d)/logl;    //(per wave number add "minus")
                         dist=floor(d);
                         s_func[k-ghost][dist]+=norma(f[0][in][jn][k+n]-f[0][i][j][k],
                                                f[1][in][jn][k+n]-f[1][i][j][k],
                                                f[2][in][jn][k+n]-f[2][i][j][k],
                                                q);
                         num_points[dist]++;
                         }
        for(i=0;i<=kol_masht-1&&num_points[i];i++)
                s_func[k-ghost][i]=(s_func[k-ghost][i]/num_points[i]);
	}//"for" per layers
}//struct_func

/*void main()
{

//fsf = prepout("strfunc.dat");
//fprintf(fsf,"{");

long i,j,k;

 double lx = 2, lz = 1, ly= 1;
 double dx = lx/n1; dy = ly/Ny; dz=lz/Nz;
  //initial filling of arrays
  for(i=0;i<=Nx+1;i++)
	for(j=0;j<=Ny+1;j++)
		for(k=0;k<=Nz+1;k++)
				 {
				 vx[0][i][j][k]=Noise*((double)rand()-RAND_MAX/2)/RAND_MAX
								+ 1 - 4./Nz/Nz*(k-(1+Nz)/2.)*(k-(1+Nz)/2.);
				 vy[0][i][j][k]=Noise*((double)rand()-RAND_MAX/2)/RAND_MAX;
				 vz[0][i][j][k]=Noise*((double)rand()-RAND_MAX/2)/RAND_MAX;
				 }
	 struct_func(2,0);
//fprintf(fsf,"}");
fclose(fsf);
 }*/

