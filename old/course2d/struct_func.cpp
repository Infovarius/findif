#define LEVEL extern

#include "head.h"

/*extern double lx,ly,lz,
		 vx[2][Nx+2][Ny+2][Nz+2],
		 vy[2][Nx+2][Ny+2][Nz+2],
		 vz[2][Nx+2][Ny+2][Nz+2];*/
FILE *fsf;
#define min_d min(min(dx[0],dx[1]),dx[2])

void struct_func(double ****f,int q)
       //structural function of order q and multiplier lambda getting from max_okr^3 cube
       // are putting in array s_func
{
int i,j,k,l,m,n,dist,wrong=0;
double d,rx,ry,rz,logl;
logl=log(lambda);
for(k=ghost;k<mm3;k+=1)
	{
        for(i=0;i<=kol_masht-1;i++) s_func[k-ghost][i] = 0;
        for(i=ghost;i<mm1;i++)
        for(j=ghost;j<=ghost;j++)
                for(l=-max_okr;l<=max_okr;l++)
                for(m=-max_okr;m<=max_okr;m++)
		for(n=-max_okr;n<=max_okr;n++)
			 {
			 if ((l==0)&&(m==0)&&(n==0)) continue;
			 rx = abs(l)*dx[0];
                         ry = abs(m)*dx[1];
                         rz = abs(n)*dx[2];
                         d = sqrt(rx*rx+ry*ry+rz*rz);
                         d = log(d/min_d)/logl;    //(per wave number add "minus")
                         dist=floor(d);
                         if(dist>=kol_masht || dist<0) { wrong++;continue;}
                         s_func[k-ghost][dist]+=norma(f[0][i+l][j+m][k+n]-f[0][i][j][k],
                                                f[1][i+l][j+m][k+n]-f[1][i][j][k],
                                                f[2][i+l][j+m][k+n]-f[2][i][j][k],
                                                q);
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
