# include "mesh.h"
# include "vars.h"

# define min_d min(dx,dy)

# define M 2         // number of edge structural points
# define koef 1      //koefficient in Fourier transform
double sinkx[M][Nx],sinky[M][Ny];
double vx[Nx+2][Ny+2],
       vy[Nx+2][Ny+2],
       vz[Nx+2][Ny+2];
double fvx[M],fvy[M],fvz[M];
double mul[M] = {1,1};     //dissipation parameters

void fourier()      //fourier transform with dissipation
{
double klim=M_PI/2/min_d;
int i,j,k;
//initial filling of sine tables
for(i=0;i<=M-1;i++)
   {
	for(j=0;j<=Nx-1;j++) sinkx[i][j] = sin(klim /(i+1)*(j+1)*dx);
	for(j=0;j<=Ny-1;j++) sinky[i][j] = sin(klim /(i+1)*(j+1)*dy);
   };
//direct transform
for(i=0;i<=M-1;i++)
   {
   fvx[i] = fvy[i] = fvz[i] = 0;
   for(j=1;j<=Nx;j++)
   	for(k=1;k<=Ny;k++)
      	{
      	fvx[i] += vx[j][k]*sinkx[i][j]*sinky[i][k];
      	fvy[i] += vy[j][k]*sinkx[i][j]*sinky[i][k];
      	fvz[i] += vz[j][k]*sinkx[i][j]*sinky[i][k];
         };
   };
//inverse transform
for(j=1;j<=Nx;j++)
  	for(k=1;k<=Ny;k++)
       for(i=0;i<=M-1;i++)
           {
           vx[j][k] += (mul[i]-1)*koef*fvx[i]*sinkx[i][j]*sinky[i][k];
           vy[j][k] += (mul[i]-1)*koef*fvy[i]*sinkx[i][j]*sinky[i][k];
           vz[j][k] += (mul[i]-1)*koef*fvz[i]*sinkx[i][j]*sinky[i][k];
           };
}
