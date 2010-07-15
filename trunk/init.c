# include "structur.c"
# include "stdlib.h"

void init(meshArr vx, meshArr vy, meshArr vz, meshArr p, meshArr nut)
{
int i,j,k;
  //parabole profil
for(i=0;i<=Nx+1;i++)
    for(j=0;j<=Ny+1;j++)
    	for(k=0;k<=Nz+1;k++)
            {
	      vx[i][j][k]=Noise*((double)rand()-RAND_MAX/2)/RAND_MAX
	      			//+ 1 - 4./Nz/Nz*(k-1-Nz/2)*(k-1-Nz/2)
                                                             ;
	      vy[i][j][k]=Noise*((double)rand()-RAND_MAX/2)/RAND_MAX;
	      vz[i][j][k]=Noise*((double)rand()-RAND_MAX/2)/RAND_MAX;
	      p[i][j][k] = p1+(i-0.5)*(p2-p1)/Nx;
	      nut[i][j][k]= 1./Re;
            }
}

# undef p1
# undef p2
# undef Noise
