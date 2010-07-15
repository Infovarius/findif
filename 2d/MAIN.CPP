# include <stdio.h>
# include <conio.h>
# include <stdlib.h>
# include <math.h>

double dt = 1e-4;                     

# include "structur.c"
# include "init.c"
# include "process.c"

main()
{
meshArr  vx[2],vy[2],vz[2],p[2],nut;          //defining

  init(vx[0],vy[0],vz[0],p[0],nut);
  process(vx,vy,vz,p,nut);

  printf("\n The End \n");
  putch(getch());
  printf("\n vx[%u][%u][k]:\n",Nx/2,Ny/2);
  for(int k=1;k<=Nz;k++)
	printf("%g\n",vx[0][Nx/2][Ny/2][k]);
  putch(getch());

}