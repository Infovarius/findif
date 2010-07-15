
double gamma=0.01;

void velocitybounder(meshArr Vx, meshArr Vy, meshArr Vz)  //boundary conditions on velocities
{
int i,j,k;
			 //sticking conditions on horizontal surfaces
 for(i=1;i<=Nx;i++)
	  for(j=1;j<=Ny;j++)
	{
	 Vx[i][j][0] = -Vx[i][j][1];
	 Vy[i][j][0] = -Vy[i][j][1];
	 Vz[i][j][0] = -Vz[i][j][1];

	 Vx[i][j][Nz+1] = -Vx[i][j][Nz];
	 Vy[i][j][Nz+1] = -Vy[i][j][Nz];
	 Vz[i][j][Nz+1] = -Vz[i][j][Nz];
	}
			 //periodic conditions on vertical surfaces
  for(i=1;i<=Nx;i++)
	  for(k=1;k<=Nz;k++)
	 {
	 Vx[i][0][k] = Vx[i][Ny][k];
	 Vy[i][0][k] = Vy[i][Ny][k];
	 Vz[i][0][k] = Vz[i][Ny][k];

	 Vx[i][Ny+1][k] = Vx[i][1][k];
	 Vy[i][Ny+1][k] = Vy[i][1][k];
	 Vz[i][Ny+1][k] = Vz[i][1][k];
	 }

			 //periodic conditions on stream surfaces
	for(j=1;j<=Ny;j++)
		for(k=1;k<=Nz;k++)
	 {
	 Vx[0][j][k] = Vx[Ny][j][k];
	 Vy[0][j][k] = Vy[Ny][j][k];
	 Vz[0][j][k] = Vz[Ny][j][k];

	 Vx[Ny+1][j][k] = Vx[1][j][k];
	 Vy[Ny+1][j][k] = Vy[1][j][k];
	 Vz[Ny+1][j][k] = Vz[1][j][k];
	 }
} //velocitybounder

void next_layer(  meshArr VxOld, meshArr VxNew,
                  meshArr VyOld, meshArr VyNew,
                  meshArr VzOld, meshArr VzNew,
                  meshArr pOld,  meshArr pNew,
                  meshArr nut     )
     //calculating next velocities
{
 int  i,j,k ;
 double diver[Nx+1][Ny+1][Nz+1];

/*------------ Step 0 Boundary condition --------------------*/

velocitybounder(VxOld, VyOld, VzOld);

/*------------ Step I (velocities without pressure) ------------------*/
for(i=1;i<=Nx;i++)
	for(j=1;j<=Ny;j++)
		for(k=1;k<=Nz;k++)
			{
	 VxNew[i][j][k] = VxOld[i][j][k] +
		dt*( ( (nut[i+1][j][k]-nut[i-1][j][k])/(2*dx)-VxOld[i][j][k] )*
								(VxOld[i+1][j][k]-VxOld[i-1][j][k])/(2*dx)+
			  ( (nut[i][j+1][k]-nut[i][j-1][k])/(2*dy)-VyOld[i][j][k] )*
								(VxOld[i][j+1][k]-VxOld[i][j-1][k])/(2*dy)+
			  ( (nut[i][j][k+1]-nut[i][j][k-1])/(2*dz)-VzOld[i][j][k] )*
								(VxOld[i][j][k+1]-VxOld[i][j][k-1])/(2*dz)+
			  ((VxOld[i+1][j][k]-2*VxOld[i][j][k]+VxOld[i-1][j][k])/(dx*dx)+
				(VxOld[i][j+1][k]-2*VxOld[i][j][k]+VxOld[i][j-1][k])/(dy*dy)+
				(VxOld[i][j][k+1]-2*VxOld[i][j][k]+VxOld[i][j][k-1])/(dz*dz))*
						  (nut[i][j][k]) );
	 VyNew[i][j][k] = VyOld[i][j][k] +
		dt*( ( (nut[i+1][j][k]-nut[i-1][j][k])/(2*dx)-VxOld[i][j][k] )*
								(VyOld[i+1][j][k]-VyOld[i-1][j][k])/(2*dx)+
			  ( (nut[i][j+1][k]-nut[i][j-1][k])/(2*dy)-VyOld[i][j][k] )*
								(VyOld[i][j+1][k]-VyOld[i][j-1][k])/(2*dy)+
			  ( (nut[i][j][k+1]-nut[i][j][k-1])/(2*dz)-VzOld[i][j][k] )*
								(VyOld[i][j][k+1]-VyOld[i][j][k-1])/(2*dz)+
			  ((VyOld[i+1][j][k]-2*VyOld[i][j][k]+VyOld[i-1][j][k])/(dx*dx)+
				(VyOld[i][j+1][k]-2*VyOld[i][j][k]+VyOld[i][j-1][k])/(dy*dy)+
				(VyOld[i][j][k+1]-2*VyOld[i][j][k]+VyOld[i][j][k-1])/(dz*dz))*
						  ( nut[i][j][k]) );
	 VzNew[i][j][k] = VzOld[i][j][k] +
		dt*( ( (nut[i+1][j][k]-nut[i-1][j][k])/(2*dx)-VxOld[i][j][k] )*
								(VzOld[i+1][j][k]-VzOld[i-1][j][k])/(2*dx)+
			  ( (nut[i][j+1][k]-nut[i][j-1][k])/(2*dy)-VyOld[i][j][k] )*
								(VzOld[i][j+1][k]-VzOld[i][j-1][k])/(2*dy)+
			  ( (nut[i][j][k+1]-nut[i][j][k-1])/(2*dz)-VzOld[i][j][k] )*
								(VzOld[i][j][k+1]-VzOld[i][j][k-1])/(2*dz)+
			  ((VzOld[i+1][j][k]-2*VzOld[i][j][k]+VzOld[i-1][j][k])/(dx*dx)+
				(VzOld[i][j+1][k]-2*VzOld[i][j][k]+VzOld[i][j-1][k])/(dy*dy)+
				(VzOld[i][j][k+1]-2*VzOld[i][j][k]+VzOld[i][j][k-1])/(dz*dz))*
						  (nut[i][j][k]) );
			}

/*------------ Step I Boundary condition --------------------*/

velocitybounder(VxNew, VyNew, VzNew);

/*------------ Step II Divergention---------------------------*/

  for(i=1;i<=Nx;i++)
	 for(j=1;j<=Ny;j++)
		 for(k=1;k<=Nz;k++)
		 {
		 diver[i][j][k] = (VxNew[i+1][j][k]-VxNew[i-1][j][k])/(2*dx)+
				 (VyNew[i][j+1][k]-VyNew[i][j-1][k])/(2*dy)+
				 (VzNew[i][j][k+1]-VzNew[i][j][k-1])/(2*dz);
		 }

/*------------ Step III (counting pressure)---------------------------*/

	for(i=1;i<=Nx;i++)
	  for(j=1;j<=Ny;j++)
			for(k=1;k<=Nz;k++)
		  pNew[i][j][k] = pOld[i][j][k] - dt*diver[i][j][k]/gamma;

/*------------ Step IV (pressure to velocities)-----------------------*/

			  //free conditions on horizontal surfaces
	for(i=1;i<=Nx;i++)
	  for(j=1;j<=Ny;j++)
		{
	  pNew[i][j][0] = pNew[i][j][1];
	  pNew[i][j][Nz+1] = pNew[i][j][Nz];
		}

				//periodic conditions on vertical surfaces
	for(i=1;i<=Nx;i++)
	  for(k=0;k<=Nz+1;k++)
		{
	  pNew[i][0][k] = pNew[i][Ny][k];
	  pNew[i][Ny+1][k] = pNew[i][1][k];
		}

				//gradient-periodic conditions on stream surfaces
	for(j=0;j<=Ny+1;j++)
	  for(k=0;k<=Nz+1;k++)
	 {
	  pNew[0][j][k]= pNew[Nx][j][k]-(p2-p1);
	  pNew[Nx+1][j][k]= pNew[1][j][k]+(p2-p1);
	 }

//correction of velocities by pressure
	for(i=1;i<=Nx;i++)
	  for(j=1;j<=Ny;j++)
			for(k=1;k<=Nz;k++)
			  {
		 VxNew[i][j][k] -= dt*(pNew[i+1][j][k]-pNew[i-1][j][k])/(2*dx);
		 VyNew[i][j][k] -= dt*(pNew[i][j+1][k]-pNew[i][j-1][k])/(2*dy);
		 VzNew[i][j][k] -= dt*(pNew[i][j][k+1]-pNew[i][j][k-1])/(2*dz);
			  }

}  //next_layer
