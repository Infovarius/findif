void extract(double VxOld[][][], double VxNL[][][], double VxRH[][][],
             double VyOld[][][], double VyNL[][][], double VyRH[][][],
             double VzOld[][][], double VzNL[][][], double VzRH[][][],
             double pOld[][][],  double pNL[][][],  double pRH[][][]   )
      //converting next-layer's VxNew and P to right-hand function in integrating
{
int i,j,k;
for(i=0;i<=Nx+1;i++)
   for(j=0;j<=Ny+1;j++)
       for(k=0;k<=Nz+1;k++)
           {
           VxRH[i][j][k] = (VxNL[i][j][k] - VxOld[i][j][k]) / dt;
           VyRH[i][j][k] = (VyNL[i][j][k] - VyOld[i][j][k]) / dt;
           VzRH[i][j][k] = (VzNL[i][j][k] - VzOld[i][j][k]) / dt;
           pRH[i][j][k]  = (pNL[i][j][k]  - pOld[i][j][k])  / dt;
           }
}

void Euler(  double VxOld[][][], double VxNew[][][],
             double VyOld[][][], double VyNew[][][],
             double VzOld[][][], double VzNew[][][],
             double pOld[][][],  double pNew[][][]     )
             //Euler's method of integrating
{
double Vx1[Nx+2][Ny+2][Nz+2],Vy1[Nx+2][Ny+2][Nz+2],
       Vz1[Nx+2][Ny+2][Nz+2],p1[Nx+2][Ny+2][Nz+2];
int i,j,k;

next_layer(VxOld,VxNew,VyOld,VyNew,VzOld,VzNew,pOld,pNew);
extract(&VxOld,&VxNew,&Vx1,&VyOld,&VyNew,&Vy1,&VzOld,&VzNew,&Vz1,&pOld,&pNew,&p1);
for(i=0;i<=Nx+1;i++)
   for(j=0;j<=Ny+1;j++)
       for(k=0;k<=Nz+1;k++)
           {
           VxNew[i][j][k] = VxOld[i][j][k] + Vx1[i][j][k] * dt;
           VyNew[i][j][k] = VyOld[i][j][k] + Vy1[i][j][k] * dt;
           VzNew[i][j][k] = VzOld[i][j][k] + Vz1[i][j][k] * dt;
           pNew[i][j][k]  = pOld[i][j][k]  + p1[i][j][k]  * dt;
           }
}

void RungeKutt( double VxOld[][][], double VxNew[][][],
                double VyOld[][][], double VyNew[][][],
                double VzOld[][][], double VzNew[][][],
                double pOld[][][],  double pNew[][][]    )
             //Runge-Kutta's method of integrating
{
double Vx1[Nx+2][Ny+2][Nz+2],Vx2[Nx+2][Ny+2][Nz+2],Vx3[Nx+2][Ny+2][Nz+2],Vx4[Nx+2][Ny+2][Nz+2],
       Vy1[Nx+2][Ny+2][Nz+2],Vy2[Nx+2][Ny+2][Nz+2],Vy3[Nx+2][Ny+2][Nz+2],Vy4[Nx+2][Ny+2][Nz+2],
       Vz1[Nx+2][Ny+2][Nz+2],Vz2[Nx+2][Ny+2][Nz+2],Vz3[Nx+2][Ny+2][Nz+2],Vz4[Nx+2][Ny+2][Nz+2],
       p1[Nx+2][Ny+2][Nz+2],p2[Nx+2][Ny+2][Nz+2],p3[Nx+2][Ny+2][Nz+2],p4[Nx+2][Ny+2][Nz+2];

int i,j,k;

next_layer(VxOld,VxNew,VyOld,VyNew,VzOld,VzNew,pOld,pNew);
extract(VxOld,VxNew,Vx1,VyOld,VyNew,Vy1,VzOld,VzNew,Vz1,pOld,pNew,p1);
for(i=0;i<=Nx+1;i++)
   for(j=0;j<=Ny+1;j++)
       for(k=0;k<=Nz+1;k++)
           {
           Vx1[i][j][k] = VxOld[i][j][k] + Vx1[i][j][k] * dt/2;
           Vy1[i][j][k] = VyOld[i][j][k] + Vy1[i][j][k] * dt/2;
           Vz1[i][j][k] = VzOld[i][j][k] + Vz1[i][j][k] * dt/2;
           p1[i][j][k]  = pOld[i][j][k]  + p1[i][j][k]  * dt/2;
           }
next_layer(Vx1,VxNew,Vy1,VyNew,Vz1,VzNew,p1,pNew);
extract(Vx1,VxNew,Vx2,Vy1,VyNew,Vy2,Vz1,VzNew,Vz2,p1,pNew,p2);
for(i=0;i<=Nx+1;i++)
   for(j=0;j<=Ny+1;j++)
       for(k=0;k<=Nz+1;k++)
           {
           Vx2[i][j][k] = VxOld[i][j][k] + Vx2[i][j][k] * dt/2;
           Vy2[i][j][k] = VyOld[i][j][k] + Vy2[i][j][k] * dt/2;
           Vz2[i][j][k] = VzOld[i][j][k] + Vz2[i][j][k] * dt/2;
           p2[i][j][k]  = pOld[i][j][k]  + p2[i][j][k]  * dt/2;
           }
next_layer(Vx2,VxNew,Vy2,VyNew,Vz2,VzNew,p2,pNew);
extract(Vx2,VxNew,Vx3,Vy2,VyNew,Vy2,Vz2,VzNew,Vz3,p2,pNew,p3);
for(i=0;i<=Nx+1;i++)
   for(j=0;j<=Ny+1;j++)
       for(k=0;k<=Nz+1;k++)
           {
           Vx3[i][j][k] = VxOld[i][j][k] + Vx3[i][j][k] * dt;
           Vy3[i][j][k] = VyOld[i][j][k] + Vy3[i][j][k] * dt;
           Vz3[i][j][k] = VzOld[i][j][k] + Vz3[i][j][k] * dt;
           p3[i][j][k]  = pOld[i][j][k]  + p3[i][j][k]  * dt;
           }
next_layer(Vx3,VxNew,Vy3,VyNew,Vz3,VzNew,p3,pNew);
extract(Vx3,VxNew,Vx4,Vy3,VyNew,Vy4,Vz3,VzNew,Vz4,p3,pNew,p4);
for(i=0;i<=Nx+1;i++)
   for(j=0;j<=Ny+1;j++)
       for(k=0;k<=Nz+1;k++)
           {
           VxNew[i][j][k] = VxOld[i][j][k]+
                (Vx1[i][j][k] + Vx4[i][j][k] + (Vx1[i][j][k]+Vx4[i][j][k]) *2) * dt/6;
           VyNew[i][j][k] = VyOld[i][j][k]+
                (Vy1[i][j][k] + Vy4[i][j][k] + (Vy1[i][j][k]+Vy4[i][j][k]) *2) * dt/6;
           VzNew[i][j][k] = VzOld[i][j][k]+
                (Vz1[i][j][k] + Vz4[i][j][k] + (Vz1[i][j][k]+Vz4[i][j][k]) *2) * dt/6;
           pNew[i][j][k] = pOld[i][j][k]+
                (p1[i][j][k] + p4[i][j][k] + (p1[i][j][k]+p4[i][j][k]) *2) * dt/6;
           }
}

void step_of_time(  double VxOld[][][], double VxNew[][][],
                    double VyOld[][][], double VyNew[][][],
                    double VzOld[][][], double VzNew[][][],
                    double pOld[][][],  double pNew[][][]      )
                 //one of the time steps of calculation
{
Euler(VxOld,VyOld,VzOld,VxNew,VyNew,VzNew,pOld,pNew);
//Runge(VxOld,VyOld,VzOld,VxNew,VyNew,VzNew,pOld,pNew);

} //step_of_time
