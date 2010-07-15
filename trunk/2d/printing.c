const char *NameNuFile = "nut.dat";
const char *NameVFile  = "vv.dat";

FILE *fileopen(const char *x, int mode)  //opening of file to ff
{
FILE *ff;
if(mode) ff=fopen(x,"a");
else if ((ff = fopen (x,"w+"))==NULL)
	 {
		printf ("Can't open file %s !\n",x);
		exit(-1);
	 }
return(ff);
}

# define min_d min(min(dx,dy),dz)
# define maxvel  100
# define SAFETY  0.2

void check(meshArr vx[2],meshArr vy[2],meshArr vz[2],meshArr p[2],meshArr nut,int ns,
                double *dt,bool *konez)
//whether is unstability and correcting time step
{
 int i,j,k;
// double epsp=0, epsvx=0, epsvy=0, epsvz=0;
 double vxmax=0, vymax=0, vzmax=0, pmax=0, numax=0;
 int n = ns%2;
 for(i=1;i<=Nx;i++)
    for(j=1;j<=Ny;j++)
        for(k=1;k<=Nz;k++)
            {
            vxmax = max(vxmax, fabs(vx[n][i][j][k])                );
            vymax = max(vymax, fabs(vy[n][i][j][k])                );
            vzmax = max(vzmax, fabs(vz[n][i][j][k])                );
            if(i<Nx)
              pmax  = max(pmax,  fabs(p[n][i][j][k]-p[n][i+1][j][k]));
            numax = max(numax, fabs(nut[i][j][k])                  );
 //           epsvx = max(epsvx, fabs(vx[1][i][j][k]-vx[0][i][j][k]) );
 //           epsvy = max(epsvy, fabs(vy[1][i][j][k]-vy[0][i][j][k]) );
 //           epsvz = max(epsvz, fabs(vz[1][i][j][k]-vz[0][i][j][k]) );
            }

//finding unstability
	 if(vxmax>maxvel || vymax>maxvel ||vzmax>maxvel)
			*konez=true;
//correction time step
         if(*dt>SAFETY*min_d/modul(vxmax,vymax,vzmax))
               *dt = SAFETY*min_d/modul(vxmax,vymax,vzmax);     //on nonlinear term
         if(*dt>modul(vxmax,vymax,vzmax)*lx/pmax)
               *dt = SAFETY*modul(vxmax,vymax,vzmax)*dx/pmax;   //on pressure gradient
         if(*dt>pow(min_d,2)/numax)
               *dt = SAFETY*pow(min_d,2)/numax;                 //on viscosity

         if(*dt<1e-20) *konez=true;
         if(*dt>0.1) *dt = 0.1;
}//check

#undef min_d
#undef SAFETY
#undef maxvel

void printing(meshArr vx[2],meshArr vy[2],meshArr vz[2],meshArr p[2],meshArr nut,
                int ns,double tm,double dt)
{
 FILE  *fv, *fnu, *fsf;

 int i,j,k;
 int n0 = (ns + 1)%2 , ind = ns%2;
 double epsp=0, epsvx=0, epsvy=0, epsvz=0;
 double vxmax=0, vymax=0, vzmax=0;
 double avervx[Nz+2], avernu[Nz+2];

   for(i=1;i<=Nx;i++)
      for(j=1;j<=Ny;j++)
	  for(k=1;k<=Nz;k++)
	       {
		 if(vxmax < fabs(vx[ind][i][j][k])) vxmax = fabs(vx[ind][i][j][k]);
		 if(vymax < fabs(vy[ind][i][j][k])) vymax = fabs(vy[ind][i][j][k]);
		 if(vzmax < fabs(vz[ind][i][j][k])) vzmax = fabs(vz[ind][i][j][k]);
		 if(epsvx < fabs(vx[ind][i][j][k]-vx[n0][i][j][k]))
		       epsvx = fabs(vx[ind][i][j][k]-vx[n0][i][j][k]);
		 if(epsvy < fabs(vy[ind][i][j][k]-vy[n0][i][j][k]))
		       epsvy = fabs(vy[ind][i][j][k]-vy[n0][i][j][k]);
		 if(epsvz < fabs(vz[ind][i][j][k]-vz[n0][i][j][k]))
		       epsvz = fabs(vz[ind][i][j][k]-vz[n0][i][j][k]);
		 if(epsp  < fabs(p[ind][i][j][k]-p[n0][i][j][k]))
		       epsp  = fabs(p[ind][i][j][k]-p[n0][i][j][k]);
	      }

for(k=1;k<=Nz;k++)
	{
	  avervx[k] = avernu[k] = 0;
	  for(i=1;i<=Nx;i++)
	      for(j=1;j<=Ny;j++)
		 {
		   avervx[k] += vx[ind][i][j][k];
                   avernu[k] += nut[i][j][k];
		  }
	   avervx[k] /= Nx*Ny;
	   avernu[k] /= Nx*Ny;
	 }

clrscr();
printf("\r %9f    %e %e %e %e\n",tm,epsvx,epsvy,epsvz,epsp);
printf("              %e %e %e %e\n",vxmax,vymax,vzmax,dt);

//putting velocities to file
        fv = fileopen(NameVFile,ns);
	fprintf(fv,"{");
	for(k=1;k<=Nz-1;k++)
		fprintf(fv,"%0.5f,",vx[ind][Nx/2][Ny/2][k]);
	fprintf(fv,"%0.5f}\n",vx[ind][Nx/2][Ny/2][Nz]);
        fclose(fv);
//putting viscosities to file
        fnu = fileopen(NameNuFile,ns);
	fprintf(fnu,"{");
	for(k=1;k<=Nz-1;k++)
		fprintf(fnu,"%0.10f,",avernu[k]);
	fprintf(fnu,"%0.10f}\n",avernu[Nz]);
        fclose(fnu);
	//putting structural functions to file
/*	struct_func(2,ind);
	fprintf(fsf,"{%0.6f",s_func[1][0]);
	for(k=1;k<=Nx/2+Ny/2+Nz-1&&num_points[k];k++)
		fprintf(fsf,",%0.6f",s_func[1][k]);
	fprintf(fsf,"}\n");*/

 //it's a pitty to lose calculated eps's and max's!
} //printing


