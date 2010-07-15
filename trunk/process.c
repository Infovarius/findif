bool konez;
//# include "integrating.c" in perspective
#include "model.c"


# include "printing.c"
double tmmax = 500;

void process(meshArr vx[2],meshArr vy[2],meshArr vz[2],meshArr p[2],meshArr nut)
{

konez=false;
double tm = 0;

long ns = 0;
int n0 = ns%2 ,n1 = (ns+1)%2;
while(tm<tmmax && ! konez)
   {
   n0 = ns%2;
   n1 = (ns + 1)%2;
   next_layer(vx[n0],vx[n1],vy[n0],vy[n1],vz[n0],vz[n1],p[n0],p[n1],nut);
//   nut_by_flux(n0);
   if((int)((tm+dt/2)/STEP)-(int)((tm-dt/2)/STEP))
	 {
	 printing(vx,vy,vz,p,nut,ns,tm,dt);
	 }
   check(vx,vy,vz,p,nut,ns, &dt,&konez);
   if(kbhit()&&getch()=='q') konez=true;

   ns++;
   tm+=dt;
   }


}


