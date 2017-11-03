#define LEVEL extern

#include "head.h"

#define MaxTicks 100
#define TicksOutput 20

struct
       {long double sum_time,begin;
        long kol;
        char active;
        char name[20];
       }
       Tick[MaxTicks];

void init_tick(int num_stage,char *inname)
{
 strcpy(Tick[num_stage].name, inname);
 Tick[num_stage].active = 0;
}

inline void start_tick(int num_stage)
{
// strcpy(Tick[num_stage].name, inname);
 Tick[num_stage].active = 1;
 Tick[num_stage].begin = MPI_Wtime();
}

inline void finish_tick(int num_stage)
{
// if(!Tick[num_stage].active) return;
 Tick[num_stage].sum_time += MPI_Wtime() - Tick[num_stage].begin;
 Tick[num_stage].kol++;
 Tick[num_stage].active = 0;
}

void init_timers()
{
 int i;
 for(i=0;i<MaxTicks;i++)
    {
    Tick[i].kol = 0;
    Tick[i].active = 0;
    Tick[i].sum_time = 0;
    }
 init_tick(0,"0.init");
 init_tick(1,"1.BC_inside");
 init_tick(2,"2.none");
 init_tick(3,"3.pde_inside");
 init_tick(4,"4.timestep_inside");
 init_tick(5,"5.nut_inside");
 init_tick(6,"6.none");
 init_tick(7,"7.none");
 init_tick(8,"8.none");
 init_tick(9,"9.dr_inside");
 init_tick(10,"10.d2cr_inside");
 init_tick(11,"11.1st_pde_wrap");
 init_tick(12,"12.timestep_wrap");
 init_tick(13,"13.rest_cycle");
 init_tick(14,"14.array_exch");
 init_tick(15,"15.pde_dervs");
 init_tick(16,"16.pde_eqs");
 init_tick(17,"17.1st_derv");
 init_tick(18,"18.2nd_derv");
 init_tick(19,"19.mix_derv");
 init_tick(20,"20.none");
}

void print_CPU_usage()
{
 int i;
 double sum1=0,sum2=0;
 FILE *f;
 f = fileopen("CPUusage.dat",1);
 for(i=0;i<MaxTicks;i++)
   sum1 += Tick[i].sum_time;
 sum2 = MPI_Wtime()-time_begin;
 fprintf(f,"\n Sum of times  = %9.4lf;\n Total runtime = %9.4lf\n",sum1,sum2);
 for(i=0;i<=TicksOutput;i++)
   {
   if(Tick[i].name[0]) fprintf(f,"%s timer (per %d times) =\t",Tick[i].name,Tick[i].kol);
                  else fprintf(f,"%2dth timer (per %d times) =\t",i,Tick[i].kol);
   fprintf(f,"%10.4Lf=\t",Tick[i].sum_time);
   fprintf(f,"%9.4Lf%%",100*Tick[i].sum_time/sum2);
   fprintf(f,"\n");
   }
 fileclose(f);
}
