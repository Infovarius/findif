#define LEVEL extern

#include "head.h"

#define MaxTicks 100
struct
       {long double sum_time,begin;
        long kol;
        char active;
        char name[15];
       }
       Tick[MaxTicks];

void start_tick(int num_stage,char *inname)
{
 strcpy(Tick[num_stage].name, inname);
 Tick[num_stage].active = 1;
 Tick[num_stage].begin = MPI_Wtime();
}

void finish_tick(int num_stage)
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
 for(i=0;i<=30;i++)
   {
   if(Tick[i].name[0]) fprintf(f,"%s timer (per %d times) =\t",Tick[i].name,Tick[i].kol);
                  else fprintf(f,"%2dth timer (per %d times) =\t",i,Tick[i].kol);
   fprintf(f,"%10.4Lf=\t",Tick[i].sum_time);
   fprintf(f,"%9.4Lf%%",100*Tick[i].sum_time/sum2);
   fprintf(f,"\n");
   }
 fileclose(f);
}