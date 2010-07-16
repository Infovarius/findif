//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#define LEVEL extern
#include "run.h"
#include "head.h"
#include "kaskad.h"
#pragma package(smart_init)
//---------------------------------------------------------------------------
//   Important: Methods and properties of objects in VCL can only be
//   used in a method called using Synchronize, for example:
//
//      Synchronize(UpdateCaption);
//
//   where UpdateCaption could look like:
//
//      void __fastcall Unit1::UpdateCaption()
//      {
//        Form1->Caption = "Updated in a thread";
//      }
//---------------------------------------------------------------------------

__fastcall CalcProcess::CalcProcess(bool CreateSuspended)
        : TThread(CreateSuspended)
{
 Priority = tpIdle;
}
//---------------------------------------------------------------------------
void __fastcall CalcProcess::Execute()
{
   double dttry, dtdid, dtnext, PulsEn;
   int i,j,k,l;
   const min_dump=200;
   AnsiString TempScreen;
   MainWindow->ChangeStatus("Идут вычисления...","Начаты вычисления");
   MainWindow->BeginWork->Caption="Пауза";
   time(&time_begin);
   t_cur = 0;
   count = 0;  enter = 0;
   nmessage("work has begun",t_cur);
   PulsEn=check(f);
   printing(f,dtdid,t_cur,count,PulsEn);

   dtnext=1e-3;
   dump(n1,n2,n3,Re,f,nut,t_cur,count);
   action = ActNone;

        /*----------- MAIN ITERATIONS --------------*/
   while (t_cur < Ttot && !razlet && action!=ActQuit) {

        pde(t_cur, f, df);
        dttry=dtnext;
        timestep(f, df, t_cur, f1, dttry, &dtdid, &dtnext);
        nut_by_flux(f,dtdid);
        if(int(min_dump*t_cur/Ttot)-int(min_dump*(t_cur+dtdid)/Ttot))
                   dump(n1,n2,n3,Re,f,nut,t_cur,count);
        t_cur+=dtdid;
        count++;
        if (count%CheckStep==0)
            PulsEn=check(f);
        if (count%OutStep==0)
            {
            if (count%CheckStep!=0)
                PulsEn=check(f);
            printing(f1,dtdid,t_cur,count,PulsEn);
            };
        for(l=0;l<nvar;l++)
        for(i=ghost;i<mm1;i++)
        for(j=ghost;j<mm2;j++)
        for(k=ghost;k<mm3;k++)
           f[l][i][j][k]=f1[l][i][j][k];
        if(action)
             {
                switch (action) {
                        case ActDump  : { dump(n1,n2,n3,Re,f,nut,t_cur,count);
                                          action = ActNone;}       break;
                        case ActQuit  : { dump(n1,n2,n3,Re,f,nut,t_cur,count);
                                          nrerror("You asked to exit. Here you are...",t_cur);
                                          TempScreen = "You asked to exit. Here you are...";
                                          MainWindow->Screen->Lines->Append(TempScreen);
                                          RunStatus = StNone;
                                          MainWindow->BeginWork->Caption = "Начать счет";
                                          return;
                                   }; break;
                    /*    case ActPause : { MainWindow->BeginWork->Caption="Продолжить счет";
                                          RunStatus = StPause;
                                          MainWindow->ButtonRecontinue->Show();
                                          MainWindow->ChangeStatus("Идут вычисления...","Вычисления приостановлены");
                                          do ; while(action!=ActContinue && action!=ActBegin);
                                          MainWindow->BeginWork->Caption="Пауза";
                                          MainWindow->ButtonRecontinue->Hide();
                                          MainWindow->ChangeStatus("Идут вычисления...","Продолжение вычислений");
                                          RunStatus = StRun;
                                          action = ActNone;
                                        }; break;*/
                        }
             }
   }
 MainWindow->Screen->Lines->Append("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
 MainWindow->Screen->Lines->Append(" ");

  if(t_cur>Ttot&&!razlet)
       { nmessage("work is succesfully done",t_cur);
         TempScreen.printf("work is succesfully done at modelling time=%g",t_cur);
         MainWindow->Screen->Lines->Append(TempScreen);}
       else
         { nrerror("this is break of scheme",t_cur);
           TempScreen.printf("this is break of scheme at modelling time=%g",t_cur);
           dump(n1,n2,n3,Re,f,nut,t_cur,count);
           MainWindow->Screen->Lines->Append(TempScreen);}
 at_the_end_of_work();          
}
//---------------------------------------------------------------------------
