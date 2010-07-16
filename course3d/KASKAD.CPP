//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#define LEVEL
#include "kaskad.h"
#include "head.h"
#include "run.h"
#include <time.h>
#include <conio.h>
#include <system.hpp>
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TMainWindow *MainWindow;
CalcProcess *proc;
TThreadPriority prioritet;
bool input_done=false;
//---------------------------------------------------------------------------
__fastcall TMainWindow::TMainWindow(TComponent* Owner)
        : TForm(Owner)
{
 ChangeStatus("","Запущена программа, спасибо");

}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::FormCreate(TObject *Sender)
{
   Re=10000;
   maschtab=1000;
   Gamma=1e-3;
   l1=3;                    CheckStep=10;
   l2=1;                    OutStep = 20;
   l3=1;

   nvar=4;
   n1=5;
   n2=5;
   n3=11;
   Ns=25;
   approx=7;                      //derivatives approximation order
   ghost=(approx-1)/2;            //radius of approx sample
   dx[0]=l1/n1;
   dx[1]=l2/n2;
   dx[2]=l3/n3;
   p1 = 8*l1/(l3*Re) ; p2 = 0;


   m1 = n1+2*ghost;
   m2 = n2+2*ghost;
   m3 = n3+2*ghost;
   mm1 = ghost+n1;
   mm2 = ghost+n2;
   mm3 = ghost+n3;

   t_cur=0; Ttot=100;
   count=0; enter = 0;

   NameNuFile = "nut.dat";
   NameVFile  = "vv.dat";
   NameDumpFile = "dump.dat";
   NameEnergyFile = "energy.dat";
   NameKaskadVarFile = "kaskvar.dat";
   NameMessageFile = "message.dat";
   NameErrorFile = "error.err";
   NameInitFile = "init.dat";

   s_func = alloc_mem_2f(n3+2,kol_masht);
   f  =alloc_mem_4f(nvar, m1, m2, m3);   //f[3]-pressure,f[0..2]-v(vector)
   f1 =alloc_mem_4f(nvar, m1, m2, m3);
   df =alloc_mem_4f(nvar, m1, m2, m3);
   df2=alloc_mem_4f(nvar, m1, m2, m3);
   df3=alloc_mem_4f(nvar, m1, m2, m3);
   df4=alloc_mem_4f(nvar, m1, m2, m3);
   df5=alloc_mem_4f(nvar, m1, m2, m3);
   nut=alloc_mem_3f(m1, m2, m3);
   init_shell();

   Noise = 0.1;
   NoiseNorm = 0.05;

   prioritet = tpIdle;
   RunStatus = StNone;
   fileopen("error.err",0);
   RefreshForm();
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::RefreshForm()
{
 AnsiString Temp;
 Temp.printf("%g",Re);           EditRe->Text = Temp;
 Temp.printf("%g",maschtab);     EditMaschtab->Text = Temp;
 Temp.printf("%d",OutStep);      EditOutStep->Text = Temp;
 Temp.printf("%d",CheckStep);    EditCheckStep->Text = Temp;
   EditNu->Text = NameNuFile;
   EditV->Text  = NameVFile;
   EditDump->Text = NameDumpFile;
   EditEnergy->Text = NameEnergyFile;
   EditKaskvar->Text = NameKaskadVarFile;
   EditMessage->Text = NameMessageFile;
   EditError->Text = NameErrorFile;
   EditInit->Text = NameInitFile;
 Temp.printf("%g",Noise);        EditNoise->Text = Temp;
 Temp.printf("%g",NoiseNorm);    EditNorm->Text = Temp;
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::FormDestroy(TObject *Sender)
{
   free_mem_2f(s_func,n3+2,kol_masht);
   free_mem_4f(f  ,nvar, m1, m2, m3);
   free_mem_4f(f1 ,nvar, m1, m2, m3);
   free_mem_4f(df ,nvar, m1, m2, m3);
   free_mem_4f(df2,nvar, m1, m2, m3);
   free_mem_4f(df3,nvar, m1, m2, m3);
   free_mem_4f(df4,nvar, m1, m2, m3);
   free_mem_4f(df5,nvar, m1, m2, m3);
   free_mem_3f(nut, m1, m2, m3);
   erase_shell();
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::BeginWorkClick(TObject *Sender)
{
 switch (RunStatus)
 {
   case StNone:  { RunStatus = StRun;
                   ChangeStatus("Загрузка задания...","Сейчас будут начаты вычисления");
                   MainWindow->EditRe->Enabled = false;
                   MainWindow->EditRe->Hint = "Число не может быть изменено во время счета";
                   MainWindow->EditRe->ShowHint = false;
                   action = ActBegin;
                   if(!input_done) init_conditions(f,Re);
                   proc = new CalcProcess(false);
                   proc->Priority = prioritet;
//                   ButtonRecontinue->Visible = false;
                 }
   break;
   case StPause: { RunStatus = StRun;
                   action = ActContinue;
//                   ButtonRecontinue->Visible = false;
                 }
   case StRun  : { RunStatus = StPause;
                   action = ActPause;
//                   ButtonRecontinue->Visible = true;
                 }
   break;
 }
}
// ------------------------Main Block-----------------------------------------
void __fastcall TMainWindow::RecontinueClick(TObject *Sender)
{
   ChangeStatus("Перезагрузка задания...","Вычисления будут начаты сначала");
   RunStatus = StRun;
   action = ActBegin;
   init_conditions(f,Re);
   delete proc;
   proc = new CalcProcess(false);
   proc->Priority = prioritet;
   ChangeStatus("Идут вычисления...","Вычисления начаты сначала");
}
//---------------------------------------------------------------------------
//                           other handlers
void __fastcall TMainWindow::EditMessageKeyUp(TObject *Sender, WORD &Key,
      TShiftState Shift)
{
  NameMessageFile = MainWindow->EditMessage->Text.c_str();
}
//---------------------------------------------------------------------------

void __fastcall TMainWindow::EditEnergyKeyUp(TObject *Sender, WORD &Key,
      TShiftState Shift)
{
  NameEnergyFile = MainWindow->EditEnergy->Text.c_str();
}
//---------------------------------------------------------------------------

void __fastcall TMainWindow::EditVKeyUp(TObject *Sender, WORD &Key,
      TShiftState Shift)
{
  NameVFile = MainWindow->EditV->Text.c_str();
}
//---------------------------------------------------------------------------

void __fastcall TMainWindow::EditNuKeyUp(TObject *Sender, WORD &Key,
      TShiftState Shift)
{
  NameNuFile = MainWindow->EditNu->Text.c_str();
}
//---------------------------------------------------------------------------

void __fastcall TMainWindow::EditKaskvarKeyUp(TObject *Sender, WORD &Key,
      TShiftState Shift)
{
  NameKaskadVarFile = MainWindow->EditKaskvar->Text.c_str();
}
//---------------------------------------------------------------------------

void __fastcall TMainWindow::EditErrorKeyUp(TObject *Sender, WORD &Key,
      TShiftState Shift)
{
  NameErrorFile = MainWindow->EditError->Text.c_str();
}
//---------------------------------------------------------------------------

void __fastcall TMainWindow::EditDumpChange(TObject *Sender)
{
  NameDumpFile = MainWindow->EditDump->Text.c_str();
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::EditInitKeyUp(TObject *Sender, WORD &Key,
      TShiftState Shift)
{
  NameInitFile = MainWindow->EditInit->Text.c_str();
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::InputInitClick(TObject *Sender)
{
  int code=0;
  long dummy_count;
  double dummy_Re,dummy_t;
  if (RadioInitVel->Checked)
         code=init_data(&n1,&n2,&n3,&dummy_Re,f,nut,&dummy_t,&dummy_count);
         else code=init_data(&n1,&n2,&n3,&Re,f,nut,&t_cur,&count);
  RefreshForm();
  if ( code > 0)
     nrerror("Не удался ввод начальных данных. Код ошибки = %d",code);
   else input_done = true;
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::DoDumpClick(TObject *Sender)
{
  dump(n1,n2,n3,Re,f,nut,t_cur,count);
}
//---------------------------------------------------------------------------

void __fastcall TMainWindow::ExitClick(TObject *Sender)
{
 if(RunStatus == StRun || RunStatus == StPause)
            { ChangeStatus("Удаление задания...","Вычисления будут прерваны");
              action = ActQuit;
              if(proc->WaitFor()) delete proc;
              RunStatus = StNone;
              MainWindow->EditRe->EnableRun-time error at t=0.0535:
You asked to exit. Here you are...
...now exiting to system...
                                                                                                                                                                                                                                                                                                                                                                                                                                  --------------------------------------------------
__fastcall TMainWindow::TMainWindow(TComponent* Owner)
        : TForm(Owner)
{
 ChangeStatus("","Запущена программа, спасибо");

}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::FormCreate(TObject *Sender)
{
   Re=10000;
   maschtab=1000;
   Gamma=1e-3;
   l1=3;                    CheckStep=10;
   l2=1;                    OutStep = 20;
   l3=1;

   nvar=4;
   n1=10;
   n2=5;
   n3=10;
   Ns=25;
   approx=7;                      //derivatives approximation order
   ghost=(approx-1)/2;            //radius of approx sample
   dx[0]=l1/n1;
   dx[1]=l2/n2;
   dx[2]=l3/n3;
   p1 = 8*l1/(l3*Re) ; p2 = 0;


   m1 = n1+2*ghost;
   m2 = n2+2*ghost;
   m3 = n3+2*ghost;
   mm1 = ghost+n1;
   mm2 = ghost+n2;
   mm3 = ghost+n3;

   t_cur=0; Ttot=100;
   count=0; enter = 0;

   NameNuFile = "nut.dat";
   NameVFile  = "vv.dat";
   NameDumpFile = "dump.dat";
   NameEnergyFile = "energy.dat";
   NameKaskadVarFile = "kaskvar.dat";
   NameMessageFile = "message.dat";
   NameErrorFile = "error.err";
   NameInitFile = "init.dat";

   s_func = alloc_mem_2f(n3+2,kol_masht);
   f  =alloc_mem_4f(nvar, m1, m2, m3);   //f[3]-pressure,f[0..2]-v(vector)
   f1 =alloc_mem_4f(nvar, m1, m2, m3);
   df =alloc_mem_4f(nvar, m1, m2, m3);
   df2=alloc_mem_4f(nvar, m1, m2, m3);
   df3=alloc_mem_4f(nvar, m1, m2, m3);
   df4=alloc_mem_4f(nvar, m1, m2, m3);
   df5=alloc_mem_4f(nvar, m1, m2, m3);
   nut=alloc_mem_3f(m1, m2, m3);
   init_shell();

   Noise = 0.1;
   NoiseNorm = 0.05;

   prioritet = tpIdle;
   RunStatus = StNone;
   fileopen("error.err",0);
   RefreshForm();
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::RefreshForm()
{
 AnsiString Temp;
 Temp.printf("%g",Re);           EditRe->Text = Temp;
 Temp.printf("%g",maschtab);     EditMaschtab->Text = Temp;
 Temp.printf("%d",OutStep);      EditOutStep->Text = Temp;
 Temp.printf("%d",CheckStep);    EditCheckStep->Text = Temp;
   EditNu->Text = NameNuFile;
   EditV->Text  = NameVFile;
   EditDump->Text = NameDumpFile;
   EditEnergy->Text = NameEnergyFile;
   EditKaskvar->Text = NameKaskadVarFile;
   EditMessage->Text = NameMessageFile;
   EditError->Text = NameErrorFile;
   EditInit->Text = NameInitFile;
 Temp.printf("%g",Noise);        EditNoise->Text = Temp;
 Temp.printf("%g",NoiseNorm);    EditNorm->Text = Temp;
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::FormDestroy(TObject *Sender)
{
   free_mem_2f(s_func,n3+2,kol_masht);