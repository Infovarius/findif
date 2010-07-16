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
   n3=10;
   m1 = n1+2*ghost;
   mm1 = ghost+n1;
   m2 = n2+2*ghost;
   mm2 = ghost+n2;
   m3 = n3+2*ghost;
   mm3 = ghost+n3;

   Ns=25;
   approx=7;                      //derivatives approximation order
   ghost=(approx-1)/2;            //radius of approx sample
   dx[0]=l1/n1;
   dx[1]=l2/n2;
   dx[2]=l3/n3;
   p1 = 8*l1/(l3*Re) ; p2 = 0;

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

//   Noise = 0.1;
//   NoiseNorm = 0.0;

   prioritet = tpIdle;
   RunStatus = StNone;
   fileopen("error.err",0);
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
  if ( (code=init_data(&n1,&n2,&n3,&Re,f,nut,&t_cur,&count)) > 0)
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
              MainWindow->EditRe->Enabled = true;
              MainWindow->EditRe->ShowHint = false;
            }
     else MainWindow->Close();
}
//---------------------------------------------------------------------------

void __fastcall TMainWindow::EditCheckStepChanging(TObject *Sender)
{
 CheckStep = MainWindow->EditCheckStep->Text.ToInt();
}
//---------------------------------------------------------------------------

void __fastcall TMainWindow::EditOutStepChanging(TObject *Sender)
{
 OutStep = MainWindow->EditOutStep->Text.ToInt();
}
//---------------------------------------------------------------------------

void __fastcall TMainWindow::EditReKeyUp(TObject *Sender, WORD &Key,
      TShiftState Shift)
{
 Re = MainWindow->EditRe->Text.ToDouble();
 if(Re<1) Re=1.;
 p1 = 8*l1/(l3*Re);
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::EditNoiseKeyUp(TObject *Sender, WORD &Key,
      TShiftState Shift)
{
 Noise = MainWindow->EditNoise->Text.ToDouble();
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::CheckNoiseClick(TObject *Sender)
{
 MainWindow->EditNoise->Enabled = CheckNoise->Checked;
}
//---------------------------------------------------------------------------

void __fastcall TMainWindow::CheckNormClick(TObject *Sender)
{
 MainWindow->EditNorm->Enabled = CheckNorm->Checked;
}
//---------------------------------------------------------------------------

void __fastcall TMainWindow::EditNormChange(TObject *Sender)
{
 NoiseNorm = MainWindow->EditNorm->Text.ToDouble();
}
//---------------------------------------------------------------------------

void __fastcall TMainWindow::ClearScrClick(TObject *Sender)
{
 MainWindow->Screen->Lines->Clear();         
}
//---------------------------------------------------------------------------

void __fastcall TMainWindow::EditMaschtabChange(TObject *Sender)
{
  maschtab = MainWindow->EditMaschtab->Text.ToDouble();        
}
//---------------------------------------------------------------------------

void __fastcall TMainWindow::CheckNutClick(TObject *Sender)
{
 MainWindow->EditMaschtab->Enabled = CheckNut->Checked;
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::SelPriorChange(TObject *Sender)
{
switch (MainWindow->SelPrior->Text[1])
 {
   case '0' : prioritet = tpIdle;   break;
   case '1' : prioritet = tpLowest;   break;
   case '2' : prioritet = tpLower;   break;
   case '3' : prioritet = tpNormal;   break;
   case '4' : prioritet = tpHigher;   break;
   case '5' : prioritet = tpHighest;   break;
 }
 if(proc) proc->Priority = prioritet;
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::ChangeStatus(char *str1,char *str2)
{
  StatusBarCalc->Panels->BeginUpdate();
  try
  {
   if(str1[0]!='\0') MainWindow->StatusBarCalc->Panels->Items[0]->Text=str1;
   if(str2[0]!='\0') MainWindow->StatusBarCalc->Panels->Items[1]->Text=str2;
  }
  __finally
  {
    StatusBarCalc->Panels->EndUpdate();
  }
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::EditNxChange(TObject *Sender)
{
 n1 = MainWindow->EditNx->Text.ToIntDef(10);
 m1 = n1+2*ghost;
 mm1 = ghost+n1;
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::EditNyChange(TObject *Sender)
{
 n2 = MainWindow->EditNy->Text.ToIntDef(10);
 m2 = n2+2*ghost;
 mm2 = ghost+n2;
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::EditNzChange(TObject *Sender)
{
 n3 = MainWindow->EditNz->Text.ToIntDef(20);
 m3 = n3+2*ghost;
 mm3 = ghost+n3;
}
//---------------------------------------------------------------------------

