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
//                            initial actions and procedures
//---------------------------------------------------------------------------
__fastcall TMainWindow::TMainWindow(TComponent* Owner)
        : TForm(Owner)
{
 ChangeStatus("","Запущена программа, спасибо");

}
//---------------------------------------------------------------------------
void allocate_memory()
{
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
}
//---------------------------------------------------------------------------
void free_memory(double **s_func,double ****f,double ****f1,double ****df,double ****df2,
                     double ****df3,double ****df4,double ****df5,double ***nut)
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
void __fastcall TMainWindow::FormCreate(TObject *Sender)
{
   #include "init_vars.h"

   Re=10000;
   maschtab=1000;
   Gamma=1e-3;
   l1=1;                    CheckStep=10;
   l2=1;                    OutStep = 20;
   l3=1;

   nvar=4;
   n1=5;
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

   t_cur=0; Ttot=1000;
   count=0; enter = 0;

   num_dump = 0;
   allocate_memory();
   Noise = 0.1;
   NoiseNorm = 0.05;

   prioritet = tpIdle;
   RunStatus = StNone;
//   fileopen("error.err",0);
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
 Temp.printf("%d",n1);           EditNx->Text = Temp;
 Temp.printf("%d",n2);           EditNy->Text = Temp;
 Temp.printf("%d",n3);           EditNz->Text = Temp;
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::FormDestroy(TObject *Sender)
{
   free_memory(s_func,f,f1,df,df2,df3,df4,df5,nut);
}
//---------------------------------------------------------------------------
//                         Main Block
void __fastcall TMainWindow::BeginWorkClick(TObject *Sender)
{
 switch (RunStatus)
 {
   case StNone:  { RunStatus = StRun;
                   ChangeStatus("Загрузка задания...","Сейчас будут начаты вычисления");
                   EditRe->Enabled = false;
                   EditRe->Hint = "Число не может быть изменено во время счета";
                   EditRe->ShowHint = false;
                   EditNx->Enabled = false;
                   EditNy->Enabled = false;
                   EditNz->Enabled = false;
                   action = ActBegin;
                   if(!input_done)
                     { init_conditions(f,Re);
                       t_cur = 0;  count = 0;  enter = 0;
                     }  
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
//---------------------------------------------------------------------------
void __fastcall TMainWindow::RecontinueClick(TObject *Sender)
{
   ChangeStatus("Перезагрузка задания...","Вычисления будут начаты сначала");
   RunStatus = StRun;
   action = ActBegin;
   init_conditions(f,Re);
   delete proc;
   t_cur = 0;  count = 0;  enter = 0;
   proc = new CalcProcess(false);
   proc->Priority = prioritet;
   ChangeStatus("Идут вычисления...","Вычисления начаты сначала");
}
//---------------------------------------------------------------------------
void at_the_end_of_work()
{
     MainWindow->EditRe->Enabled = true;
     MainWindow->EditRe->ShowHint = false;
     MainWindow->EditNx->Enabled = true;
     MainWindow->EditNy->Enabled = true;
     MainWindow->EditNz->Enabled = true;
  RunStatus = StNone;
  action = ActQuit;
  MainWindow->BeginWork->Caption = "Начать счет";
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::ExitClick(TObject *Sender)
{
 if(RunStatus == StRun || RunStatus == StPause)
            { ChangeStatus("Удаление задания...","Вычисления будут прерваны");
              action = ActQuit;
              if(proc->WaitFor()) delete proc;
              RunStatus = StNone;
              at_the_end_of_work();
            }
     else MainWindow->Close();
}
//---------------------------------------------------------------------------
//                           other handlers
void __fastcall TMainWindow::EditMessageKeyUp(TObject *Sender, WORD &Key,
      TShiftState Shift)
{
  NameMessageFile = EditMessage->Text;
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::EditEnergyKeyUp(TObject *Sender, WORD &Key,
      TShiftState Shift)
{
  NameEnergyFile = EditEnergy->Text;
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::EditVKeyUp(TObject *Sender, WORD &Key,
      TShiftState Shift)
{
  NameVFile = EditV->Text;
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::EditNuKeyUp(TObject *Sender, WORD &Key,
      TShiftState Shift)
{
  NameNuFile = EditNu->Text;
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::EditKaskvarKeyUp(TObject *Sender, WORD &Key,
      TShiftState Shift)
{
  NameKaskadVarFile = EditKaskvar->Text;
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::EditErrorKeyUp(TObject *Sender, WORD &Key,
      TShiftState Shift)
{
  NameErrorFile = EditError->Text;
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::EditDumpChange(TObject *Sender)
{
  NameDumpFile = EditDump->Text;
  num_dump = 0;
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::EditInitKeyUp(TObject *Sender, WORD &Key,
      TShiftState Shift)
{
  NameInitFile = EditInit->Text;
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::InputInitClick(TObject *Sender)
{
  int code=0;
  long dummy_count;
  double dummy_Re,dummy_t;
  if (RadioInitVel->Checked) //from beginning
         { code=init_data(&dummy_Re,&dummy_t,&dummy_count);
           t_cur = 0;  count = 0;  enter = 0; }
         else code=init_data(&Re,&t_cur,&count);
  RefreshForm();
  if ( code > 0)
     nrerror("Не удался ввод начальных данных. Код ошибки = "/*,code*/,0);
   else input_done = true;
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::DoDumpClick(TObject *Sender)
{
  dump(n1,n2,n3,Re,f,nut,t_cur,count);
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::EditCheckStepChanging(TObject *Sender)
{
 CheckStep = MainWindow->EditCheckStep->Text.ToIntDef(10);
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::EditOutStepChanging(TObject *Sender)
{
 OutStep = MainWindow->EditOutStep->Text.ToIntDef(20);
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::EditReKeyUp(TObject *Sender, WORD &Key,
      TShiftState Shift)
{
 Re = atof(MainWindow->EditRe->Text.c_str());
 if(Re<1) Re=1.;
 p1 = 8*l1/(l3*Re);
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::EditNoiseKeyUp(TObject *Sender, WORD &Key,
      TShiftState Shift)
{
 Noise = atof(MainWindow->EditNoise->Text.c_str());
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
 NoiseNorm = atof(MainWindow->EditNorm->Text.c_str());
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::ClearScrClick(TObject *Sender)
{
 MainWindow->Screen->Lines->Clear();
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::EditMaschtabChange(TObject *Sender)
{
  maschtab = atof(MainWindow->EditMaschtab->Text.c_str());
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::CheckNutClick(TObject *Sender)
{
// MainWindow->EditMaschtab->Enabled = CheckNut->Checked;
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
void __fastcall TMainWindow::ChangeStatus(AnsiString str1,AnsiString str2)
{
  StatusBarCalc->Panels->BeginUpdate();
  try
  {
   if(!str1.IsEmpty()) MainWindow->StatusBarCalc->Panels->Items[0]->Text=str1;
   if(!str2.IsEmpty()) MainWindow->StatusBarCalc->Panels->Items[1]->Text=str2;
  }
  __finally
  {
    StatusBarCalc->Panels->EndUpdate();
  }
}
//---------------------------------------------------------------------------
void recreate_arrays(int n1new,int n2new,int n3new)
{
  free_memory(s_func,f,f1,df,df2,df3,df4,df5,nut);
  n1 = n1new;
  n2 = n2new;
  n3 = n3new;
  m1 = n1+2*ghost;
  m2 = n2+2*ghost;
  m3 = n3+2*ghost;
  mm1 = n1+ghost;
  mm2 = n2+ghost;
  mm3 = n3+ghost;
  dx[0]=l1/n1;
  dx[1]=l2/n2;
  dx[2]=l3/n3;
  allocate_memory();
}

void __fastcall TMainWindow::ButtonSaveNiClick(TObject *Sender)
{
int n1,n2,n3;
  n1 = EditNx->Text.ToIntDef(5);
  if(n1<5) n1 = 5;
  n2 = EditNy->Text.ToIntDef(5);
  if(n2<5) n2 = 5;
  n3 = EditNz->Text.ToIntDef(5);
  if(n3<5) n3 = 5;
  recreate_arrays(n1,n2,n3);
}
//---------------------------------------------------------------------------

