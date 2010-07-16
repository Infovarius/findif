//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#define LEVEL
#include "kaskad.h"
#include "head.h"
#include <time.h>
#include <conio.h>
#include <system.hpp>
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma link "CSPIN"
#pragma resource "*.dfm"
TMainWindow *MainWindow;
//---------------------------------------------------------------------------
__fastcall TMainWindow::TMainWindow(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::FormCreate(TObject *Sender)
{
   Re=10000;
   Gamma=1e-3;
   l1=3;                    CheckStep=10;
   l2=1;                    OutStep = 20;
   l3=1;

   nvar=4;
   n1=10;
   n2=10;
   n3=20;
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
}
//---------------------------------------------------------------------------
// ------------------------Main Block-----------------------------------------
void __fastcall TMainWindow::BeginWorkClick(TObject *Sender)
{
   double dttry, dtdid, dtnext, PulsEn;
   int i,j,k,l;
   BeginWork->Caption="Пауза";
   time(&time_begin);
   t_cur = 0;
   count = 0;  enter = 0;
   nmessage("work has begun",t_cur);
   init_conditions(f,Re);
   PulsEn=check(f);
   printing(f,dtdid,t_cur,count,PulsEn);

   dtnext=1e-3;
   dump(n1,n2,n3,Re,f,nut,t_cur,count);

        /*----------- MAIN ITERATIONS --------------*/
   while (t_cur < Ttot && !razlet) {
        pde(t_cur, f, df);
        dttry=dtnext;
        timestep(f, df, t_cur, f1, dttry, &dtdid, &dtnext);
        nut_by_flux(f,dtdid);
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
    /*    if(kbhit())
             {
                switch (getch()) {
                        case 'd' : dump(n1,n2,n3,Re,f,nut,t_cur,count); break;
                        case 'q' : { dump(n1,n2,n3,Re,f,nut,t_cur,count);
                                     nrerror("You asked to exit. Here you are...",t_cur);
                                   }
                        case 'p' : return;
                        }
             }*/
   }
  BeginWork->Caption="Начать счет";
  if(t_cur>Ttot&&!razlet) nmessage("work is succesfully done",t_cur);
       else nrerror("this is break of scheme",t_cur);
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::ContinueWorkClick(TObject *Sender)
{
   double dttry, dtdid, dtnext, PulsEn;
   int i,j,k,l;
   time(&time_begin);
   nmessage("work has begun",t_cur);
   init_conditions(f,Re);
   PulsEn=check(f);
   printing(f,dtdid,t_cur,count,PulsEn);

   dtnext=1e-3;
   dump(n1,n2,n3,Re,f,nut,t_cur,count);

        /*----------- MAIN ITERATIONS --------------*/
   while (t_cur < Ttot && !razlet) {
        pde(t_cur, f, df);
        dttry=dtnext;
        timestep(f, df, t_cur, f1, dttry, &dtdid, &dtnext);
        nut_by_flux(f,dtdid);
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
    /*    if(kbhit())
             {
                switch (getch()) {
                        case 'd' : dump(n1,n2,n3,Re,f,nut,t_cur,count); break;
                        case 'q' : { dump(n1,n2,n3,Re,f,nut,t_cur,count);
                                     nrerror("You asked to exit. Here you are...",t_cur);
                                   }
                        case 'p' : return;
                        }
             }*/
   }
  if(t_cur>Ttot&&!razlet) nmessage("work is succesfully done",t_cur);
       else nrerror("this is break of scheme",t_cur);
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
}
//---------------------------------------------------------------------------
void __fastcall TMainWindow::DoDumpClick(TObject *Sender)
{
  dump(n1,n2,n3,Re,f,nut,t_cur,count);
}
//---------------------------------------------------------------------------

void __fastcall TMainWindow::ExitClick(TObject *Sender)
{
 MainWindow->Close();
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

