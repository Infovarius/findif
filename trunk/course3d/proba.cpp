//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "proba.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm1 *Form1;
//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Button1Click(TObject *Sender)
{
TMyThread *SecondProcess = new TMyThread(true); // create but don’t run
SecondProcess->Priority = tpHigher; // set the priority lower than normal
SecondProcess->Resume(); // now start the thread running
while(a>0){
 SecondProcess-
 }
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Button2KeyPress(TObject *Sender, char &Key)
{
 a=0;
}
//---------------------------------------------------------------------------
