//---------------------------------------------------------------------------

#ifndef kaskadH
#define kaskadH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <Tabs.hpp>
#include <ComCtrls.hpp>
#include <Tabnotbk.hpp>
#include <Menus.hpp>
#include "CSPIN.h"
//---------------------------------------------------------------------------
class TMainWindow : public TForm
{
__published:	// IDE-managed Components
        TTabbedNotebook *Notebook;
        TButton *BeginWork;
        TGroupBox *GroupOutput;
        TLabel *Label1;
        TEdit *EditMessage;
        TLabel *Label2;
        TEdit *EditEnergy;
        TLabel *Label3;
        TEdit *EditNu;
        TLabel *Label4;
        TEdit *EditKaskvar;
        TLabel *Label5;
        TEdit *EditError;
        TLabel *Label6;
        TEdit *EditDump;
        TLabel *Label7;
        TEdit *EditV;
        TGroupBox *GroupInit;
        TLabel *Label8;
        TEdit *EditInit;
        TButton *InputInit;
        TMemo *Screen;
        TButton *DoDump;
        TButton *Exit;
        TCheckBox *CheckScreenOutput;
        TCheckBox *CheckFileOutput;
        TCheckBox *CheckClear;
        TLabel *Label9;
        TLabel *Label10;
        TLabel *Label11;
        TLabel *Label12;
        TLabel *Label13;
        TLabel *Label14;
        TEdit *EditRe;
        TLabel *Label15;
        TButton *ButtonContinue;
        TGroupBox *GroupInitVariant;
        TRadioButton *RadioInit1;
        TRadioButton *RadioInit2;
        TGroupBox *GroupNoise;
        TCheckBox *CheckNoise;
        TEdit *EditNoise;
        TLabel *Label16;
        TCheckBox *CheckNorm;
        TLabel *Label17;
        TEdit *EditNorm;
        TButton *ClearScr;
        TUpDown *UpDownCheck;
        TEdit *EditCheckStep;
        TEdit *EditOutStep;
        TUpDown *UpDownOut;
        TGroupBox *GroupTurbulent;
        TCheckBox *CheckNut;
        TLabel *Label18;
        TEdit *EditMaschtab;
        void __fastcall FormCreate(TObject *Sender);
        void __fastcall FormDestroy(TObject *Sender);
        void __fastcall EditMessageKeyUp(TObject *Sender, WORD &Key,
          TShiftState Shift);
        void __fastcall EditEnergyKeyUp(TObject *Sender, WORD &Key,
          TShiftState Shift);
        void __fastcall EditVKeyUp(TObject *Sender, WORD &Key,
          TShiftState Shift);
        void __fastcall EditNuKeyUp(TObject *Sender, WORD &Key,
          TShiftState Shift);
        void __fastcall EditKaskvarKeyUp(TObject *Sender, WORD &Key,
          TShiftState Shift);
        void __fastcall EditErrorKeyUp(TObject *Sender, WORD &Key,
          TShiftState Shift);
        void __fastcall EditDumpChange(TObject *Sender);
        void __fastcall EditInitKeyUp(TObject *Sender, WORD &Key,
          TShiftState Shift);
        void __fastcall DoDumpClick(TObject *Sender);
        void __fastcall InputInitClick(TObject *Sender);
        void __fastcall BeginWorkClick(TObject *Sender);
        void __fastcall ExitClick(TObject *Sender);
        void __fastcall EditReKeyUp(TObject *Sender, WORD &Key,
          TShiftState Shift);
        void __fastcall ContinueWorkClick(TObject *Sender);
        void __fastcall EditNoiseKeyUp(TObject *Sender, WORD &Key,
          TShiftState Shift);
        void __fastcall CheckNoiseClick(TObject *Sender);
        void __fastcall CheckNormClick(TObject *Sender);
        void __fastcall EditNormChange(TObject *Sender);
        void __fastcall ClearScrClick(TObject *Sender);
        void __fastcall EditOutStepChanging(TObject *Sender);
        void __fastcall EditCheckStepChanging(TObject *Sender);
        void __fastcall EditMaschtabChange(TObject *Sender);
        void __fastcall CheckNutClick(TObject *Sender);
private:	// User declarations
public:		// User declarations
        __fastcall TMainWindow(TComponent* Owner);
};
  
//---------------------------------------------------------------------------
extern PACKAGE TMainWindow *MainWindow;
//---------------------------------------------------------------------------
#endif
