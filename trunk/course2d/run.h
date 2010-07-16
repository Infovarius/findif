//---------------------------------------------------------------------------

#ifndef runH
#define runH
//---------------------------------------------------------------------------
#include <Classes.hpp>
//---------------------------------------------------------------------------
class CalcProcess : public TThread
{
private:
protected:
        void __fastcall Execute();
public:
        __fastcall CalcProcess(bool CreateSuspended);
};
//---------------------------------------------------------------------------
#endif
