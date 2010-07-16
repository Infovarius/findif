//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop
USERES("course.res");
USEFORM("kaskad.cpp", MainWindow);
USEUNIT("ioutil.cpp");
USEUNIT("UTIL.CPP");
USEUNIT("GOY.CPP");
USEUNIT("pde.cpp");
USEUNIT("RKAD.CPP");
USEUNIT("struct_func.cpp");
USE("course.todo", ToDo);
USEUNIT("run.cpp");
//---------------------------------------------------------------------------
WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int)
{
        try
        {
                 Application->Initialize();
                 Application->Title = "Каскадно-сеточный метод";
                 Application->CreateForm(__classid(TMainWindow), &MainWindow);
                 Application->Run();
        }
        catch (Exception &exception)
        {
                 Application->ShowException(&exception);
        }
        return 0;
}
//---------------------------------------------------------------------------
