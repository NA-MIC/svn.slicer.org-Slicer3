// .NAME vtkSlicerApplication
// .SECTION Description
// Contains slicer's style, application and collection of associated guis.
//


#ifndef __vtkSlicerApplication_h
#define __vtkSlicerApplication_h

#include "vtkSlicerBaseGUIWin32Header.h"
#include "vtkKWApplication.h"
#include "vtkSlicerGUILayout.h"
#include "vtkSlicerGUICollection.h"
#include "vtkSlicerTheme.h"

class vtkSlicerModuleGUI;

// Description:
// Contains slicer's style, application and collection of associated guis.
//
class VTK_SLICER_BASE_GUI_EXPORT vtkSlicerApplication : public vtkKWApplication
{
 public:
    static vtkSlicerApplication* New ( );
    vtkTypeRevisionMacro ( vtkSlicerApplication, vtkKWApplication );

    vtkGetObjectMacro ( MainLayout, vtkSlicerGUILayout );
    vtkGetObjectMacro ( SlicerTheme, vtkSlicerTheme );
    
    vtkGetObjectMacro ( ModuleGUICollection, vtkSlicerGUICollection );

    // Description:
    // This method collects GUIs added to Slicer.
    virtual void AddModuleGUI ( vtkSlicerModuleGUI *gui );
    virtual vtkSlicerModuleGUI* GetModuleGUIByName ( const char *name );

    // Description:
    // These methods manage windows associated with the application
    virtual void CloseAllWindows ( ) ;
    
    // Description:
    // Sets application behavior.
    virtual void ConfigureApplication ( );
    virtual int StartApplication ( );
    // Description:
    // installs rules to specify look & feel.
    virtual void InstallTheme ( vtkKWTheme *theme );


 protected:
    vtkSlicerApplication ( );
    ~vtkSlicerApplication ( );

    vtkSlicerGUILayout *MainLayout;
    vtkSlicerTheme *SlicerTheme;
    
    // Description:
    // Collections of GUIs
    vtkSlicerGUICollection *ModuleGUICollection;
    
 private:
    vtkSlicerApplication ( const vtkSlicerApplication& ); // Not implemented.
    void operator = ( const vtkSlicerApplication& ); //Not implemented.
}; 

#endif
