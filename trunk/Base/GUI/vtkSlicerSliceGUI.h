// .NAME vtkSlicerSliceGUI 
// .SECTION Description
// Individual Slice GUI and mediator functions for slicer3.
// Contains a SliceViewer, a Slice Controller, a pointer to
// SliceLogic and a pointer to a MRMLSliceNode.


#ifndef __vtkSlicerSliceGUI_h
#define __vtkSlicerSliceGUI_h

#include "vtkSlicerBaseGUIWin32Header.h"
#include "vtkSlicerComponentGUI.h"
#include "vtkSlicerSliceViewer.h"
#include "vtkSlicerSliceController.h"
#include "vtkSlicerSliceLogic.h"
#include "vtkMRMLSliceNode.h"

class vtkObject;
class vtkKWFrame;

// Description:
// This class implements Slicer's Slice GUI.
//
class VTK_SLICER_BASE_GUI_EXPORT vtkSlicerSliceGUI : public vtkSlicerComponentGUI
{
 public:
    static vtkSlicerSliceGUI* New (  );
    vtkTypeRevisionMacro ( vtkSlicerSliceGUI, vtkSlicerComponentGUI );
    void PrintSelf (ostream& os, vtkIndent indent);

    // Description:
    // Get methods for three default SliceGUIs
    // (Each SliceGUI contains a SliceViewerWidget,
    // SliceControllerWidget, a SliceLogic pointer and
    // a SliceNode pointer.)
    vtkGetObjectMacro ( SliceViewer, vtkSlicerSliceViewer );
    vtkGetObjectMacro ( SliceController, vtkSlicerSliceController );
    vtkGetMacro ( ControllerStyle, int );
    vtkSetMacro ( ControllerStyle, int );
    vtkGetObjectMacro ( Logic, vtkSlicerSliceLogic );
    vtkGetObjectMacro ( SliceNode, vtkMRMLSliceNode );
    
    // Description:
    // API for setting SliceNode, SliceLogic and
    // for both setting and observing them.
    void SetMRMLNode ( vtkMRMLSliceNode *node )
        { this->SetMRML ( vtkObjectPointer( &this->SliceNode), node ); }
    void SetAndObserveMRMLNode ( vtkMRMLSliceNode *node )
        { this->SetMRML ( vtkObjectPointer( &this->SliceNode), node ); }
    void SetModuleLogic ( vtkSlicerSliceLogic *logic )
        { this->SetLogic ( vtkObjectPointer (&this->Logic), logic ); }
    void SetAndObserveModuleLogic ( vtkSlicerSliceLogic *logic )
        { this->SetLogic ( vtkObjectPointer (&this->Logic), logic ); }

    // Description:
    // Build the SlicesGUI's UIPanel and three main SliceGUIs 
    virtual void BuildGUI ( vtkKWFrame *f );

    // Description:
    // Add/Remove Observers on UIPanel widgets and SliceGUIs.
    virtual void AddGUIObservers ( );
    virtual void RemoveGUIObservers ( );
    
    // Description:
    // Processes all events raised by the logic
    virtual void ProcessLogicEvents ( vtkObject *caller, unsigned long event, void *callData );
    // Description:
    // Processes all events raised by the GUI
    virtual void ProcessGUIEvents ( vtkObject *caller, unsigned long event, void *callData );
    // Description:
    // Processes all events raised by MRML
    virtual void ProcessMRMLEvents ( vtkObject *caller, unsigned long event, void *callData );
    
    // Description:
    // Functions that define and undefine module-specific behaviors.
    virtual void Enter ( );
    virtual void Exit ( );
    
 protected:
    vtkSlicerSliceGUI ( );
    ~vtkSlicerSliceGUI ( );

    // Description:
    // Three slice widgets by default.
    vtkSlicerSliceViewer *SliceViewer;
    vtkSlicerSliceController *SliceController;
    int ControllerStyle;
    vtkSlicerSliceLogic *Logic;
    vtkMRMLSliceNode *SliceNode;

 private:
    vtkSlicerSliceGUI ( const vtkSlicerSliceGUI& ); // Not implemented.
    void operator = ( const vtkSlicerSliceGUI& ); //Not implemented.
}; 

#endif
