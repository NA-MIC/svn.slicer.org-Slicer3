#ifndef __ULTRASOUND_EXAMPLE_GUI_H_
#define __ULTRASOUND_EXAMPLE_GUI_H_

#include "vtkSlicerApplicationGUI.h"
#include "vtkUltrasoundModule.h"

class vtkUltrasoundModuleLogic;
class vtkKWCheckButton;


class VTK_ULTRASOUNDMODULE_EXPORT vtkUltrasoundModuleGUI : public vtkSlicerApplicationGUI
{
    vtkTypeRevisionMacro(vtkUltrasoundModuleGUI, vtkSlicerApplicationGUI);
    static vtkUltrasoundModuleGUI *New();


protected:
    vtkUltrasoundModuleGUI();
    virtual ~vtkUltrasoundModuleGUI();

    // Description:
    // Process events generated by Logic
    virtual void ProcessLogicEvents ( vtkObject *caller, unsigned long event, void *callData );

    /// GUI part
    virtual void BuildGUI ( );
    // This method releases references and key-bindings,
    // and optionally removes observers.
    virtual void TearDownGUI ( );

    // Description:
    // Methods for adding module-specific key bindings and
    // removing them.
    virtual void CreateModuleEventBindings ( );
    virtual void ReleaseModuleEventBindings ( );

    // Description:
    // Add obsereves to GUI widgets
    virtual void AddGUIObservers ( );

    // Description:
    // Remove obsereves to GUI widgets
    virtual void RemoveGUIObservers ( );
    virtual void RemoveMRMLNodeObservers ( );
    virtual void RemoveLogicObservers ( );


    // Description:
    // Process events generated by GUI widgets
    virtual void ProcessGUIEvents ( vtkObject *caller, unsigned long event,
        void *callData );

    // Description:
    // Process events generated by MRML
    virtual void ProcessMRMLEvents ( vtkObject *caller, unsigned long event,
        void *callData);

    /// Logic part
    vtkGetObjectMacro(Logic, vtkUltrasoundModuleLogic);
    // this does not work vtkSetObjectMacro(Logic, vtkVolumeRenderingCudaModuleLogic);
    virtual void SetLogic(vtkUltrasoundModuleLogic *logic)
    {
        this->Logic=logic;
    }

    virtual void CreateWidget();

    // Description:
    // Methods describe behavior at module enter and exit.
    virtual void Enter ( );
    virtual void Exit ( );


    void PrintSelf(ostream& os, vtkIndent indent);

protected:
        vtkMRMLVolumeNode* AddVolumeNode(const char* volumeNodeName);


private:
    vtkUltrasoundModuleGUI operator=(const vtkUltrasoundModuleGUI&);
    vtkUltrasoundModuleGUI(const vtkUltrasoundModuleGUI&);

private:
    vtkKWCheckButton*           cb_Scanning;

    vtkUltrasoundModuleLogic*   Logic;
};

#endif /* __ULTRASOUND_EXAMPLE_GUI_H_ */
