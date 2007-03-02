// .NAME vtkWFEngineModuleGUI
// .SECTION Description
// GUI for the WFEngine module


#ifndef __vtkWFEngineModuleGUI_h
#define __vtkWFEngineModuleGUI_h

#include "vtkWFEngineModuleWin32Header.h"
#include "vtkSlicerModuleGUI.h"
#include "vtkSlicerModuleLogic.h"
#include "vtkSlicerNodeSelectorWidget.h"
#include "vtkWFEngineModuleLogic.h"

#include <map>
//#include <string>

// Description:
// This class implements Slicer's WFEngineModule GUI
// VTK_WFENGINEMODULE_EXPORT

class vtkKWMultiColumnList;
class vtkKWMyWizardWidget;
class vtkKWPushButtonSet;
class vtkKWWidget;
class vtkMRMLWFEngineModuleNode;
class vtkWFStepHandler;
class vtkSlicerParameterWidget;

//BTX
namespace WFEngine
{
    namespace nmWFStepObject
    {
        class WFStepObject;
    }
}
//ETX

class VTK_WFENGINEMODULE_EXPORT vtkWFEngineModuleGUI : public vtkSlicerModuleGUI
{
 public:
    // Description:
    // Usual vtk class functions
    static vtkWFEngineModuleGUI* New (  );
    vtkTypeRevisionMacro ( vtkWFEngineModuleGUI, vtkSlicerModuleGUI );
    void PrintSelf ( ostream& os, vtkIndent indent );
    
    // Description: Get/Set MRML node
    vtkGetObjectMacro (WFEngineModuleNode, vtkMRMLWFEngineModuleNode);
    virtual void SetWFEngineModuleNode(vtkMRMLWFEngineModuleNode* node);
    
    // Description:
    // Get methods on class members ( no Set methods required. )
    vtkGetObjectMacro ( Logic, vtkWFEngineModuleLogic);
        
    void SetModuleLogic ( vtkWFEngineModuleLogic *logic )
    { this->SetLogic ( vtkObjectPointer (&this->Logic), logic ); }
    void SetAndObserveModuleLogic ( vtkWFEngineModuleLogic *logic )
    { this->SetAndObserveLogic ( vtkObjectPointer (&this->Logic), logic ); }

    // Description:
    // This method builds the WFEngineModule module GUI
    virtual void BuildGUI ( ) ;

    // Description:
    // Add/Remove observers on widgets in the GUI
    virtual void AddGUIObservers ( );
    virtual void RemoveGUIObservers ( );

    // Description:
    // Class's mediator methods for processing events invoked by
    // either the Logic, MRML or GUI.
    virtual void ProcessLogicEvents ( vtkObject *caller, unsigned long event, void *callData );
    virtual void ProcessGUIEvents ( vtkObject *caller, unsigned long event, void *callData );
    virtual void ProcessMRMLEvents ( vtkObject *caller, unsigned long event, void *callData );

    // Description:
    // Methods describe behavior at module enter and exit.
    virtual void Enter ( );
    virtual void Exit ( );
    
    virtual void SaveState();
    
    virtual void workStepGUICallBack();       
    
    // Descritpion:
    // Event that handles all events declared in the workflow xml file
    //BTX
    enum{
        WorkflowHandleEvent = 10000
    };
    
    const char* getStepInputValueByName(std::string name); 
    //ETX
    
protected:
    vtkWFEngineModuleGUI ( );
    virtual ~vtkWFEngineModuleGUI ( );

    // Module logic and mrml pointers
    vtkWFEngineModuleLogic *Logic;
    
    int ConnectToWFEngine();
    virtual void addWorkflowToList(const char* workflowName, const char* fileName, int date);
    
    //Callback Functions
    static void mclDWSelectionChangedCallback(vtkObject* obj, unsigned long,void*, void*);
    static void loadBtnPushCmdCallback(vtkObject* obj, unsigned long,void*, void*);
    static void nextTransitionCallback(vtkObject* obj, unsigned long,void*, void*);
    static void backTransitionCallback(vtkObject* obj, unsigned long,void*, void*);
    static void closeBtnPushCmdCallback(vtkObject* obj, unsigned long,void*, void*);
    
    static void widgetChangedCallback(vtkObject* obj, unsigned long,void*, void*);
    
    virtual void closeWorkflow();
    
    virtual void createWizard();
    
    virtual void UpdateMRML();
    
    virtual void UpdateWorkflowList();
    
    virtual void UpdateParameter();
    
    virtual void UpdateGUI();
//BTX
    virtual void workStepValidationCallBack(WFEngine::nmWFStepObject::WFStepObject* nextWS);
//    virtual void initializeTCLConditions(WFEngine::nmWFStepObject::WFStepObject* curWS);    
//ETX    

    void deleteWizardWidgetContainer();
private:
    int m_curStepID;
    int m_selectedWF;
    
    //BTX
    WFEngine::nmWFStepObject::WFStepObject *m_curWFStep;
    std::map<std::string, std::string> *m_curNameToValueMap;
    //ETX
    
    vtkKWMyWizardWidget *m_curWizWidg;
    vtkKWMultiColumnList *m_mclDW;
    vtkKWPushButtonSet *m_pbtnSet;
    vtkSlicerModuleCollapsibleFrame *m_wizFrame;
    vtkSlicerParameterWidget *m_curParameterWidgets;    
    bool m_inTransition;
    
    vtkMRMLWFEngineModuleNode *WFEngineModuleNode;
    
    vtkWFStepHandler *m_wfStepHandler;
    
    vtkWFEngineModuleGUI ( const vtkWFEngineModuleGUI& ); // Not implemented.
    void operator = ( const vtkWFEngineModuleGUI& ); //Not implemented.
};


#endif
