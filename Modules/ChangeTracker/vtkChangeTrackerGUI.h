#ifndef __vtkChangeTrackerGUI_h
#define __vtkChangeTrackerGUI_h

#include "vtkSlicerModuleGUI.h"
#include "vtkChangeTracker.h"
#include "vtkSlicerSliceLogic.h"
#include "vtkKWScale.h"

class vtkChangeTrackerLogic;
class vtkMRMLChangeTrackerNode;
class vtkKWWizardWidget;
class vtkChangeTrackerFirstScanStep;
class vtkChangeTrackerROIStep;
class vtkChangeTrackerSegmentationStep;
class vtkChangeTrackerTypeStep;
class vtkChangeTrackerAnalysisStep;
class vtkImageData;
class vtkMRMLROINode;

class VTK_CHANGETRACKER_EXPORT vtkChangeTrackerGUI : 
  public vtkSlicerModuleGUI
{
public:
  static vtkChangeTrackerGUI *New();
  vtkTypeMacro(vtkChangeTrackerGUI,vtkSlicerModuleGUI);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description: Get the categorization of the module.
  const char *GetCategory() const { return "Wizards"; }

  // Description: 
  // Get/Set logic node
  vtkGetObjectMacro(Logic, vtkChangeTrackerLogic);
  virtual void SetLogic(vtkChangeTrackerLogic*);
  

  // Description:
  // Set the logic pointer from parent class pointer.
  // Overloads implementation in vtkSlicerModulesGUI
  // to allow loadable modules.
  virtual void SetModuleLogic ( vtkSlicerLogic *logic )
  {
    this->SetLogic(reinterpret_cast<vtkChangeTrackerLogic*> (logic)); 
  }

  // Description: 
  // Get/Set MRML node
  vtkGetObjectMacro(Node, vtkMRMLChangeTrackerNode);
  virtual void SetNode(vtkMRMLChangeTrackerNode*);

  // Description: 
  // Get wizard widget
  vtkGetObjectMacro(WizardWidget, vtkKWWizardWidget);
  // vtkGetObjectMacro(AnatomicalStructureStep, vtkChangeTrackerAnatomicalStructureStep);

  // Description:
  // Events that this module GUI will observe. CLIENT MUST DELETE!
  virtual vtkIntArray* NewObservableEvents();

  // Description:
  // Create widgets
  virtual void BuildGUI();

  // Description:
  // Delete Widgets
  virtual void TearDownGUI();

  // Description:
  // Add observers to GUI widgets
  virtual void AddGUIObservers();
  
  // Description:
  // Remove observers to GUI widgets
  virtual void RemoveGUIObservers();

  // Description:
  // Remove observers to MRML node
  virtual void RemoveMRMLNodeObservers();

  // Description:
  // Remove observers to Logic
  virtual void RemoveLogicObservers();
  
  // Description:
  // Pprocess events generated by Logic
  virtual void ProcessLogicEvents( vtkObject *caller, unsigned long event,
                                   void *callData);

  // Description:
  // Pprocess events generated by GUI widgets
  virtual void ProcessGUIEvents( vtkObject *caller, unsigned long event,
                                 void *callData);

  // Description:
  // Process events generated by MRML
  // This function is automatically called 
  virtual void ProcessMRMLEvents( vtkObject *caller, unsigned long event, 
                                  void *callData);
  // Description:
  // Describe behavior at module startup and exit.
  virtual void Enter();
  virtual void Exit();

  // Description: The name of the Module - this is used to 
  // construct the proc invocations
  vtkGetStringMacro(ModuleName);
  vtkSetStringMacro(ModuleName);

  bool GetModuleEntered() { return ModuleEntered;};

  // Description: set an observer by number (work around
  // limitation in kwwidgets tcl wrapping)
  unsigned long AddObserverByNumber(vtkObject *observee, unsigned long event);

  // Description:
  // Updates parameters values in MRML node based on GUI widgets 
  void UpdateMRML();
  
  vtkChangeTrackerFirstScanStep* GetFirstScanStep() {return this->FirstScanStep;}
  void UpdateNode();

  vtkGetObjectMacro(SegmentationStep,vtkChangeTrackerSegmentationStep);

  void SliceLogicRemove();
  void SliceLogicDefine();

  // Description:
  // accessor
  vtkGetObjectMacro(SliceLogic, vtkSlicerSliceLogic);
  vtkGetObjectMacro(SliceController_OffsetScale, vtkKWScale); 

  void PropagateVolumeSelection();

  void ResetPipeline();

  void ObserveMRMLROINode(vtkMRMLROINode* roi);

protected:
   static void SliceLogicCallback(vtkObject *caller, unsigned long event, void *clientData, void *callData );
private:
  vtkChangeTrackerGUI();
  ~vtkChangeTrackerGUI();
  vtkChangeTrackerGUI(const vtkChangeTrackerGUI&);
  void operator=(const vtkChangeTrackerGUI&);

  // Description:
  // Updates GUI widgets based on parameters values in MRML node
  void UpdateGUI();

  // Description:
  // Updates registration progress on the status bar of the main application. 
  virtual void UpdateRegistrationProgress();

  void SliceLogicRemoveGUIObserver();

  vtkChangeTrackerLogic       *Logic;
  vtkMRMLChangeTrackerNode    *Node;
  
  char *ModuleName;

  // Description:
  // The wizard widget and steps
  vtkKWWizardWidget                      *WizardWidget;
  vtkChangeTrackerFirstScanStep        *FirstScanStep;
  vtkChangeTrackerROIStep              *ROIStep;
  vtkChangeTrackerSegmentationStep     *SegmentationStep;
  vtkChangeTrackerTypeStep             *TypeStep;
  vtkChangeTrackerAnalysisStep         *AnalysisStep;

  vtkSlicerSliceLogic *SliceLogic;
  vtkKWScale *SliceController_OffsetScale;
  vtkCallbackCommand *SliceLogicCallbackCommand;

  // Wizard step cannot observe MRML events
  vtkMRMLROINode *roiNode;

  bool ModuleEntered;
};

#endif
