#ifndef __vtkTumorGrowthGUI_h
#define __vtkTumorGrowthGUI_h

#include "vtkSlicerModuleGUI.h"
#include "vtkTumorGrowth.h"

class vtkTumorGrowthLogic;
class vtkMRMLTumorGrowthNode;
class vtkKWWizardWidget;
class vtkTumorGrowthFirstScanStep;
class vtkTumorGrowthROIStep;
class vtkTumorGrowthSegmentationStep;
class vtkTumorGrowthSecondScanStep;
class vtkTumorGrowthAnalysisStep;

class VTK_TUMORGROWTH_EXPORT vtkTumorGrowthGUI : 
  public vtkSlicerModuleGUI
{
public:
  static vtkTumorGrowthGUI *New();
  vtkTypeMacro(vtkTumorGrowthGUI,vtkSlicerModuleGUI);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description: Get the categorization of the module.
  //   const char *GetCategory() const
  //     { return "Segmentation"; }

  // Description: 
  // Get/Set logic node
  vtkGetObjectMacro(Logic, vtkTumorGrowthLogic);
  virtual void SetLogic(vtkTumorGrowthLogic*);
  
  // Description: 
  // Get/Set MRML node
  vtkGetObjectMacro(Node, vtkMRMLTumorGrowthNode);
  virtual void SetNode(vtkMRMLTumorGrowthNode*);

  // Description: 
  // Get wizard widget
  vtkGetObjectMacro(WizardWidget, vtkKWWizardWidget);
  // vtkGetObjectMacro(AnatomicalStructureStep, vtkTumorGrowthAnatomicalStructureStep);

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
  virtual void Enter(){};
  virtual void Exit(){};

  // Description: The name of the Module - this is used to 
  // construct the proc invocations
  vtkGetStringMacro(ModuleName);
  vtkSetStringMacro(ModuleName);

  // Description: set an observer by number (work around
  // limitation in kwwidgets tcl wrapping)
  unsigned long AddObserverByNumber(vtkObject *observee, unsigned long event);

  // Description:
  // Updates parameters values in MRML node based on GUI widgets 
  void UpdateMRML();
  
  vtkTumorGrowthFirstScanStep* GetFirstScanStep() {return this->FirstScanStep;}

protected:

private:
  vtkTumorGrowthGUI();
  ~vtkTumorGrowthGUI();
  vtkTumorGrowthGUI(const vtkTumorGrowthGUI&);
  void operator=(const vtkTumorGrowthGUI&);

  // Description:
  // Updates GUI widgets based on parameters values in MRML node
  void UpdateGUI();

  // Description:
  // Updates registration progress on the status bar of the main application. 
  virtual void UpdateRegistrationProgress();

  vtkTumorGrowthLogic       *Logic;
  vtkMRMLTumorGrowthNode    *Node;
  
  char *ModuleName;

  // Description:
  // The wizard widget and steps
  vtkKWWizardWidget                      *WizardWidget;
  vtkTumorGrowthFirstScanStep        *FirstScanStep;
  vtkTumorGrowthROIStep              *ROIStep;
  vtkTumorGrowthSegmentationStep     *SegmentationStep;
  vtkTumorGrowthSecondScanStep       *SecondScanStep;
  vtkTumorGrowthAnalysisStep         *AnalysisStep;

  // Description:
  // Populate the logic with testing data, load some volumes
  virtual void PopulateTestingData();
};

#endif
