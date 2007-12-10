#ifndef VTKSLICERVOLUMERENDERINGCUDA_H_
#define VTKSLICERVOLUMERENDERINGCUDA_H_

#include "vtkSlicerModuleGui.h"

class VTK_VOLUMERENDERINGCUDAMODULE_EXPORT vtkVolumeRenderingCudaGUI : vtkSlicerModuleGui
{
 public:
  static vtkVolumeRenderingCudaModuleGUI* New();
  vtkTypeMacro(vtkVolumeRenderingCudaModuleGui, vtkSlicerModuleGUI);


  /// Logic part
  vtkGetObjectMacro(Logic, vtkVolumeRenderingCudaModuleLogic);
  vtkSetObjectMacro(Logic, vtkVolumeRenderingCudaModuleLogic);

  // Description:
  // Process events generated by Logic
  virtual void ProcessLogicEvents ( vtkObject *caller, unsigned long event,
                                    void *callData ){};


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
  
  
  // Description:
  // Methods describe behavior at module enter and exit.
  virtual void Enter ( );
  virtual void Exit ( );
  
  
  void PrintSelf(ostream& os, vtkIndent indent);
  
 protected:
  vtkVolumeRenderingCuda();
  virtual ~vtkVolumeRenderingCuda();

  vtkVolumeRenderingCudaGUI(const vtkVolumeRenderingCudaGUI&); // not implemented
  void operator=(const vtkVolumeRenderingCudaGUI&); // not implemented

  // Description:
  // Pointer to the module's logic class
  vtkVolumeRenderingModuleLogic *Logic;
  vtkVolumeRenderingViewerWidget *ViewerWidget;
};

#endif /*VTKSLICERVOLUMERENDERINGCUDA_H_*/
