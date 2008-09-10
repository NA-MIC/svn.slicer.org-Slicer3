/*=auto=========================================================================

  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See Doc/copyright/copyright.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: vtkSlicerTractographyFiducialSeedingGUI.h,v $
  Date:      $Date: 2006/03/19 17:12:29 $
  Version:   $Revision: 1.3 $

=========================================================================auto=*/
#ifndef __vtkSlicerTractographyFiducialSeedingGUI_h
#define __vtkSlicerTractographyFiducialSeedingGUI_h

#include "vtkSlicerBaseGUIWin32Header.h"
#include "vtkSlicerModuleGUI.h"

#include "vtkMRMLScene.h"
#include "vtkMRMLFiducialListNode.h"

#include "vtkSlicerTractographyFiducialSeeding.h"
#include "vtkMRMLTractographyFiducialSeedingNode.h"

class vtkKWFrame;
class vtkSlicerNodeSelectorWidget;
class vtkKWCheckButton;
class vtkKWMenuButtonWithLabel;
class vtkKWScaleWithLabel;

class VTK_FIDUCIALSEEDING_EXPORT vtkSlicerTractographyFiducialSeedingGUI : public vtkSlicerModuleGUI
{
  public:
  static vtkSlicerTractographyFiducialSeedingGUI *New();
  vtkTypeMacro(vtkSlicerTractographyFiducialSeedingGUI, vtkSlicerModuleGUI);
  void PrintSelf(ostream& os, vtkIndent indent);

    // Description: 
    // Get the categorization of the module.
    const char *GetCategory() const
        { return "Tractography.Seeding"; }

  // Description:
  // Events that this module GUI will observe. CLIENT MUST DELETE!
  virtual vtkIntArray* NewObservableEvents();

  // Description:
  // Create widgets
  virtual void BuildGUI ( );

  // Description:
  // Module initialization
  virtual void Init ( );

  // Description:
  // Add obsereves to GUI widgets
  virtual void AddGUIObservers ( );
  
  // Description:
  // Remove obsereves to GUI widgets
  virtual void RemoveGUIObservers ( );
  
  // Description:
  // Pprocess events generated by Logic
  virtual void ProcessLogicEvents ( vtkObject *caller, unsigned long event,
                                  void *callData ){};

  // Description:
  // Pprocess events generated by GUI widgets
  virtual void ProcessGUIEvents ( vtkObject *caller, unsigned long event,
                                  void *callData );

  // Description:
  // Pprocess events generated by MRML
  virtual void ProcessMRMLEvents ( vtkObject *caller, unsigned long event, 
                                  void *callData);
  // Description:
  // Describe behavior at module startup and exit.
  virtual void Enter ( ){};
  virtual void Exit ( ){};

  // Description:
  // Override this method so that nothing will happen when called 
  // by the LoadableModule code.
   virtual void SetModuleLogic ( vtkSlicerLogic* ) { };
  
  // Description:
  // Type of anisotropy used to stop tractography.
  vtkGetStringMacro(StoppingMode);
  vtkSetStringMacro(StoppingMode);
  
  // If StoppingMode criterion becomes smaller than this number,
  // tracking stops.
  vtkGetMacro(StoppingThreshold,vtkFloatingPointType);
  vtkSetMacro(StoppingThreshold,vtkFloatingPointType);

  // Show warning or not
  vtkBooleanMacro(OverwritePolyDataWarning, int);
  vtkGetMacro(OverwritePolyDataWarning, int);
  vtkSetMacro(OverwritePolyDataWarning, int);

   // Description:
  // Set / get the maximum length of the hyperstreamline expressed as absolute
  // distance (i.e., arc length) value.
  vtkSetClampMacro(MaximumPropagationDistance,double,0.0,VTK_DOUBLE_MAX);
  vtkGetMacro(MaximumPropagationDistance,double);

  // Description:
  // Create tracts from fiducuials
  void CreateTracts();

  void SetVolumeSelector(vtkMRMLNode *node);
  void SetFiducialSelector(vtkMRMLNode *node);
  void SetOutFiberSelector(vtkMRMLNode *node);

  // Description: Get/Set MRML node
  vtkGetObjectMacro (TractographyFiducialSeedingNode, vtkMRMLTractographyFiducialSeedingNode);

  void AddFiducialListNodeObserver(vtkMRMLFiducialListNode *n);
  
protected:
  vtkSlicerTractographyFiducialSeedingGUI();
  virtual ~vtkSlicerTractographyFiducialSeedingGUI();
  vtkSlicerTractographyFiducialSeedingGUI(const vtkSlicerTractographyFiducialSeedingGUI&);
  void operator=(const vtkSlicerTractographyFiducialSeedingGUI&);

  // Description:
  // Updates GUI widgets based on parameters values in MRML node
  void UpdateGUI();

  // Description:
  // Updates parameters values in MRML node based on GUI widgets 
  void UpdateMRML();

  void RegisterTractographyFiducialSeedingNode();
  int RegisteredNode;

  char* StoppingMode;
  vtkFloatingPointType StoppingThreshold;
  double MaximumPropagationDistance;
  
  int OverwritePolyDataWarning;

  vtkSlicerNodeSelectorWidget* VolumeSelector;
  vtkSlicerNodeSelectorWidget* FiducialSelector;
  vtkSlicerNodeSelectorWidget* OutFiberSelector;
  
  vtkKWMenuButtonWithLabel *StoppingModeMenu;
  vtkKWScaleWithLabel *StoppingValueScale;
  vtkKWScaleWithLabel *StoppingCurvatureScale;
  vtkKWScaleWithLabel *IntegrationStepLengthScale;

  vtkKWScaleWithLabel *RegionSizeScale;
  vtkKWScaleWithLabel *RegionSampleSizeScale;
  
  vtkKWCheckButton *SeedButton;

  vtkSlicerNodeSelectorWidget* TractographyFiducialSeedingNodeSelector;


  void SetFiducialListNode(vtkMRMLFiducialListNode* FiducialListNode);

  vtkMRMLFiducialListNode* FiducialListNode;
  vtkMRMLTractographyFiducialSeedingNode *TractographyFiducialSeedingNode;

  int UpdatingMRML;
  int UpdatingGUI;
};

#endif

