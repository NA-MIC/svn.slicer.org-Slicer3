/*=auto=========================================================================

Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

See Doc/copyright/copyright.txt
or http://www.slicer.org/copyright/copyright.txt for details.

Program:   3D Slicer
Module:    $RCSfile: vtkVolumeRenderingGUI.h,v $
Date:      $Date: 2006/03/19 17:12:29 $
Version:   $Revision: 1.3 $

=========================================================================auto=*/
#ifndef __vtkVolumeRenderingGUI_h
#define __vtkVolumeRenderingGUI_h

#include "vtkSlicerModuleGUI.h"
#include "vtkVolumeRendering.h"
#include "vtkVolumeRenderingLogic.h"

#include "vtkMRMLVolumeRenderingNode.h"
#include "vtkSlicerNodeSelectorWidget.h"
#include "vtkSlicerNodeSelectorVolumeRenderingWidget.h"

#include "vtkKWLabel.h"
#include "vtkKWHistogram.h"
#include "vtkKWEntryWithLabel.h"
#include "vtkKWTkUtilities.h"
#include "vtkMRMLVolumeRenderingSelectionNode.h"

#include <string>

class vtkSlicerVolumeTextureMapper3D;
class vtkFixedPointVolumeRayCastMapper;
class vtkSlicerVRHelper;
class vtkSlicerVolumePropertyWidget;
class VTK_SLICERVOLUMERENDERING_EXPORT vtkVolumeRenderingGUI :public vtkSlicerModuleGUI
{
public:

    static vtkVolumeRenderingGUI *New();
    vtkTypeMacro(vtkVolumeRenderingGUI,vtkSlicerModuleGUI);

    void PrintSelf(ostream& os, vtkIndent indent);

    // Description: Get/Set module logic
    vtkGetObjectMacro (Logic, vtkVolumeRenderingLogic);
    virtual void SetLogic(vtkVolumeRenderingLogic *log)
    {
        this->Logic=log;
        this->Logic->RegisterNodes();
    }

    // Description: Implements vtkSlicerModuleGUI::SetModuleLogic for this GUI
    virtual void SetModuleLogic(vtkSlicerLogic *logic)
    {
      this->SetLogic( dynamic_cast<vtkVolumeRenderingLogic*> (logic) );
    }

    //vtkSetObjectMacro (Logic, vtkVolumeRenderingLogic);
    // Description:
    // Create widgets
    virtual void BuildGUI ( );

    // Description:
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
    // Process events generated by Logic
    virtual void ProcessLogicEvents ( vtkObject *caller, unsigned long event,
        void *callData ){};

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


    // Description:
    // Get/Set the main slicer viewer widget, for picking
    vtkGetObjectMacro(ViewerWidget, vtkSlicerViewerWidget);
    virtual void SetViewerWidget(vtkSlicerViewerWidget *viewerWidget);

    // Description:
    // Get/Set the slicer interactorstyle, for picking
    vtkGetObjectMacro(InteractorStyle, vtkSlicerViewerInteractorStyle);
    virtual void SetInteractorStyle(vtkSlicerViewerInteractorStyle *interactorStyle);

    vtkSetMacro(PipelineInitialized,int);
    vtkGetMacro(PipelineInitialized,int);
    vtkBooleanMacro(PipelineInitialized,int);

    // Description:
    // Get methods on class members ( no Set methods required. )
    vtkGetObjectMacro (PB_HideSurfaceModels,vtkKWPushButton);
    vtkGetObjectMacro (PB_CreateNewVolumeRenderingNode,vtkKWPushButton);
    vtkGetObjectMacro (NS_ImageData,vtkSlicerNodeSelectorWidget);
    vtkGetObjectMacro (NS_VolumeRenderingDataSlicer,vtkSlicerNodeSelectorVolumeRenderingWidget);
    vtkGetObjectMacro (NS_VolumeRenderingDataScene,vtkSlicerNodeSelectorVolumeRenderingWidget);
    vtkGetObjectMacro (EWL_CreateNewVolumeRenderingNode,vtkKWEntryWithLabel);
    vtkGetObjectMacro (DetailsFrame,vtkSlicerModuleCollapsibleFrame);
    vtkGetObjectMacro (CurrentNode,vtkMRMLVolumeRenderingNode);
    vtkGetObjectMacro (Presets, vtkMRMLScene);
    vtkGetObjectMacro (Helper, vtkSlicerVRHelper);

protected:
    vtkVolumeRenderingGUI();
    ~vtkVolumeRenderingGUI();
    vtkVolumeRenderingGUI(const vtkVolumeRenderingGUI&);//not implemented
    void operator=(const vtkVolumeRenderingGUI&);//not implemented

    // Description:
    // Updates GUI widgets based on parameters values in MRML node
    void UpdateGUI();

    // Description:
    // Updates parameters values in MRML node based on GUI widgets 
    void UpdateMRML();

    // Description:
    // GUI elements

    // Description:
    // Pointer to the module's logic class
    vtkVolumeRenderingLogic *Logic;

    vtkMRMLVolumeRenderingSelectionNode *SelectionNode;

    // Description:
    // A pointer back to the viewer widget, useful for picking
    vtkSlicerViewerWidget *ViewerWidget;

    // Description:
    // A poitner to the interactor style, useful for picking
    vtkSlicerViewerInteractorStyle *InteractorStyle;

    int PipelineInitialized;//0=no,1=Yes
    void InitializePipelineNewCurrentNode();
    void InitializePipelineFromMRMLScene();
    void InitializePipelineFromSlicer();
    void InitializePipelineFromImageData();
    void LabelMapInitializePipelineNewCurrentNode();
    void LabelMapInitializePipelineFromMRMLScene();
    void LabelMapInitializePipelineFromSlicer();
    void LabelMapInitializePipelineFromImageData();

    //OWN GUI Elements

    //Frame Save/Load
    vtkKWPushButton *PB_HideSurfaceModels;
    vtkKWPushButton *PB_CreateNewVolumeRenderingNode;
    vtkSlicerNodeSelectorWidget *NS_ImageData;
    //BTX
    std::string PreviousNS_ImageData;
    std::string PreviousNS_VolumeRenderingSlicer;
    std::string PreviousNS_VolumeRenderingDataScene;
    //ETX
    vtkSlicerNodeSelectorVolumeRenderingWidget *NS_VolumeRenderingDataSlicer;
    vtkSlicerNodeSelectorVolumeRenderingWidget *NS_VolumeRenderingDataScene;
    vtkKWEntryWithLabel *EWL_CreateNewVolumeRenderingNode;

    //Frame Details
    vtkSlicerModuleCollapsibleFrame *DetailsFrame;
    
    //Other members
    vtkMRMLVolumeRenderingNode  *CurrentNode;
    vtkMRMLScene *Presets;

    void PackSvpGUI(void);
    void UnpackSvpGUI(void);
    vtkSlicerVRHelper *Helper;
    //0 means grayscale, 1 means LabelMap
    int HelperNumber;
    
    int UpdatingGUI;
    int ProcessingGUIEvents;
    int ProcessingMRMLEvents;
};

#endif
