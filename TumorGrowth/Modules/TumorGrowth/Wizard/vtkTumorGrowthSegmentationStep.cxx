#include "vtkTumorGrowthSegmentationStep.h"

#include "vtkTumorGrowthGUI.h"
#include "vtkMRMLTumorGrowthNode.h"

#include "vtkKWWizardWidget.h"
#include "vtkKWWizardWorkflow.h"
#include "vtkKWRange.h"
#include "vtkKWFrameWithLabel.h"
#include "vtkKWLabel.h"
#include "vtkKWEntry.h"
#include "vtkTumorGrowthLogic.h"
#include "vtkSlicerApplicationGUI.h"
#include "vtkSlicerSliceControllerWidget.h"
#include "vtkKWScale.h"
#include "vtkImageAccumulate.h"
#include "vtkImageThreshold.h"
#include "vtkSlicerVolumesLogic.h" 
#include "vtkSlicerVolumesGUI.h"
#include "vtkSlicerApplication.h"
#include "vtkImageIslandFilter.h"

#include "vtkVolumeTextureMapper3D.h"
#include "vtkPiecewiseFunction.h"
#include "vtkColorTransferFunction.h"
#include "vtkVolumeProperty.h"
#include "vtkVolume.h"
#include "vtkSlicerModelsLogic.h"
#include "vtkKWScaleWithEntry.h"

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkTumorGrowthSegmentationStep);
vtkCxxRevisionMacro(vtkTumorGrowthSegmentationStep, "$Revision: 1.2 $");

//----------------------------------------------------------------------------
vtkTumorGrowthSegmentationStep::vtkTumorGrowthSegmentationStep()
{
  this->SetName("3/4. Identify Tumor in First Scan"); 
  this->SetDescription("Move slider to outline boundary of tumor."); 
  this->WizardGUICallbackCommand->SetCallback(vtkTumorGrowthSegmentationStep::WizardGUICallback);

  this->ThresholdFrame = NULL;
  this->ThresholdRange = NULL;
  this->ThresholdLabel = NULL;

  this->PreSegment = NULL;
  this->PreSegmentNode = NULL;
  this->SegmentNode = NULL;

  this->PreSegment_Render_Mapper = NULL;
  this->PreSegment_Render_BandPassFilter = NULL;
  this->PreSegment_Render_ColorMapping = NULL;
  this->PreSegment_Render_VolumeProperty = NULL;
  this->PreSegment_Render_Volume = NULL;
  this->PreSegment_Render_OrientationMatrix = NULL; 

  this->SliceLogic      = NULL;
  this->SliceController_OffsetScale = NULL;
}

//----------------------------------------------------------------------------
vtkTumorGrowthSegmentationStep::~vtkTumorGrowthSegmentationStep()
{
 
  if (this->ThresholdFrame)
    {
    this->ThresholdFrame->Delete();
    this->ThresholdFrame = NULL;
    }

  if (this->ThresholdRange)
    {
    this->ThresholdRange->Delete();
    this->ThresholdRange = NULL;
    }

  if (this->ThresholdLabel)
    {
    this->ThresholdLabel->Delete();
    this->ThresholdLabel = NULL;
    }

  this->PreSegmentScan1Remove();
  this->SegmentScan1Remove();
  this->SliceLogicRemove();
}

//----------------------------------------------------------------------------
void vtkTumorGrowthSegmentationStep::ShowUserInterface()
{
  // ----------------------------------------
  // Display Super Sampled Volume 
  // ----------------------------------------
  
  // cout << "vtkTumorGrowthSegmentationStep::ShowUserInterface()" << endl;

  vtkMRMLTumorGrowthNode* node = this->GetGUI()->GetNode();
  int intMin, intMax;

  if (node) { 
    vtkMRMLVolumeNode *volumeNode = vtkMRMLVolumeNode::SafeDownCast(node->GetScene()->GetNodeByID(node->GetScan1_SuperSampleRef()));
    if (volumeNode) {
      vtkSlicerApplicationLogic *applicationLogic = this->GetGUI()->GetLogic()->GetApplicationLogic();
      applicationLogic->GetSelectionNode()->SetActiveVolumeID(volumeNode->GetID());
      applicationLogic->PropagateVolumeSelection();

      double *imgRange = volumeNode->GetImageData()->GetPointData()->GetScalars()->GetRange();

      intMin = int(imgRange[0]);
      intMax = int(imgRange[1]);
    } else {
      intMin = 0;
      intMax = 0;      
    } 
  } else {
      intMin = 0;
      intMax = 0;      
  }
  //  cout << "What h" << endl;
  // ----------------------------------------
  // Build GUI 
  // ----------------------------------------
  this->vtkTumorGrowthStep::ShowUserInterface();
  this->Frame->SetLabelText("Identify Tumor");
  this->Script("pack %s -side top -anchor nw -fill x -padx 0 -pady 2", this->Frame->GetWidgetName());

  if (!this->ThresholdFrame)
    {
    this->ThresholdFrame = vtkKWFrame::New();
    }
  if (!this->ThresholdFrame->IsCreated())
    {
      this->ThresholdFrame->SetParent(this->Frame->GetFrame());
     this->ThresholdFrame->Create();
  }


  if (!this->ThresholdLabel)
    {
    this->ThresholdLabel = vtkKWLabel::New();
    }
  if (!this->ThresholdLabel->IsCreated())
  {
    this->ThresholdLabel->SetParent(this->ThresholdFrame);
    this->ThresholdLabel->Create();
    this->ThresholdLabel->SetText("Threshold:");
  }

  if (!this->ThresholdRange)
    {
    this->ThresholdRange = vtkKWRange::New();
    }
  if (!this->ThresholdRange->IsCreated())
  {
    this->ThresholdRange->SetParent(this->ThresholdFrame);
    this->ThresholdRange->Create();
    this->ThresholdRange->SymmetricalInteractionOff();
    this->ThresholdRange->SetCommand(this, "ThresholdRangeChangedCallback"); 
    this->ThresholdRange->SetWholeRange(intMin, intMax); 
    this->ThresholdRange->SetResolution(1);

  }

  this->Script("pack %s -side top -anchor nw -padx 0 -pady 3",this->ThresholdFrame->GetWidgetName()); 
  this->Script("pack %s %s -side left -anchor nw -padx 2 -pady 0",this->ThresholdLabel->GetWidgetName(),this->ThresholdRange->GetWidgetName());

  this->CreateGridButton(); 
  // ----------------------------------------
  // Show segmentation 
  // ----------------------------------------
  this->PreSegmentScan1Define();

  {
    vtkMRMLTumorGrowthNode *mrmlNode = this->GetGUI()->GetNode();
    double min, max;
    if (mrmlNode && (mrmlNode->GetSegmentThresholdMin() > -1) && (mrmlNode->GetSegmentThresholdMax() > -1)) {
      min =  mrmlNode->GetSegmentThresholdMin();
      max =  mrmlNode->GetSegmentThresholdMax();
    } else {
      min = (intMax - intMin)/2.0;
      max = intMax;
    }
    this->ThresholdRange->SetRange(min,max);

    // Necesary in order to transfere results from above lines  
    this->ThresholdRangeChangedCallback(min, max);
    // this->TransitionCallback();   
  }
   

  // Show Reference Image 1 in the 3D Slicer Viewer
  this->SliceLogicDefine(); 

 
  // Kilian
  // this->TransitionCallback();   
  //  ~/Slicer/Slicer3/Base/Logic/vtkSlicersModelsLogic Clone 
  // -> copy slice modules
  //  -> look in ~/Slicer/Slicer3/Base/Logic/vtkSlicerSliceLogic 
  //    -> CreateSlicerModel 
  //       -> do a copy 
  // 
  // -> Look up pid 
  // like vtkSlicersVolumeLogic 
  // 
  // 
  // use ddd 
  // or 
  // use gdb
  // -> attach pid 
  // after it loads
  // -> c 
  // after it gives a seg fault just say 
  // -> where
  
  // cout << "End Show user interface" << endl;
}

void vtkTumorGrowthSegmentationStep::PreSegmentScan1Remove() {

  // cout << "vtkTumorGrowthSegmentationStep::PreSegmentScan1Remove() Start " << this->PreSegmentNode << endl;
  if (this->PreSegmentNode) {
    this->GetGUI()->GetMRMLScene()->RemoveNode(this->PreSegmentNode);  
    this->PreSegmentNode = NULL;
  } 
  if (this->PreSegment_Render_Volume) {
    this->GetGUI()->GetApplicationGUI()->GetViewerWidget()->GetMainViewer()->RemoveViewProp(this->PreSegment_Render_Volume);
    this->PreSegment_Render_Volume = NULL; 
  }

  if (this->PreSegment) {
    this->PreSegment->Delete();
    this->PreSegment = NULL;
  }

  if (this->PreSegment_Render_Mapper) {
    this->PreSegment_Render_Mapper->Delete();
    this->PreSegment_Render_Mapper = NULL;
  }
  if (this->PreSegment_Render_BandPassFilter) {
    this->PreSegment_Render_BandPassFilter->Delete();
    this->PreSegment_Render_BandPassFilter = NULL;
  }
  if (this->PreSegment_Render_ColorMapping) {
    this->PreSegment_Render_ColorMapping->Delete();
    this->PreSegment_Render_ColorMapping = NULL;
  }
  if (this->PreSegment_Render_VolumeProperty) {
    this->PreSegment_Render_VolumeProperty->Delete();
    this->PreSegment_Render_VolumeProperty = NULL;
  }
  if (this->PreSegment_Render_Volume) {
    this->PreSegment_Render_Volume->Delete();
    this->PreSegment_Render_Volume = NULL;
  }
  if (this->PreSegment_Render_OrientationMatrix) {
    this->PreSegment_Render_OrientationMatrix->Delete();
    this->PreSegment_Render_OrientationMatrix = NULL; 
  }
  // cout << "vtkTumorGrowthSegmentationStep::PreSegmentScan1Remove() End " << endl;
}

void vtkTumorGrowthSegmentationStep::SetPreSegment_Render_BandPassFilter(double min, double max) {
  // cout <<  "SetPreSegment_Render_BandPassFilter " << value << endl;

  vtkMRMLTumorGrowthNode *mrmlNode = this->GetGUI()->GetNode();
  if (!mrmlNode) return;
  // 3D Render 
  vtkMRMLVolumeNode *volumeNode = vtkMRMLVolumeNode::SafeDownCast(mrmlNode->GetScene()->GetNodeByID(mrmlNode->GetScan1_SuperSampleRef()));
  if (!volumeNode) return;
  double* imgRange  = volumeNode->GetImageData()->GetPointData()->GetScalars()->GetRange();

  this->PreSegment_Render_BandPassFilter->RemoveAllPoints();
  this->PreSegment_Render_BandPassFilter->AddPoint(imgRange[0], 0.0);
  this->PreSegment_Render_BandPassFilter->AddPoint(min - 1, 0.0);
  this->PreSegment_Render_BandPassFilter->AddPoint(min, 1);
  this->PreSegment_Render_BandPassFilter->AddPoint(max, 1);
  if (max < imgRange[1]) { 
    this->PreSegment_Render_BandPassFilter->AddPoint(max + 1, 0);
    if (max+1 < imgRange[1]) { 
      this->PreSegment_Render_BandPassFilter->AddPoint(imgRange[1], 0);
    }
  }
}

void vtkTumorGrowthSegmentationStep::PreSegmentScan1Define() {

  // ---------------------------------
  // Initialize Function
  // ---------------------------------
  vtkMRMLTumorGrowthNode* Node      =  this->GetGUI()->GetNode();
  if (!Node) return;
  vtkMRMLVolumeNode *volumeNode = vtkMRMLVolumeNode::SafeDownCast(Node->GetScene()->GetNodeByID(Node->GetScan1_SuperSampleRef()));
  if (!volumeNode) return;
  vtkMRMLVolumeNode *spgrNode = vtkMRMLVolumeNode::SafeDownCast(Node->GetScene()->GetNodeByID(Node->GetScan1_Ref()));
  if (!spgrNode) return;
  if (!this->ThresholdRange) return;

  vtkSlicerApplication      *application      =  vtkSlicerApplication::SafeDownCast(this->GetApplication()); 
  vtkSlicerApplicationGUI   *applicationGUI   = this->GetGUI()->GetApplicationGUI();
  vtkSlicerApplicationLogic *applicationLogic = this->GetGUI()->GetLogic()->GetApplicationLogic();
  vtkSlicerVolumesLogic     *volumesLogic      = (vtkSlicerVolumesGUI::SafeDownCast(application->GetModuleGUIByName("Volumes")))->GetLogic();

  vtkMRMLScene *mrmlScene = this->GetGUI()->GetMRMLScene();

  if (this->PreSegment || this->PreSegmentNode) this->PreSegmentScan1Remove();

  // ---------------------------------
  // Define LabelMap 
  // ---------------------------------

  this->PreSegment = vtkImageThreshold::New(); 
  int range[2] = {int(this->ThresholdRange->GetRange()[0]),int(this->ThresholdRange->GetRange()[1])}; 
  vtkTumorGrowthLogic::DefinePreSegment(volumeNode->GetImageData(),range,this->PreSegment);

  // ---------------------------------
  // show SPGR in 3D Viewer 
  // ------------------------------


  // ---------------------------------
  // show segmentation in Slice view 
  // ------------------------------
  this->PreSegmentNode = volumesLogic->CreateLabelVolume(Node->GetScene(),volumeNode, "TG_Scan1_PreSegment");
  this->PreSegmentNode->SetAndObserveImageData(this->PreSegment->GetOutput());
  applicationGUI->GetMainSliceGUI0()->GetSliceController()->GetForegroundSelector()->SetSelected(this->PreSegmentNode);
  applicationGUI->GetMainSliceGUI1()->GetSliceController()->GetForegroundSelector()->SetSelected(this->PreSegmentNode);
  applicationGUI->GetMainSliceGUI2()->GetSliceController()->GetForegroundSelector()->SetSelected(this->PreSegmentNode);
  applicationGUI->GetSlicesControlGUI()->GetSliceFadeScale()->SetValue(0.6);
  applicationLogic->PropagateVolumeSelection();
  // ------------------------------------
  // Show Segmentation through 3D Volume Rendering
  // ------------------------------------
  
  this->PreSegment_Render_Mapper = vtkVolumeTextureMapper3D::New();
  this->PreSegment_Render_Mapper->SetInput(volumeNode->GetImageData());

  double* imgRange  = volumeNode->GetImageData()->GetPointData()->GetScalars()->GetRange();
  // cout << "vtkTumorGrowthSegmentationStep::PreSegmentScan1Define()" << 100 << endl;

  this->PreSegment_Render_BandPassFilter = vtkPiecewiseFunction::New();
  this->SetPreSegment_Render_BandPassFilter(range[0],range[1]);

  this->PreSegment_Render_ColorMapping = vtkColorTransferFunction::New();
  this->PreSegment_Render_ColorMapping->AddRGBPoint( imgRange[0] , 0.8, 0.8, 0.0 );
  this->PreSegment_Render_ColorMapping->AddRGBPoint( imgRange[1] , 0.8, 0.8, 0.0 );

  this->PreSegment_Render_VolumeProperty = vtkVolumeProperty::New();
  this->PreSegment_Render_VolumeProperty->SetShade(1);
  this->PreSegment_Render_VolumeProperty->SetAmbient(0.3);
  this->PreSegment_Render_VolumeProperty->SetDiffuse(0.6);
  this->PreSegment_Render_VolumeProperty->SetSpecular(0.5);
  this->PreSegment_Render_VolumeProperty->SetSpecularPower(40.0);
  this->PreSegment_Render_VolumeProperty->SetScalarOpacity(this->PreSegment_Render_BandPassFilter);
  this->PreSegment_Render_VolumeProperty->SetColor( this->PreSegment_Render_ColorMapping );
  this->PreSegment_Render_VolumeProperty->SetInterpolationTypeToLinear();

  this->PreSegment_Render_OrientationMatrix = vtkMatrix4x4::New();
  volumeNode->GetIJKToRASMatrix(this->PreSegment_Render_OrientationMatrix);

  this->PreSegment_Render_Volume = vtkVolume::New();
  this->PreSegment_Render_Volume->SetProperty(this->PreSegment_Render_VolumeProperty);
  this->PreSegment_Render_Volume->SetMapper(this->PreSegment_Render_Mapper);
  this->PreSegment_Render_Volume->PokeMatrix(this->PreSegment_Render_OrientationMatrix);
  
  applicationGUI->GetViewerWidget()->GetMainViewer()->AddViewProp(this->PreSegment_Render_Volume);

  return;
}

void vtkTumorGrowthSegmentationStep::SegmentScan1Remove() {
  vtkMRMLTumorGrowthNode* Node = this->GetGUI()->GetNode();
  if (Node) {
    vtkMRMLVolumeNode* currentNode =  vtkMRMLVolumeNode::SafeDownCast(Node->GetScene()->GetNodeByID(Node->GetScan1_SegmentRef()));
    if (currentNode) this->GetGUI()->GetMRMLScene()->RemoveNode(currentNode); 
    Node->SetScan1_SegmentRef(NULL);
  }
  if (this->SegmentNode) {
    this->SegmentNode->Delete();
    this->SegmentNode = NULL;
  }
}

int vtkTumorGrowthSegmentationStep::SegmentScan1Define() {
  // Initialize
  if (!this->PreSegment || !this->PreSegmentNode) return 0;
  vtkMRMLTumorGrowthNode* Node = this->GetGUI()->GetNode();
  if (!Node) return 0 ;

  this->SegmentScan1Remove();

  vtkImageIslandFilter *RemoveIslands = vtkImageIslandFilter::New();
  vtkTumorGrowthLogic::DefineSegment(this->PreSegment->GetOutput(),RemoveIslands);

  // Set It up 
  vtkSlicerVolumesLogic *volumesLogic         = (vtkSlicerVolumesGUI::SafeDownCast(vtkSlicerApplication::SafeDownCast(this->GetApplication())->GetModuleGUIByName("Volumes")))->GetLogic();

  this->SegmentNode = volumesLogic->CreateLabelVolume(Node->GetScene(),this->PreSegmentNode, "TG_scan1_Segment");
  this->SegmentNode->SetAndObserveImageData(RemoveIslands->GetOutput());

  RemoveIslands->Delete(); 
  this->PreSegmentScan1Remove();

  // Added it to MRML Script
  Node->SetScan1_SegmentRef(this->SegmentNode->GetID());


  return 1;
}


//----------------------------------------------------------------------------
void vtkTumorGrowthSegmentationStep::ThresholdRangeChangedCallback(double min , double max)
{
  if (!this->ThresholdRange || !this->PreSegment) return;
  PreSegment->ThresholdBetween(min,max); 
  PreSegment->Update();
  this->PreSegmentNode->Modified();

  vtkMRMLTumorGrowthNode *mrmlNode = this->GetGUI()->GetNode();
  if (!mrmlNode) return;
  mrmlNode->SetSegmentThresholdMin(min);
  mrmlNode->SetSegmentThresholdMax(max);

  // 3D Render 
  this->SetPreSegment_Render_BandPassFilter(min,max);
  //  applicationGUI->GetViewerWidget()->GetMainViewer()->RequestRender();


  // set GUI  [$::slicer3::Application GetModuleGUIByName "TumorGrowth"]
  // set STEP [$GUI GetSegmentationStep]
  // set FILT [$STEP GetPreSegment]

  // You can also watch MRML by doing 
  // MRMLWatcher m
  // parray MRML
  // $MRML(TG_scan1_SuperSampled) Print

}

void vtkTumorGrowthSegmentationStep::SliceLogicRemove() {
  cout << "vtkTumorGrowthSegmentationStep::SliceLogicRemove" << endl;
  if (this->SliceLogic) {
     vtkSlicerApplicationLogic *applicationLogic = this->GetGUI()->GetLogic()->GetApplicationLogic();
     applicationLogic->GetSlices()->RemoveItem(this->SliceLogic);
     this->SliceLogic->Delete();
     this->SliceLogic = NULL;
  } 
  if (this->SliceController_OffsetScale) {
    this->SliceController_OffsetScale->GetWidget()->RemoveObservers(vtkKWScale::ScaleValueChangedEvent, this->WizardGUICallbackCommand);
    this->SliceController_OffsetScale = NULL;
  }
}

void vtkTumorGrowthSegmentationStep::SliceLogicDefine() {
  if (!this->SliceLogic) {
      vtkIntArray *events = vtkIntArray::New();
      events->InsertNextValue(vtkMRMLScene::NewSceneEvent);
      events->InsertNextValue(vtkMRMLScene::SceneCloseEvent);
      events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
      events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);

      this->SliceLogic = vtkSlicerSliceLogic::New ( );
      this->SliceLogic->SetName("TG");

      this->SliceLogic->SetMRMLScene ( this->GetGUI()->GetMRMLScene());
      this->SliceLogic->ProcessLogicEvents ();
      this->SliceLogic->ProcessMRMLEvents (this->GetGUI()->GetMRMLScene(), vtkCommand::ModifiedEvent, NULL);
      this->SliceLogic->SetAndObserveMRMLSceneEvents (this->GetGUI()->GetMRMLScene(), events );
      events->Delete();

      vtkSlicerApplicationLogic *applicationLogic = this->GetGUI()->GetLogic()->GetApplicationLogic();
      if (applicationLogic->GetSlices())
      {
        applicationLogic->GetSlices()->AddItem(this->SliceLogic);
      }
    } 

    this->SliceLogic->GetSliceNode()->SetSliceVisible(1);
    this->SliceLogic->GetSliceCompositeNode()->SetReferenceBackgroundVolumeID(this->GetGUI()->GetNode()->GetScan1_Ref());
    this->SliceLogic->GetSliceNode()->SetFieldOfView(250,250,1);

    // Link to slicer control pannel 
    this->SliceController_OffsetScale =  this->GetGUI()->GetApplicationGUI()->GetMainSliceGUI0()->GetSliceController()->GetOffsetScale();
    this->SliceLogic->SetSliceOffset(this->SliceController_OffsetScale->GetValue());
    this->SliceController_OffsetScale->GetWidget()->AddObserver(vtkKWScale::ScaleValueChangedEvent, this->WizardGUICallbackCommand);

    // Note : Setting things manually in TCL 
    // [[[vtkTumorGrowthROIStep ListInstances] GetSliceLogic] GetSliceCompositeNode] SetReferenceBackgroundVolumeID vtkMRMLScalarVolumeNode1
    // [[[vtkTumorGrowthROIStep ListInstances] GetSliceLogic] GetSliceNode] SetFieldOfView 200 200 1
} 

//----------------------------------------------------------------------------
void vtkTumorGrowthSegmentationStep::TransitionCallback() 
{ 
  cout << "vtkTumorGrowthSegmentationStep::TransitionCallback()  " << endl;
  this->SegmentScan1Remove();
  if (!this->SegmentScan1Define()) return; 
  vtkSlicerApplication *application   = vtkSlicerApplication::SafeDownCast(this->GetGUI()->GetApplication());
  this->GetGUI()->GetLogic()->SaveVolume(application,this->SegmentNode); 

  this->SliceLogicRemove();
 
  // Proceed to next step 
  this->GUI->GetWizardWidget()->GetWizardWorkflow()->AttemptToGoToNextStep();
}


//----------------------------------------------------------------------------
void vtkTumorGrowthSegmentationStep::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

void  vtkTumorGrowthSegmentationStep::WizardGUICallback(vtkObject *caller, unsigned long event, void *clientData, void *callData )
{
     vtkTumorGrowthSegmentationStep *self = reinterpret_cast< vtkTumorGrowthSegmentationStep *>(clientData);
    if (self) { self->ProcessGUIEvents(caller, event, callData); }


}

void vtkTumorGrowthSegmentationStep::ProcessGUIEvents(vtkObject *caller, unsigned long event, void *callData) {

  if (event == vtkKWScale::ScaleValueChangedEvent) {
    vtkKWScale *scale = vtkKWScale::SafeDownCast(caller);
    if (scale && (scale == this->SliceController_OffsetScale->GetWidget())) 
    { 
      this->SliceLogic->SetSliceOffset(this->SliceController_OffsetScale->GetValue());
    }
  }
}
 
