#include "vtkTumorGrowthStep.h"
#include "vtkTumorGrowthGUI.h"
#include "vtkMRMLTumorGrowthNode.h"

#include "vtkKWWizardWidget.h"
#include "vtkKWWizardWorkflow.h"
#include "vtkKWFrameWithLabel.h"
#include "vtkCallbackCommand.h"
#include "vtkKWPushButton.h"
#include "vtkSlicerApplication.h"
#include "vtkSlicerVolumesLogic.h" 
#include "vtkSlicerVolumesGUI.h" 
#include "vtkTumorGrowthLogic.h"

#include "vtkVolumeTextureMapper3D.h"
#include "vtkPiecewiseFunction.h"
#include "vtkColorTransferFunction.h"
#include "vtkVolumeProperty.h"
#include "vtkVolume.h"
#include "vtkSlicerSliceControllerWidget.h"

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkTumorGrowthStep);
vtkCxxRevisionMacro(vtkTumorGrowthStep, "$Revision: 1.2 $");
vtkCxxSetObjectMacro(vtkTumorGrowthStep,GUI,vtkTumorGrowthGUI);

//----------------------------------------------------------------------------
vtkTumorGrowthStep::vtkTumorGrowthStep()
{
  this->GUI = NULL;
  this->Frame           = NULL;
  this->NextStep = NULL; 
  this->WizardGUICallbackCommand = vtkCallbackCommand::New();
  this->WizardGUICallbackCommand->SetClientData(reinterpret_cast<void *>(this));
  this->GridButton = NULL;
  this->ResetButton = NULL;

  this->Render_Image = NULL;
  this->Render_Mapper = NULL;
  this->Render_Filter = NULL;
  this->Render_ColorMapping = NULL;
  this->Render_VolumeProperty = NULL;
  this->Render_Volume = NULL;
  this->Render_OrientationMatrix = NULL; 
}

//----------------------------------------------------------------------------
vtkTumorGrowthStep::~vtkTumorGrowthStep()
{
  this->SetGUI(NULL);
  if (this->Frame)
  {
    this->Frame->Delete();
    this->Frame = NULL;
  }

  if(this->WizardGUICallbackCommand) 
  {
        this->WizardGUICallbackCommand->Delete();
        this->WizardGUICallbackCommand=NULL;
  }

  if (this->GridButton)
    {
      this->GridButton->Delete();
      this->GridButton = NULL;
    }


  if (this->ResetButton)
    {
      this->ResetButton->Delete();
      this->ResetButton = NULL;
    }


  this->RenderRemove();

}

void vtkTumorGrowthStep::RenderRemove() { 
  if (this->Render_Volume) {
    this->GetGUI()->GetApplicationGUI()->GetViewerWidget()->GetMainViewer()->RemoveViewProp(this->Render_Volume);
    this->Render_Volume->Delete();
    this->Render_Volume = NULL; 
  }

  if (this->Render_Mapper) {
    this->Render_Mapper->Delete();
    this->Render_Mapper = NULL;
  }
  if (this->Render_Filter) {
    this->Render_Filter->Delete();
    this->Render_Filter = NULL;
  }
  if (this->Render_ColorMapping) {
    this->Render_ColorMapping->Delete();
    this->Render_ColorMapping = NULL;
  }
  if (this->Render_VolumeProperty) {
    this->Render_VolumeProperty->Delete();
    this->Render_VolumeProperty = NULL;
  }

  if (this->Render_OrientationMatrix) {
    this->Render_OrientationMatrix->Delete();
    this->Render_OrientationMatrix = NULL; 
  }
  this->Render_Image = NULL;
}

//----------------------------------------------------------------------------
void vtkTumorGrowthStep::HideUserInterface()
{
  this->Superclass::HideUserInterface();

  if (this->GetGUI())
    {
    this->GetGUI()->GetWizardWidget()->ClearPage();
    }
}

//----------------------------------------------------------------------------
void vtkTumorGrowthStep::Validate()
{
  this->Superclass::Validate();

  vtkKWWizardWorkflow *wizard_workflow = 
    this->GetGUI()->GetWizardWidget()->GetWizardWorkflow();

  wizard_workflow->PushInput(vtkKWWizardStep::GetValidationSucceededInput());
  wizard_workflow->ProcessInputs();
}

//----------------------------------------------------------------------------
int vtkTumorGrowthStep::CanGoToSelf()
{
  return this->Superclass::CanGoToSelf() || 1;
}

void vtkTumorGrowthStep::ShowUserInterface()
{
  this->Superclass::ShowUserInterface();
  
  if (this->NextStep) { this->NextStep->RemoveResults(); }


  vtkKWWizardWidget *wizard_widget = this->GetGUI()->GetWizardWidget();
     wizard_widget->GetCancelButton()->SetEnabled(0);
  vtkKWWidget *parent = wizard_widget->GetClientArea();

  if (!this->Frame)
  {
    this->Frame = vtkKWFrameWithLabel::New();
  }
  if (!this->Frame->IsCreated())
    {
    this->Frame->SetParent(parent);
    this->Frame->Create();
    this->Frame->AllowFrameToCollapseOff();
  }

  wizard_widget->NextButtonVisibilityOff();
  wizard_widget->CancelButtonVisibilityOn();
  wizard_widget->GetCancelButton()->SetText("Next >");
  wizard_widget->GetCancelButton()->SetCommand(this, "TransitionCallback");
  wizard_widget->GetCancelButton()->EnabledOn();

  // Does not work 
  // wizard_widget->GetBackButton()->SetCommand(this, "TransitionToPreviousStep");
  // OK Button only is shown at the end  
}
//----------------------------------------------------------------------------
void vtkTumorGrowthStep::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}


void  vtkTumorGrowthStep::GridCallback() {
  vtkMRMLTumorGrowthNode* Node = this->GetGUI()->GetNode();
  if (!Node) return;
 
  vtkMRMLScalarVolumeNode* currentNode =  vtkMRMLScalarVolumeNode::SafeDownCast(Node->GetScene()->GetNodeByID(Node->GetGrid_Ref()));
  if (currentNode) {
    this->GridRemove();
    this->GridButton->SetReliefToRidge();
  }
  else if (this->GridDefine()) {
    this->GridButton->SetReliefToSunken();
  }
  this->GetGUI()->PropagateVolumeSelection();
}

void vtkTumorGrowthStep::CreateGridButton() {
  // Grid Button 
  if (!this->GridButton) {
     this->GridButton = vtkKWPushButton::New();
  }

  if (!this->GridButton->IsCreated()) {
    vtkKWWizardWidget *wizard_widget = this->GetGUI()->GetWizardWidget();
    this->GridButton->SetParent(wizard_widget->GetCancelButton()->GetParent());
    this->GridButton->Create();
    this->GridButton->SetWidth(wizard_widget->GetCancelButton()->GetWidth());
    this->GridButton->SetCommand(this, "GridCallback"); 
    this->GridButton->SetText("Grid");
  }
  this->Script("pack %s -side left -anchor nw -expand n -padx 0 -pady 2", this->GridButton->GetWidgetName()); 

  // Button is hold down if Grid already exists 
  vtkMRMLTumorGrowthNode* Node = this->GetGUI()->GetNode();
  if (!Node) return;

  vtkMRMLScalarVolumeNode* currentNode =  vtkMRMLScalarVolumeNode::SafeDownCast(Node->GetScene()->GetNodeByID(Node->GetGrid_Ref()));
  if (currentNode) {
    this->GridButton->SetReliefToSunken(); 
  }
}
  
void vtkTumorGrowthStep::GridRemove() {
  vtkMRMLTumorGrowthNode* Node = this->GetGUI()->GetNode();
  if (Node) {
    vtkMRMLScalarVolumeNode* currentNode =  vtkMRMLScalarVolumeNode::SafeDownCast(Node->GetScene()->GetNodeByID(Node->GetGrid_Ref()));
    if (currentNode) this->GetGUI()->GetMRMLScene()->RemoveNode(currentNode); 
    vtkSlicerApplicationLogic *applicationLogic = this->GetGUI()->GetLogic()->GetApplicationLogic();
    applicationLogic->GetSelectionNode()->SetReferenceActiveLabelVolumeID(NULL);
    applicationLogic->PropagateVolumeSelection( 0 );
    Node->SetGrid_Ref(NULL);
  }
}


int vtkTumorGrowthStep::GridDefine() {
  // Initialize
  this->GridRemove();

  vtkMRMLTumorGrowthNode* Node = this->GetGUI()->GetNode();
  if (!Node) return 0 ;

  vtkMRMLScene* mrmlScene       =  Node->GetScene();
  vtkMRMLNode* mrmlScan1Node  =  mrmlScene->GetNodeByID(Node->GetScan1_Ref());
  vtkMRMLVolumeNode* volumeNode =  vtkMRMLVolumeNode::SafeDownCast(mrmlScan1Node);
  if (!volumeNode) return 0;
  
  vtkSlicerApplication    *application   = vtkSlicerApplication::SafeDownCast(this->GetApplication());
  vtkSlicerVolumesGUI     *volumesGUI    = vtkSlicerVolumesGUI::SafeDownCast(application->GetModuleGUIByName("Volumes")); 
  vtkSlicerVolumesLogic   *volumesLogic  = volumesGUI->GetLogic();
  vtkMRMLScalarVolumeNode *GridNode      = volumesLogic->CreateLabelVolume(mrmlScene,volumeNode, "TG_Grid");
  Node->SetGrid_Ref(GridNode->GetID());

  vtkSlicerSliceControllerWidget *ControlWidget = this->GetGUI()->GetApplicationGUI()->GetMainSliceGUI0()->GetSliceController();
  double oldOffset = ControlWidget->GetOffsetScale()->GetValue(); 
  vtkSlicerApplicationLogic *applicationLogic = this->GetGUI()->GetLogic()->GetApplicationLogic();
  applicationLogic->GetSelectionNode()->SetReferenceActiveLabelVolumeID(GridNode->GetID());
  applicationLogic->PropagateVolumeSelection( 0 ); 
  ControlWidget->GetOffsetScale()->SetValue(oldOffset); 

  return 1;
}

void vtkTumorGrowthStep::CreateResetButton() {
  // Grid Button 
  if (!this->ResetButton) {
     this->ResetButton = vtkKWPushButton::New();
  }

  if (!this->ResetButton->IsCreated()) {
    vtkKWWizardWidget *wizard_widget = this->GetGUI()->GetWizardWidget();
    this->ResetButton->SetParent(wizard_widget->GetCancelButton()->GetParent());
    this->ResetButton->Create();
    this->ResetButton->SetWidth(wizard_widget->GetCancelButton()->GetWidth());
    this->ResetButton->SetCommand(this->GetGUI(), "PropagateVolumeSelection"); 
    this->ResetButton->SetText("Reset 3D Viewer");
  }
  this->Script("pack %s -side left -anchor nw -expand n -padx 0 -pady 2", this->ResetButton->GetWidgetName()); 
}

//void vtkTumorGrowthStep::SliceLogicDefine() {
//  vtkTumorGrowthGUI *GUI  = this->GetGUI();
//  GUI->SliceLogicDefine();
//  
//  GUI->GetSliceLogic()->GetSliceCompositeNode()->SetReferenceBackgroundVolumeID(GUI->GetNode()->GetScan1_Ref());
//  GUI->GetSliceLogic()->GetSliceNode()->SetFieldOfView(250,250,1);
//  GUI->GetSliceLogic()->SetSliceOffset(GUI->GetSliceController_OffsetScale()->GetValue());
//} 


/// For Rendering results
void vtkTumorGrowthStep::SetRender_BandPassFilter(double min, double max) {
  // cout <<  "SetPreSegment_Render_BandPassFilter " << value << endl;
  double* imgRange  =   this->Render_Image->GetPointData()->GetScalars()->GetRange();
  this->Render_Filter->RemoveAllPoints();
  this->Render_Filter->AddPoint(imgRange[0], 0.0);
  this->Render_Filter->AddPoint(min - 1, 0.0);
  this->Render_Filter->AddPoint(min, 1);
  this->Render_Filter->AddPoint(max, 1);
  if (max < imgRange[1]) { 
    this->Render_Filter->AddPoint(max + 1, 0);
    if (max+1 < imgRange[1]) { 
      this->Render_Filter->AddPoint(imgRange[1], 0);
    }
  }
}

void vtkTumorGrowthStep::SetRender_HighPassFilter(double min) {
  // cout <<  "SetPreSegment_Render_BandPassFilter " << value << endl;
  double* imgRange  =   this->Render_Image->GetPointData()->GetScalars()->GetRange();
  this->Render_Filter->RemoveAllPoints();
  this->Render_Filter->AddPoint(imgRange[0], 0.0);
  this->Render_Filter->AddPoint(min - 1, 0.0);
  this->Render_Filter->AddPoint(min, 1);
}


void vtkTumorGrowthStep::CreateRender(vtkMRMLVolumeNode *volumeNode, float colorR, float colorG, float colorB ) {
  this->RenderRemove();
  if (!volumeNode) return;

  this->Render_Image = volumeNode->GetImageData();
  
  // set PROP [[vtkTumorGrowthAnalysisStep ListInstances] GetRender_Mapper]
  this->Render_Mapper = vtkVolumeTextureMapper3D::New();
  this->Render_Mapper->SetInput(this->Render_Image);

  double* imgRange  = this->Render_Image->GetPointData()->GetScalars()->GetRange();

  this->Render_Filter = vtkPiecewiseFunction::New();
  // this->SetRender_BandPassFilter(range[0],range[1]);

  this->Render_ColorMapping = vtkColorTransferFunction::New();
  this->Render_ColorMapping->AddRGBPoint( imgRange[0] , colorR, colorG, colorB);
  this->Render_ColorMapping->AddRGBPoint( imgRange[1] , colorR, colorG, colorB );

  // set PROP [[vtkTumorGrowthAnalysisStep ListInstances] GetRender_VolumeProperty]

  this->Render_VolumeProperty = vtkVolumeProperty::New();
  this->Render_VolumeProperty->SetShade(1);
  this->Render_VolumeProperty->SetAmbient(0.3);
  this->Render_VolumeProperty->SetDiffuse(0.6);
  this->Render_VolumeProperty->SetSpecular(0.5);
  this->Render_VolumeProperty->SetSpecularPower(40.0);
  this->Render_VolumeProperty->SetScalarOpacity(this->Render_Filter);
  this->Render_VolumeProperty->SetColor( this->Render_ColorMapping );
  this->Render_VolumeProperty->SetInterpolationTypeToLinear();

  this->Render_OrientationMatrix = vtkMatrix4x4::New();
  volumeNode->GetIJKToRASMatrix(this->Render_OrientationMatrix);

  this->Render_Volume = vtkVolume::New();
  this->Render_Volume->SetProperty(this->Render_VolumeProperty);
  this->Render_Volume->SetMapper(this->Render_Mapper);
  this->Render_Volume->PokeMatrix(this->Render_OrientationMatrix);
  
  this->GetGUI()->GetApplicationGUI()->GetViewerWidget()->GetMainViewer()->AddViewProp(this->Render_Volume);
}
