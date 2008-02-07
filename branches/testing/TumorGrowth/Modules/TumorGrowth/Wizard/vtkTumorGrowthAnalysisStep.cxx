#include "vtkTumorGrowthAnalysisStep.h"

#include "vtkTumorGrowthGUI.h"
#include "vtkMRMLTumorGrowthNode.h"

#include "vtkKWWizardWidget.h"
#include "vtkKWWizardWorkflow.h"
#include "vtkKWThumbWheel.h"

#include "vtkKWFrameWithLabel.h"
#include "vtkKWLabel.h"
#include "vtkKWEntry.h"
#include "vtkSlicerApplicationLogic.h"
#include "vtkTumorGrowthLogic.h"
#include "vtkSlicerSliceControllerWidget.h"
#include "vtkKWScale.h"
#include "vtkSlicerApplication.h"
#include "vtkKWPushButton.h"
#include "vtkKWMessageDialog.h"

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkTumorGrowthAnalysisStep);
vtkCxxRevisionMacro(vtkTumorGrowthAnalysisStep, "$Revision: 1.2 $");

//----------------------------------------------------------------------------
vtkTumorGrowthAnalysisStep::vtkTumorGrowthAnalysisStep()
{
  this->SetName("Analysis"); 
  this->SetDescription("Analysis of Tumor Growth"); 
  this->WizardGUICallbackCommand->SetCallback(vtkTumorGrowthAnalysisStep::WizardGUICallback);

  this->SensitivityScale = NULL;
  this->GrowthLabel = NULL;

  this->ButtonsSave = NULL;
  this->ButtonsSnapshot = NULL;
  this->FrameButtons = NULL;
}

//----------------------------------------------------------------------------
vtkTumorGrowthAnalysisStep::~vtkTumorGrowthAnalysisStep()
{

  if (this->ButtonsSave)
    {
    this->ButtonsSave->Delete();
    this->ButtonsSave  = NULL;
    }

  if (this->ButtonsSnapshot)
    {
    this->ButtonsSnapshot->Delete();
    this->ButtonsSnapshot  = NULL;
    }

  if (this->FrameButtons)
    {
    this->FrameButtons->Delete();
    this->FrameButtons  = NULL;
    }
  
  if (this->SensitivityScale)
    {
    this->SensitivityScale->Delete();
    this->SensitivityScale = NULL;
    }
  if (this->GrowthLabel) 
    {
      this->GrowthLabel->Delete();
      this->GrowthLabel = NULL;
    }


}

//----------------------------------------------------------------------------
void vtkTumorGrowthAnalysisStep::AddGUIObservers() 
{
  // cout << "vtkTumorGrowthROIStep::AddGUIObservers()" << endl; 
  // Make sure you do not add the same event twice - need to do it bc of wizrd structure
  if (this->ButtonsSnapshot && (!this->ButtonsSnapshot->HasObserver(vtkKWPushButton::InvokedEvent, this->WizardGUICallbackCommand))) 
    {
      this->ButtonsSnapshot->AddObserver(vtkKWPushButton::InvokedEvent, this->WizardGUICallbackCommand );  
    } 

  if (this->ButtonsSave && (!this->ButtonsSave->HasObserver(vtkKWPushButton::InvokedEvent, this->WizardGUICallbackCommand))) 
    {
      this->ButtonsSave->AddObserver(vtkKWPushButton::InvokedEvent, this->WizardGUICallbackCommand );  
    } 
}

void vtkTumorGrowthAnalysisStep::RemoveGUIObservers() 
{
  // cout << "vtkTumorGrowthAnalysisStep::RemoveGUIObservers" << endl;
  if (this->ButtonsSnapshot) 
    {
      this->ButtonsSnapshot->RemoveObservers(vtkKWPushButton::InvokedEvent, this->WizardGUICallbackCommand);  
    }

  if (this->ButtonsSave) 
  {
      this->ButtonsSave->RemoveObservers(vtkKWPushButton::InvokedEvent, this->WizardGUICallbackCommand);  
  }
}

void vtkTumorGrowthAnalysisStep::WizardGUICallback(vtkObject *caller, unsigned long event, void *clientData, void *callData )
{
  // cout << "void vtkTumorGrowthAnalysisStep::WizardGUICallback" << endl;
    vtkTumorGrowthAnalysisStep *self = reinterpret_cast<vtkTumorGrowthAnalysisStep *>(clientData);
    if (self) { self->ProcessGUIEvents(caller, event, callData); }


}


void vtkTumorGrowthAnalysisStep::ProcessGUIEvents(vtkObject *caller, unsigned long event, void *callData) {

  // cout << "vtkTumorGrowthAnalysisStep::ProcessGUIEvents" << endl;

  if (event == vtkKWPushButton::InvokedEvent) {
    vtkKWPushButton *button = vtkKWPushButton::SafeDownCast(caller);
    if (this->ButtonsSnapshot && (button == this->ButtonsSnapshot)) 
    { 
      cout << "Not implemented yet Snapshot" << endl;
    }
    else if (this->ButtonsSave && (button == this->ButtonsSave)) 
    { 
      vtkMRMLTumorGrowthNode* node = this->GetGUI()->GetNode();
      if (node) {
    // Save Data 
        vtkMRMLVolumeNode *volumeAnalysisNode = vtkMRMLVolumeNode::SafeDownCast(node->GetScene()->GetNodeByID(node->GetAnalysis_Ref()));
        if (volumeAnalysisNode) {
          vtkTumorGrowthLogic *Logic = this->GetGUI()->GetLogic();
          int oldFlag = Logic->GetSaveVolumeFlag();
      Logic->SetSaveVolumeFlag(1);      
      Logic->SaveVolume(vtkSlicerApplication::SafeDownCast(this->GetGUI()->GetApplication()),volumeAnalysisNode);
      Logic->SetSaveVolumeFlag(oldFlag);      
    }
        // Save MRML 
        node->GetScene()->SetRootDirectory(node->GetWorkingDir());

        std::string fileName(node->GetWorkingDir());
        fileName.append("/Data.mrml");
        node->GetScene()->SetURL(fileName.c_str());

        // Saves file  
        node->GetScene()->Commit();

    std::string infoMsg("Saved Data to ");
    infoMsg.append(node->GetWorkingDir());

        vtkKWMessageDialog::PopupMessage(this->GetGUI()->GetApplication(), this->GetGUI()->GetApplicationGUI()->GetMainSlicerWindow(),
                                         "Tumor Growth",infoMsg.c_str(), vtkKWMessageDialog::OkDefault);

      } else {
    this->GetGUI()->GetApplicationGUI()->ProcessSaveSceneAsCommand();
      }

    }
  }
}

//----------------------------------------------------------------------------
void vtkTumorGrowthAnalysisStep::ShowUserInterface()
{
  // ----------------------------------------
  // Display Analysis Volume 
  // ----------------------------------------  
  vtkMRMLTumorGrowthNode* node = this->GetGUI()->GetNode();
  if (node) { 
    vtkMRMLVolumeNode *volumeSampleNode = vtkMRMLVolumeNode::SafeDownCast(node->GetScene()->GetNodeByID(node->GetScan1_SuperSampleRef()));
    vtkMRMLVolumeNode *volumeAnalysisNode = vtkMRMLVolumeNode::SafeDownCast(node->GetScene()->GetNodeByID(node->GetAnalysis_Ref()));
    if (volumeSampleNode && volumeAnalysisNode) {
      vtkSlicerApplicationLogic *applicationLogic = this->GetGUI()->GetLogic()->GetApplicationLogic();
      applicationLogic->GetSelectionNode()->SetActiveVolumeID(volumeSampleNode->GetID());

      vtkSlicerApplicationGUI *applicationGUI     = this->GetGUI()->GetApplicationGUI();

      double oldSliceSetting[3];
      oldSliceSetting[0] = double(applicationGUI->GetMainSliceGUI0()->GetSliceController()->GetOffsetScale()->GetValue());
      oldSliceSetting[1] = double(applicationGUI->GetMainSliceGUI1()->GetSliceController()->GetOffsetScale()->GetValue());
      oldSliceSetting[2] = double(applicationGUI->GetMainSliceGUI2()->GetSliceController()->GetOffsetScale()->GetValue());

      applicationGUI->GetMainSliceGUI0()->GetSliceController()->GetForegroundSelector()->SetSelected(volumeAnalysisNode);
      applicationGUI->GetMainSliceGUI1()->GetSliceController()->GetForegroundSelector()->SetSelected(volumeAnalysisNode);
      applicationGUI->GetMainSliceGUI2()->GetSliceController()->GetForegroundSelector()->SetSelected(volumeAnalysisNode);
      applicationGUI->GetSlicesControlGUI()->GetSliceFadeScale()->SetValue(0.6);

      applicationLogic->PropagateVolumeSelection();

      // Return to original slice position 
      applicationGUI->GetMainSliceGUI0()->GetSliceController()->GetOffsetScale()->SetValue(oldSliceSetting[0]);
      applicationGUI->GetMainSliceGUI1()->GetSliceController()->GetOffsetScale()->SetValue(oldSliceSetting[1]);
      applicationGUI->GetMainSliceGUI2()->GetSliceController()->GetOffsetScale()->SetValue(oldSliceSetting[2]);
    } 
  }

  // ----------------------------------------
  // Build GUI 
  // ----------------------------------------

  this->vtkTumorGrowthStep::ShowUserInterface();

  this->Frame->SetLabelText("Tumor Growth");
  this->Script("pack %s -side top -anchor nw -fill x -padx 0 -pady 2", this->Frame->GetWidgetName());

  if (!this->SensitivityScale)
    {
    this->SensitivityScale = vtkKWThumbWheel::New();
    }
  if (!this->SensitivityScale->IsCreated())
  {
    this->SensitivityScale->SetParent(this->Frame->GetFrame());
    this->SensitivityScale->Create();
    this->SensitivityScale->SetRange(0.0,1.0);
    this->SensitivityScale->SetMinimumValue(0.0);
    this->SensitivityScale->ClampMinimumValueOn(); 
    this->SensitivityScale->SetMaximumValue(1.0);
    this->SensitivityScale->ClampMaximumValueOn(); 
    this->SensitivityScale->SetResolution(0.75);
    this->SensitivityScale->SetLinearThreshold(1);
    this->SensitivityScale->SetThumbWheelSize (TUMORGROWTH_WIDGETS_SLIDER_WIDTH,TUMORGROWTH_WIDGETS_SLIDER_HEIGHT);
    this->SensitivityScale->DisplayEntryOn();
    this->SensitivityScale->DisplayLabelOn();
    this->SensitivityScale->GetLabel()->SetText("Sensitivity");
    this->SensitivityScale->SetCommand(this,"SensitivityChangedCallback");
    this->SensitivityScale->DisplayEntryAndLabelOnTopOff(); 
    this->SensitivityScale->SetBalloonHelpString("The further the wheel is turned to the right the more robust the result");

    // this->SensitivityScale->GetEntry()->SetCommandTriggerToAnyChange();
  }

  // Initial value 
  vtkMRMLTumorGrowthNode *mrmlNode = this->GetGUI()->GetNode();
  if (mrmlNode) this->SensitivityScale->SetValue(mrmlNode->GetAnalysis_Sensitivity());

  this->Script( "pack %s -side top -anchor nw -padx 2 -pady 2", this->SensitivityScale->GetWidgetName());


  if (!this->GrowthLabel)
    {
    this->GrowthLabel = vtkKWLabel::New();
    }
  if (!this->GrowthLabel->IsCreated())
  {
    this->GrowthLabel->SetParent(this->Frame->GetFrame());
    this->GrowthLabel->Create();
  }
  this->Script( "pack %s -side top -anchor nw -padx 2 -pady 2", this->GrowthLabel->GetWidgetName());

  // Define buttons 
  if (!this->FrameButtons)
  {
    this->FrameButtons = vtkKWFrameWithLabel::New();
  }
  if (!this->FrameButtons->IsCreated())
  {
      vtkKWWizardWidget *wizard_widget = this->GetGUI()->GetWizardWidget();
      this->FrameButtons->SetParent(wizard_widget->GetClientArea());
      this->FrameButtons->Create();
      this->FrameButtons->SetLabelText("Save Analysis");
      this->FrameButtons->AllowFrameToCollapseOff();
  }
  this->Script("pack %s -side top -anchor nw -fill x -padx 0 -pady 2", this->FrameButtons->GetWidgetName());

 if (!this->ButtonsSnapshot) {
    this->ButtonsSnapshot = vtkKWPushButton::New();
  }

  if (!this->ButtonsSnapshot->IsCreated()) {
    this->ButtonsSnapshot->SetParent(this->FrameButtons->GetFrame());
    this->ButtonsSnapshot->Create();
    this->ButtonsSnapshot->SetWidth(TUMORGROWTH_MENU_BUTTON_WIDTH);
    this->ButtonsSnapshot->SetText("Snapshot");
    this->ButtonsSnapshot->EnabledOff();
  }

  if (!this->ButtonsSave) {
    this->ButtonsSave = vtkKWPushButton::New();
  }
  if (!this->ButtonsSave->IsCreated()) {
    this->ButtonsSave->SetParent(this->FrameButtons->GetFrame());
    this->ButtonsSave->Create();
    this->ButtonsSave->SetWidth(TUMORGROWTH_MENU_BUTTON_WIDTH);
    this->ButtonsSave->SetText("Data");
  }

  this->Script("pack %s %s -side left -anchor nw -expand n -padx 2 -pady 2", 
                this->ButtonsSnapshot->GetWidgetName(),this->ButtonsSave->GetWidgetName());

  {
    vtkKWWizardWidget *wizard_widget = this->GetGUI()->GetWizardWidget();
    // wizard_widget->GetOKButton()->SetText("Run");
    wizard_widget->GetCancelButton()->SetText("OK"); 
    wizard_widget->GetCancelButton()->SetCommand(this, "ResetPipelineCallback");
    wizard_widget->GetCancelButton()->EnabledOn();
    wizard_widget->OKButtonVisibilityOff();

  }

  // Show results 
  this->SensitivityChangedCallback(0.0);
  this->AddGUIObservers();
}


//----------------------------------------------------------------------------
void vtkTumorGrowthAnalysisStep::SensitivityChangedCallback(double value)
{
  // Sensitivity has changed because of user interaction
  vtkMRMLTumorGrowthNode *mrmlNode = this->GetGUI()->GetNode();
  if (!this->SensitivityScale || !mrmlNode || !this->GrowthLabel ) return;
  mrmlNode->SetAnalysis_Sensitivity(this->SensitivityScale->GetValue());
  double Growth = this->GetGUI()->GetLogic()->MeassureGrowth(vtkSlicerApplication::SafeDownCast(this->GetGUI()->GetApplication()));
  // show here 
  char TEXT[1024];
  // cout << "---------- " << Growth << " " << mrmlNode->GetSuperSampled_VoxelVolume() << " " << mrmlNode->GetSuperSampled_RatioNewOldSpacing() << endl;;
  sprintf(TEXT,"Growth: %.3f mm^3 (%d Voxels)", Growth*mrmlNode->GetSuperSampled_VoxelVolume(),int(Growth*mrmlNode->GetSuperSampled_RatioNewOldSpacing()));

  this->GrowthLabel->SetText(TEXT);
  // Show updated results 
  vtkMRMLVolumeNode *analysisNode = vtkMRMLVolumeNode::SafeDownCast(mrmlNode->GetScene()->GetNodeByID(mrmlNode->GetAnalysis_Ref()));
  if (analysisNode) analysisNode->Modified();
}

//----------------------------------------------------------------------------
void vtkTumorGrowthAnalysisStep::ResetPipelineCallback() 
{
  // Sensitivity has changed because of user interaction 
  vtkKWWizardWidget *wizard_widget = this->GetGUI()->GetWizardWidget();
  vtkKWWizardWorkflow *wizard_workflow = wizard_widget->GetWizardWorkflow();
  // Go Back to the beginning - you can also make this more generale by first getting the number of states 
  // and then doing a loop 
  wizard_workflow->AttemptToGoToPreviousStep();
  wizard_workflow->AttemptToGoToPreviousStep();
  wizard_workflow->AttemptToGoToPreviousStep();
  wizard_workflow->AttemptToGoToPreviousStep();
}

//----------------------------------------------------------------------------
void vtkTumorGrowthAnalysisStep::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

void  vtkTumorGrowthAnalysisStep::RemoveResults()  { 
    vtkMRMLTumorGrowthNode* Node = this->GetGUI()->GetNode();
    if (!Node) return;
    {
       vtkMRMLVolumeNode* currentNode =  vtkMRMLVolumeNode::SafeDownCast(Node->GetScene()->GetNodeByID(Node->GetAnalysis_Ref()));
       if (currentNode) { this->GetGUI()->GetMRMLScene()->RemoveNode(currentNode); }
    }
}
