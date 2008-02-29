#include <string>
#include <iostream>
#include <sstream>

#include "vtkObjectFactory.h"
#include "vtkTumorGrowthLogic.h"
#include "vtkTumorGrowth.h"
#include "vtkMRMLScene.h"
#include "vtkMRMLTumorGrowthNode.h"
#include "vtkImageClip.h"
#include "vtkImageChangeInformation.h"
#include "vtkImageResample.h"
#include "vtkSlicerApplication.h"
#include "vtkImageMathematics.h"
#include "vtkImageSumOverVoxels.h"
#include "vtkImageIslandFilter.h"
#include "vtkSlicerVolumesLogic.h"
#include "vtkSlicerVolumesGUI.h"
#include "vtkRigidRegistrator.h"
//#include "vtkSlicerApplication.h"

#define ERROR_NODE_VTKID 0

//----------------------------------------------------------------------------
vtkTumorGrowthLogic* vtkTumorGrowthLogic::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = 
    vtkObjectFactory::CreateInstance("vtkTumorGrowthLogic");
  if(ret)
    {
    return (vtkTumorGrowthLogic*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkTumorGrowthLogic;
}


//----------------------------------------------------------------------------
vtkTumorGrowthLogic::vtkTumorGrowthLogic()
{
  this->ModuleName = NULL;

  this->ProgressCurrentAction = NULL;
  this->ProgressGlobalFractionCompleted = 0.0;
  this->ProgressCurrentFractionCompleted = 0.0;

  //this->DebugOn();
  this->TumorGrowthNode = NULL; 
  // this->LocalTransform = NULL; 
  // this->GlobalTransform = NULL; 

  this->Analysis_Intensity_Mean      = 0.0;
  this->Analysis_Intensity_Variance  = 0.0;
  this->Analysis_Intensity_Threshold = 0.0;

  this->Analysis_Intensity_Final          = NULL;
  this->Analysis_Intensity_ROINegativeBin = NULL;
  this->Analysis_Intensity_ROIPositiveBin = NULL;
  this->Analysis_Intensity_ROIBinReal     = NULL;
  this->Analysis_Intensity_ROITotal       = NULL;

  // if set to zero then SaveVolume will not do anything 
  this->SaveVolumeFlag = 1;  
}


//----------------------------------------------------------------------------
vtkTumorGrowthLogic::~vtkTumorGrowthLogic()
{
  vtkSetMRMLNodeMacro(this->TumorGrowthNode, NULL);
  this->SetProgressCurrentAction(NULL);
  this->SetModuleName(NULL);

  // if (this->LocalTransform) {
  //   this->LocalTransform->Delete();
  //   this->LocalTransform = NULL;
  // }
  // 
  // if (this->GlobalTransform) {
  //   this->GlobalTransform->Delete();
  //   this->GlobalTransform = NULL;
  // }
  
  if (this->Analysis_Intensity_Final) {
    this->Analysis_Intensity_Final->Delete();
    this->Analysis_Intensity_Final= NULL;
  }

  if (this->Analysis_Intensity_ROINegativeBin) {
    this->Analysis_Intensity_ROINegativeBin->Delete();
    this->Analysis_Intensity_ROINegativeBin = NULL;
  }

  if (this->Analysis_Intensity_ROIPositiveBin) {
    this->Analysis_Intensity_ROIPositiveBin->Delete();
    this->Analysis_Intensity_ROIPositiveBin = NULL;
  }

  if (this->Analysis_Intensity_ROIBinReal) {
    this->Analysis_Intensity_ROIBinReal->Delete();
    this->Analysis_Intensity_ROIBinReal = NULL;
  }

  if (this->Analysis_Intensity_ROITotal) {
    this->Analysis_Intensity_ROITotal->Delete();
    this->Analysis_Intensity_ROITotal = NULL;
  }

}

//----------------------------------------------------------------------------
void vtkTumorGrowthLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  // !!! todo
}

// For AG
// vtkGeneralTransform* vtkTumorGrowthLogic::CreateGlobalTransform() 
// vtkTransform* vtkTumorGrowthLogic::CreateGlobalTransform() 
// {
//   this->GlobalTransform = vtkTransform::New();
//   return this->GlobalTransform;
// }
// 
// // For AG
// // vtkGeneralTransform* vtkTumorGrowthLogic::CreateLocalTransform() 
// vtkTransform* vtkTumorGrowthLogic::CreateLocalTransform() 
// {
//   this->LocalTransform = vtkTransform::New();
//   return this->LocalTransform;
// }

int vtkTumorGrowthLogic::CheckROI(vtkMRMLVolumeNode* volumeNode) {
  if (!volumeNode || !this->TumorGrowthNode) return 0;

  int* dimensions = volumeNode->GetImageData()->GetDimensions();

  for (int i = 0 ; i < 3 ; i++) {
    if ((this->TumorGrowthNode->GetROIMax(i) < 0) || (this->TumorGrowthNode->GetROIMax(i) >= dimensions[i])) return 0 ;
    if ((this->TumorGrowthNode->GetROIMin(i) < 0) || (this->TumorGrowthNode->GetROIMin(i) >= dimensions[i])) return 0 ;
    if (this->TumorGrowthNode->GetROIMax(i) < this->TumorGrowthNode->GetROIMin(i)) return 0;
  }
  return 1;
}


// Give another Volume as an example - creates a volume without anything in it 
// copied from vtkMRMLScalarVolumeNode* vtkSlicerVolumesLogic::CloneVolume without deep copy
vtkMRMLScalarVolumeNode* vtkTumorGrowthLogic::CreateVolumeNode(vtkMRMLVolumeNode *volumeNode, char *name) {
  if (!this->TumorGrowthNode || !volumeNode ) 
    {
    return NULL;
    }

  // clone the display node
  vtkMRMLScalarVolumeDisplayNode *clonedDisplayNode = vtkMRMLScalarVolumeDisplayNode::New();
  clonedDisplayNode->CopyWithScene(volumeNode->GetDisplayNode());
  this->TumorGrowthNode->GetScene()->AddNode(clonedDisplayNode);

  // clone the volume node
  vtkMRMLScalarVolumeNode *clonedVolumeNode = vtkMRMLScalarVolumeNode::New();
  clonedVolumeNode->CopyWithScene(volumeNode);
  clonedVolumeNode->SetStorageNodeID(NULL);
  clonedVolumeNode->SetName(name);
  clonedVolumeNode->SetAndObserveDisplayNodeID(clonedDisplayNode->GetID());

  if (0) {
    // copy over the volume's data
    vtkImageData* clonedVolumeData = vtkImageData::New(); 
    clonedVolumeData->DeepCopy(volumeNode->GetImageData());
    clonedVolumeNode->SetAndObserveImageData( clonedVolumeData );
    clonedVolumeNode->SetModifiedSinceRead(1);
    clonedVolumeData->Delete();
  }

  // add the cloned volume to the scene
  this->TumorGrowthNode->GetScene()->AddNode(clonedVolumeNode);

  // remove references
  clonedVolumeNode->Delete();
  clonedDisplayNode->Delete();

  return (clonedVolumeNode);
 
}

void vtkTumorGrowthLogic::DeleteSuperSample(int ScanNum) {
  cout <<  "vtkTumorGrowthLogic::DeleteSuperSample " << endl;
   // Delete old attached node first 
  vtkMRMLVolumeNode* currentNode = NULL; 
  if (ScanNum ==1) {
    currentNode =  vtkMRMLVolumeNode::SafeDownCast(this->TumorGrowthNode->GetScene()->GetNodeByID(this->TumorGrowthNode->GetScan1_SuperSampleRef()));
    this->TumorGrowthNode->SetScan1_SuperSampleRef(NULL);
  } else {
    currentNode =  vtkMRMLVolumeNode::SafeDownCast(this->TumorGrowthNode->GetScene()->GetNodeByID(this->TumorGrowthNode->GetScan2_SuperSampleRef()));
    this->TumorGrowthNode->SetScan2_SuperSampleRef(NULL);
  } 
  if (currentNode) { 
    this->TumorGrowthNode->GetScene()->RemoveNode(currentNode); 
  } 
}

double vtkTumorGrowthLogic::DefineSuperSampleSize(const double inputSpacing[3], const int ROIMin[3], const int ROIMax[3]) {
    int size = ROIMax[0] - ROIMin[0] + 1;
    double TempSpacing = double(size) * inputSpacing[0] / 100.0;
    double SuperSampleSpacing = (TempSpacing < 0.3 ?  0.3 : TempSpacing);
    
    size = ROIMax[1] - ROIMin[1] + 1;
    TempSpacing = double(size) * inputSpacing[1] / 100.0;
    if (TempSpacing > SuperSampleSpacing) { SuperSampleSpacing = TempSpacing;}

    size = ROIMax[2] - ROIMin[2] + 1;
    TempSpacing = double(size) * inputSpacing[2] / 100.0;
    if (TempSpacing > SuperSampleSpacing) { SuperSampleSpacing = TempSpacing;}
    
    return SuperSampleSpacing;
}

int vtkTumorGrowthLogic::CreateSuperSampleFct(vtkImageData *input, const int ROIMin[3], const int ROIMax[3], const double SuperSampleSpacing, vtkImageData *output) {
  if (SuperSampleSpacing <= 0.0) return 1;
  // ---------------------------------
  // Just focus on region of interest
  vtkImageClip  *ROI = vtkImageClip::New();
     ROI->SetInput(input);
     ROI->SetOutputWholeExtent(ROIMin[0],ROIMax[0],ROIMin[1],ROIMax[1],ROIMin[2],ROIMax[2]); 
     ROI->ClipDataOn();   
     ROI->Update(); 

  vtkImageChangeInformation *ROIExtent = vtkImageChangeInformation::New();
     ROIExtent->SetInput(ROI->GetOutput());
     ROIExtent->SetOutputExtentStart(0,0,0); 
  ROIExtent->Update();
 
  // ---------------------------------
  // Now perform super sampling 
  vtkImageResample *ROISuperSample = vtkImageResample::New(); 
     ROISuperSample->SetDimensionality(3);
     ROISuperSample->SetInterpolationModeToLinear();
     ROISuperSample->SetInput(ROIExtent->GetOutput());
     ROISuperSample->SetAxisOutputSpacing(0,SuperSampleSpacing);
     ROISuperSample->SetAxisOutputSpacing(1,SuperSampleSpacing);
     ROISuperSample->SetAxisOutputSpacing(2,SuperSampleSpacing);
     ROISuperSample->ReleaseDataFlagOff();
  ROISuperSample->Update();

  // ---------------------------------
  // Clean up 
  output->DeepCopy(ROISuperSample->GetOutput());

  ROISuperSample->Delete();
  ROIExtent->Delete();
  ROI->Delete();
  return 0;
}

vtkMRMLScalarVolumeNode* vtkTumorGrowthLogic::CreateSuperSample(int ScanNum) {
  // ---------------------------------
  // Initialize Variables 
  if (!this->TumorGrowthNode)  return NULL;

  vtkMRMLVolumeNode* volumeNode = NULL;

  if (ScanNum > 1)  {
    volumeNode = vtkMRMLVolumeNode::SafeDownCast(this->TumorGrowthNode->GetScene()->GetNodeByID(this->TumorGrowthNode->GetScan2_GlobalRef()));
  } else {
    volumeNode = vtkMRMLVolumeNode::SafeDownCast(this->TumorGrowthNode->GetScene()->GetNodeByID(this->TumorGrowthNode->GetScan1_Ref()));
  }

  if (!this->CheckROI(volumeNode)) return NULL;

  // ---------------------------------
  // Perform Super Sampling 

  int ROIMin[3] = {this->TumorGrowthNode->GetROIMin(0), this->TumorGrowthNode->GetROIMin(1), this->TumorGrowthNode->GetROIMin(2)};
  int ROIMax[3] = {this->TumorGrowthNode->GetROIMax(0), this->TumorGrowthNode->GetROIMax(1), this->TumorGrowthNode->GetROIMax(2)};

  double SuperSampleSpacing = -1;
  if (ScanNum == 1) {
    double *Spacing = volumeNode->GetSpacing();
    SuperSampleSpacing = this->DefineSuperSampleSize(Spacing, ROIMin, ROIMax);
    double SuperSampleVol = SuperSampleSpacing*SuperSampleSpacing*SuperSampleSpacing;
    this->TumorGrowthNode->SetSuperSampled_Spacing(SuperSampleSpacing);    
    this->TumorGrowthNode->SetSuperSampled_VoxelVolume(SuperSampleVol); 
    this->TumorGrowthNode->SetSuperSampled_RatioNewOldSpacing(SuperSampleVol/(Spacing[0]*Spacing[1]*Spacing[2]));
    this->TumorGrowthNode->SetScan1_VoxelVolume(Spacing[0]*Spacing[1]*Spacing[2]);
  } else {
    SuperSampleSpacing = this->TumorGrowthNode->GetSuperSampled_Spacing();
  }

  vtkImageChangeInformation *ROISuperSampleInput = vtkImageChangeInformation::New();
     ROISuperSampleInput->SetInput(volumeNode->GetImageData());
     ROISuperSampleInput->SetOutputSpacing(volumeNode->GetSpacing());
  ROISuperSampleInput->Update();

  vtkImageData *ROISuperSampleOutput = vtkImageData::New();
  if (this->CreateSuperSampleFct(ROISuperSampleInput->GetOutput(), ROIMin, ROIMax, SuperSampleSpacing, ROISuperSampleOutput)) {
    ROISuperSampleInput->Delete();
    ROISuperSampleOutput->Delete();
    return NULL;
  }

  vtkImageChangeInformation *ROISuperSampleExtent = vtkImageChangeInformation::New();
     ROISuperSampleExtent->SetInput(ROISuperSampleOutput);
     ROISuperSampleExtent->SetOutputSpacing(1,1,1);
  ROISuperSampleExtent->Update();

  // ---------------------------------
  // Now return results and clean up 
  char VolumeOutputName[1024];
  if (ScanNum > 1) sprintf(VolumeOutputName, "TG_scan2_Global_SuperSampled");
  else sprintf(VolumeOutputName, "TG_scan1_SuperSampled");

  vtkMRMLScalarVolumeNode *VolumeOutputNode = this->CreateVolumeNode(volumeNode,VolumeOutputName);
  // VolumeOutputNode->SetAndObserveImageData(ROIExtent->GetOutput());

  VolumeOutputNode->SetAndObserveImageData(ROISuperSampleExtent->GetOutput());
  VolumeOutputNode->SetSpacing(SuperSampleSpacing,SuperSampleSpacing,SuperSampleSpacing); 

  // Compute new rjk matrix 
  // double IJKToRASDirections[3][3];
  // Set new orgin
  vtkMatrix4x4 *ijkToRAS=vtkMatrix4x4::New();
  volumeNode->GetIJKToRASMatrix(ijkToRAS);
  double newIJKOrigin[4] = {this->TumorGrowthNode->GetROIMin(0),this->TumorGrowthNode->GetROIMin(1), this->TumorGrowthNode->GetROIMin(2), 1.0 };
  double newRASOrigin[4];
  ijkToRAS->MultiplyPoint(newIJKOrigin,newRASOrigin);
  VolumeOutputNode->SetOrigin(newRASOrigin[0],newRASOrigin[1],newRASOrigin[2]);

  ROISuperSampleExtent->Delete();
  ROISuperSampleOutput->Delete();
  ROISuperSampleInput->Delete();

  return VolumeOutputNode;
}

  // In tcl
  // set GUI  [$::slicer3::Application GetModuleGUIByName "TumorGrowth"]
  //    set NODE [$GUI  GetNode]
  //    set SCENE [$NODE GetScene]
  // set VolumeOutputNode [$SCENE GetNodeByID [$NODE GetScan1_SuperSampleRef]]
  // set VolumeOutputDisplayNode [$VolumeOutputNode GetScalarVolumeDisplayNode]
  // Important files  
  // ~/Slicer/Slicer3/Base/Logic/vtkSlicerVolumesLogic
  // ~/Slicer/Slicer3/Libs/MRML/vtkMRMLVolumeNode.h
  
  // this->TumorGrowthNode->GetScene()->AddNode(VolumeOutputNode);
  // VolumeOutputNode->Delete();


void vtkTumorGrowthLogic::SourceAnalyzeTclScripts(vtkKWApplication *app) {
 char TCL_FILE[1024]; 
 // Kilian: Can we copy this over to the build directory
 // cout - later, when it works do it this way bc more 
 sprintf(TCL_FILE,"%s/Modules/TumorGrowth/tcl/TumorGrowthFct.tcl",vtksys::SystemTools::GetEnv("SLICER_HOME"));
 // sprintf(TCL_FILE,"%s/../Slicer3/Modules/TumorGrowth/tcl/TumorGrowthFct.tcl",vtksys::SystemTools::GetEnv("SLICER_HOME"));

 app->LoadScript(TCL_FILE); 
 // later do it this way 
 sprintf(TCL_FILE,"%s/Modules/TumorGrowth/tcl/TumorGrowthReg.tcl",vtksys::SystemTools::GetEnv("SLICER_HOME"));
 // sprintf(TCL_FILE,"%s/../Slicer3/Modules/TumorGrowth/tcl/TumorGrowthReg.tcl",vtksys::SystemTools::GetEnv("SLICER_HOME"));
 app->LoadScript(TCL_FILE); 
}

void vtkTumorGrowthLogic::DeleteAnalyzeOutput(vtkSlicerApplication *app) {
   // Delete old attached node first 
  if (!TumorGrowthNode) return;
  this->SourceAnalyzeTclScripts(app);

  app->Script("::TumorGrowthTcl::Scan2ToScan1Registration_DeleteOutput Global");

  vtkMRMLVolumeNode* currentNode =  vtkMRMLVolumeNode::SafeDownCast(this->TumorGrowthNode->GetScene()->GetNodeByID(this->TumorGrowthNode->GetScan2_SuperSampleRef()));
  if (currentNode) { 
    this->TumorGrowthNode->GetScene()->RemoveNode(currentNode); 
    this->TumorGrowthNode->SetScan2_SuperSampleRef(NULL);
  }

  app->Script("::TumorGrowthTcl::HistogramNormalization_DeleteOutput"); 
  app->Script("::TumorGrowthTcl::Scan2ToScan1Registration_DeleteOutput Local"); 
  app->Script("::TumorGrowthTcl::IntensityThresholding_DeleteOutput 1");
  app->Script("::TumorGrowthTcl::IntensityThresholding_DeleteOutput 2");
  app->Script("::TumorGrowthTcl::Analysis_Intensity_DeleteOutput_GUI"); 
}

int vtkTumorGrowthLogic::AnalyzeGrowth(vtkSlicerApplication *app) {
  // This is for testing how to start a tcl script 
  cout << "=== Start ANALYSIS ===" << endl;

  this->SourceAnalyzeTclScripts(app);
  
  if (0) { 
    cout << "=== 1 ===" << endl;
    app->Script("::TumorGrowthTcl::Scan2ToScan1Registration_GUI Global");


    //----------------------------------------------
    // Second step -> Save the outcome
    if (!this->TumorGrowthNode) {return 0;}
  
    this->DeleteSuperSample(2);
    vtkMRMLScalarVolumeNode *outputNode = this->CreateSuperSample(2);
    if (!outputNode) {return 0;} 
    this->TumorGrowthNode->SetScan2_SuperSampleRef(outputNode->GetID());
    this->SaveVolume(app,outputNode);

    //----------------------------------------------

    cout << "=== 2 ===" << endl;
    app->Script("::TumorGrowthTcl::HistogramNormalization_GUI"); 
    cout << "=== 3 ===" << endl;
    app->Script("::TumorGrowthTcl::Scan2ToScan1Registration_GUI Local"); 
  } else {
    if (!this->TumorGrowthNode->GetScan2_LocalRef() || !strcmp(this->TumorGrowthNode->GetScan2_LocalRef(),"")) { 
      char fileName[1024];
      sprintf(fileName,"%s/TG_scan2_Local.nhdr",this->TumorGrowthNode->GetWorkingDir());
      vtkMRMLVolumeNode* tmp = this->LoadVolume(app,fileName,0,"TG_scan2_Local");
      if (tmp) {
        this->TumorGrowthNode->SetScan2_LocalRef(tmp->GetID());
      } else {
         cout << "Error: Could not load " << fileName << endl;
         return 0;
      }
    }
  }
  if (this->TumorGrowthNode->GetAnalysis_Intensity_Flag()) { 
    cout << "=== 4 ===" << endl;
    if (!atoi(app->Script("::TumorGrowthTcl::IntensityThresholding_GUI 1"))) return 0; 
    cout << "=== 5 ===" << endl;
    if (!atoi(app->Script("::TumorGrowthTcl::IntensityThresholding_GUI 2"))) return 0; 
    cout << "=== INTENSITY ANALYSIS ===" << endl;
    if (!atoi(app->Script("::TumorGrowthTcl::Analysis_Intensity_GUI"))) return 0; 
  } 
  if (this->TumorGrowthNode->GetAnalysis_Deformable_Flag()) {
    cout << "=== DEFORMABLE ANALYSIS ===" << endl;
    if (!atoi(app->Script("::TumorGrowthTcl::Analysis_Deformable_GUI"))) return 0; 
  }
  cout << "=== End ANALYSIS ===" << endl;

  return 1;
}

void vtkTumorGrowthLogic::RegisterMRMLNodesWithScene() {
   vtkMRMLTumorGrowthNode* tmNode =  vtkMRMLTumorGrowthNode::New();
   this->GetMRMLScene()->RegisterNodeClass(tmNode);
   tmNode->Delete();
}

vtkImageThreshold* vtkTumorGrowthLogic::CreateAnalysis_Intensity_Final() {
  if (this->Analysis_Intensity_Final) { this->Analysis_Intensity_Final->Delete(); }
  this->Analysis_Intensity_Final = vtkImageThreshold::New();
  return this->Analysis_Intensity_Final;
}

vtkImageThreshold* vtkTumorGrowthLogic::CreateAnalysis_Intensity_ROINegativeBin() {
  if (this->Analysis_Intensity_ROINegativeBin) { this->Analysis_Intensity_ROINegativeBin->Delete(); }
  this->Analysis_Intensity_ROINegativeBin = vtkImageThreshold::New();
  return this->Analysis_Intensity_ROINegativeBin;
}

vtkImageThreshold* vtkTumorGrowthLogic::CreateAnalysis_Intensity_ROIPositiveBin() {
  if (this->Analysis_Intensity_ROIPositiveBin) { this->Analysis_Intensity_ROIPositiveBin->Delete(); }
  this->Analysis_Intensity_ROIPositiveBin = vtkImageThreshold::New();
  return this->Analysis_Intensity_ROIPositiveBin;
}

vtkImageMathematics* vtkTumorGrowthLogic::CreateAnalysis_Intensity_ROIBinReal() {
  if (this->Analysis_Intensity_ROIBinReal) { this->Analysis_Intensity_ROIBinReal->Delete(); }
  this->Analysis_Intensity_ROIBinReal = vtkImageMathematics::New();
  return this->Analysis_Intensity_ROIBinReal;
}

vtkImageSumOverVoxels* vtkTumorGrowthLogic::CreateAnalysis_Intensity_ROITotal() {
  if (this->Analysis_Intensity_ROITotal) { this->Analysis_Intensity_ROITotal->Delete(); }
  this->Analysis_Intensity_ROITotal = vtkImageSumOverVoxels::New();
  return this->Analysis_Intensity_ROITotal;
}

double vtkTumorGrowthLogic::GetAnalysis_Intensity_ROITotal_VoxelSum() {
  return this->Analysis_Intensity_ROITotal->GetVoxelSum();
}

vtkImageData*  vtkTumorGrowthLogic::GetAnalysis_Intensity_ROIBinReal() { 
  return (this->Analysis_Intensity_ROIBinReal ? this->Analysis_Intensity_ROIBinReal->GetOutput() : NULL);
}


double vtkTumorGrowthLogic::MeassureGrowth() {
  
  if (!this->Analysis_Intensity_Final || !this->Analysis_Intensity_ROINegativeBin || !this->Analysis_Intensity_ROIPositiveBin || !this->Analysis_Intensity_ROITotal) return -1;
  // Just for display 
  
  this->Analysis_Intensity_Final->ThresholdByUpper(this->Analysis_Intensity_Threshold); 
  this->Analysis_Intensity_Final->Update();
  this->Analysis_Intensity_ROINegativeBin->ThresholdByLower(-this->Analysis_Intensity_Threshold); 
  this->Analysis_Intensity_ROINegativeBin->Update(); 
  this->Analysis_Intensity_ROIPositiveBin->ThresholdByUpper(this->Analysis_Intensity_Threshold); 
  this->Analysis_Intensity_ROIPositiveBin->Update(); 
  this->Analysis_Intensity_ROITotal->Update(); 
  return this->Analysis_Intensity_ROITotal->GetVoxelSum(); 
}

void vtkTumorGrowthLogic::SaveVolume(vtkSlicerApplication *app, vtkMRMLVolumeNode *volNode) {
  if (!this->SaveVolumeFlag) return;  
  this->SaveVolumeForce(app,volNode);
}


void vtkTumorGrowthLogic::SaveVolumeFileName( vtkMRMLVolumeNode *volNode, char *FileName) {
  sprintf(FileName,"%s/%s.nhdr",this->TumorGrowthNode->GetWorkingDir(),volNode->GetName());
}

void vtkTumorGrowthLogic::SaveVolumeForce(vtkSlicerApplication *app, vtkMRMLVolumeNode *volNode) {
 // Initialize
 vtkSlicerVolumesGUI  *volumesGUI    = vtkSlicerVolumesGUI::SafeDownCast(app->GetModuleGUIByName("Volumes")); 
 if (!volumesGUI) return;
 vtkSlicerVolumesLogic *volumesLogic = volumesGUI->GetLogic();

 // Create Directory if necessary
 {
   char CMD[1024];
   sprintf(CMD,"file isdirectory %s",this->TumorGrowthNode->GetWorkingDir()); 
   if (!atoi(app->Script(CMD))) { 
     sprintf(CMD,"file mkdir %s",this->TumorGrowthNode->GetWorkingDir()); 
     app->Script(CMD); 
   }
 }

 {
   char fileName[1024];
   this->SaveVolumeFileName(volNode,fileName);
   cout << "vtkTumorGrowthLogic::SaveVolume: Saving File :" << fileName << endl;
   if (!volumesLogic->SaveArchetypeVolume( fileName, volNode ) )  {
     cout << "Error: Could no save file " << endl;
   }
 }
}

vtkMRMLVolumeNode* vtkTumorGrowthLogic::LoadVolume(vtkSlicerApplication *app, char* fileName, int LabelMapFlag,char* volumeName) {
   vtkSlicerVolumesGUI  *volumesGUI    = vtkSlicerVolumesGUI::SafeDownCast(app->GetModuleGUIByName("Volumes")); 
   if (!volumesGUI) return NULL;
   vtkSlicerVolumesLogic *volumesLogic = volumesGUI->GetLogic();
   // Ignore error messages - I do not know how to get around them 
   return volumesLogic->AddArchetypeVolume(fileName,0,LabelMapFlag,volumeName);
}




//----------------------------------------------------------------------------
void vtkTumorGrowthLogic::PrintResult(ostream& os, vtkSlicerApplication *app)
{  
  // vtkMRMLNode::PrintSelf(os,indent);
  if (!this->TumorGrowthNode) return;
  os  << "This file was generated by vtkMrmTumorGrowthNode " << "\n";;
  os  << "Date:      " << app->Script("date") << "\n";;

  vtkMRMLVolumeNode *VolNode = vtkMRMLVolumeNode::SafeDownCast(this->TumorGrowthNode->GetScene()->GetNodeByID(this->TumorGrowthNode->GetScan1_Ref()));
  os  << "Scan1_Ref: " <<  (VolNode && VolNode->GetStorageNode() ? VolNode->GetStorageNode()->GetFileName() : "(none)") << "\n";

  VolNode = vtkMRMLVolumeNode::SafeDownCast(this->TumorGrowthNode->GetScene()->GetNodeByID(this->TumorGrowthNode->GetScan2_Ref()));
  os  << "Scan2_Ref: " <<  (VolNode && VolNode->GetStorageNode() ? VolNode->GetStorageNode()->GetFileName() : "(none)") << "\n";
  os  << "ROI:" << endl;
  os  << "  Min: " << this->TumorGrowthNode->GetROIMin(0) << " "<< this->TumorGrowthNode->GetROIMin(1) << " "<< this->TumorGrowthNode->GetROIMin(2) <<"\n";
  os  << "  Max: " << this->TumorGrowthNode->GetROIMax(0) << " "<< this->TumorGrowthNode->GetROIMax(1) << " "<< this->TumorGrowthNode->GetROIMax(2) <<"\n";
  os  << "Threshold: [" << this->TumorGrowthNode->GetSegmentThresholdMin() <<", " << this->TumorGrowthNode->GetSegmentThresholdMax() << "]\n";
  if (this->TumorGrowthNode->GetAnalysis_Intensity_Flag()) {
    os  << "Analysis based on Intensity Pattern" << endl;
    os  << "  Sensitivity:      "<< this->TumorGrowthNode->GetAnalysis_Intensity_Sensitivity() << "\n";
    app->Script("::TumorGrowthTcl::Analysis_Intensity_UpdateThreshold_GUI"); 
    double Growth = this->MeassureGrowth(); 
    os  << "  Intensity Metric: "<<  floor(Growth*this->TumorGrowthNode->GetSuperSampled_VoxelVolume()*1000)/1000.0 << "mm" << char(179) 
       << " (" << int(Growth*this->TumorGrowthNode->GetSuperSampled_RatioNewOldSpacing()) << " Voxels)" << "\n";
  }
  if (this->TumorGrowthNode->GetAnalysis_Deformable_Flag()) {
    os  << "Analysis based on Deformable Map" << endl;
    os  << "  Segmentation Metric: "<<  floor(this->TumorGrowthNode->GetAnalysis_Deformable_SegmentationGrowth()*1000)/1000.0 << "mm" << char(179) 
       << " (" << int(this->TumorGrowthNode->GetAnalysis_Deformable_SegmentationGrowth()/this->TumorGrowthNode->GetScan1_VoxelVolume()) << " Voxels)\n";
    os  << "  Jacobian Metric:     "<<  floor(this->TumorGrowthNode->GetAnalysis_Deformable_JacobianGrowth()*1000)/1000.0 << "mm" << char(179) 
       << " (" << int(this->TumorGrowthNode->GetAnalysis_Deformable_JacobianGrowth()/this->TumorGrowthNode->GetScan1_VoxelVolume()) << " Voxels)\n";
  }
}

// works for running stuff in TCL so that you do not need to look in two windows 
void vtkTumorGrowthLogic::PrintText(char *TEXT) {
  cout << TEXT << endl;
} 

void vtkTumorGrowthLogic::DefinePreSegment(vtkImageData *INPUT, const int RANGE[2], vtkImageThreshold *OUTPUT) {
  OUTPUT->SetInValue(10);
  OUTPUT->SetOutValue(0);
  OUTPUT->SetOutputScalarTypeToShort();
  OUTPUT->SetInput(INPUT); 
  OUTPUT->ThresholdBetween(RANGE[0],RANGE[1]); 
  OUTPUT->Update();
}

void vtkTumorGrowthLogic::DefineSegment(vtkImageData *INPUT, vtkImageIslandFilter *OUTPUT) {
  OUTPUT->SetIslandMinSize(1000);
  OUTPUT->SetInput(INPUT);
  OUTPUT->SetNeighborhoodDim3D();
  OUTPUT->Update(); 
}

// Stole it from vtkEMSegmentLogic
void vtkTumorGrowthLogic::RigidRegistration(vtkMRMLVolumeNode* fixedVolumeNode, vtkMRMLVolumeNode* movingVolumeNode, 
                        vtkMRMLVolumeNode* outputVolumeNode, vtkTransform* fixedRASToMovingRASTransform, 
                        double backgroundLevel)
{
  vtkRigidRegistrator* registrator = vtkRigidRegistrator::New();

  // set fixed image ------
  registrator->SetFixedImage(fixedVolumeNode->GetImageData());
  vtkMatrix4x4* IJKToRASMatrixFixed = vtkMatrix4x4::New();
  fixedVolumeNode->GetIJKToRASMatrix(IJKToRASMatrixFixed);
  registrator->SetFixedIJKToXYZ(IJKToRASMatrixFixed);
  IJKToRASMatrixFixed->Delete();
    
  // set moving image ------
  registrator->SetMovingImage(movingVolumeNode->GetImageData());
  vtkMatrix4x4* IJKToRASMatrixMoving = vtkMatrix4x4::New();
  movingVolumeNode->GetIJKToRASMatrix(IJKToRASMatrixMoving);
  registrator->SetMovingIJKToXYZ(IJKToRASMatrixMoving);
  IJKToRASMatrixMoving->Delete();

  registrator->SetImageToImageMetricToMutualInformation();
  registrator->SetMetricComputationSamplingRatio(0.3333);
  registrator->SetNumberOfIterations(5);
  // registrator->SetNumberOfLevels(2); 

  registrator->SetTransformInitializationTypeToImageCenters();
  registrator->SetIntensityInterpolationTypeToLinear();

  try
    {
    //
    // run registration
    registrator->RegisterImages();
    fixedRASToMovingRASTransform->DeepCopy(registrator->GetTransform());

    if (outputVolumeNode != NULL)
      {
      //
      // resample moving image
      vtkTumorGrowthLogic::LinearResample(movingVolumeNode, outputVolumeNode, fixedVolumeNode, fixedRASToMovingRASTransform, backgroundLevel);
      }
    }
  catch (...)
    {
    std::cerr << "Failed to register images!!!" << std::endl;
    }
    
  //
  // clean up
  registrator->Delete();
}

void vtkTumorGrowthLogic::LinearResample (vtkMRMLVolumeNode* inputVolumeNode, vtkMRMLVolumeNode* outputVolumeNode, vtkMRMLVolumeNode* outputVolumeGeometryNode,
                    vtkTransform* outputRASToInputRASTransform, double backgroundLevel)
{
  vtkImageData* inputImageData  = inputVolumeNode->GetImageData();
  vtkImageData* outputImageData = outputVolumeNode->GetImageData();
  vtkImageData* outputGeometryData = NULL;
  if (outputVolumeGeometryNode != NULL)
    {
    outputGeometryData = outputVolumeGeometryNode->GetImageData();
    }

  vtkImageReslice* resliceFilter = vtkImageReslice::New();

  //
  // set inputs
  resliceFilter->SetInput(inputImageData);

  //
  // set geometry
  if (outputGeometryData != NULL)
    {
    resliceFilter->SetInformationInput(outputGeometryData);
    outputVolumeNode->CopyOrientation(outputVolumeGeometryNode);
    }

  //
  // setup total transform
  // ijk of output -> RAS -> XFORM -> RAS -> ijk of input
  vtkTransform* totalTransform = vtkTransform::New();
  if (outputRASToInputRASTransform != NULL)
    {
    totalTransform->DeepCopy(outputRASToInputRASTransform);
    }

  vtkMatrix4x4* outputIJKToRAS  = vtkMatrix4x4::New();
  outputVolumeNode->GetIJKToRASMatrix(outputIJKToRAS);
  vtkMatrix4x4* inputRASToIJK = vtkMatrix4x4::New();
  inputVolumeNode->GetRASToIJKMatrix(inputRASToIJK);

  totalTransform->PreMultiply();
  totalTransform->Concatenate(outputIJKToRAS);
  totalTransform->PostMultiply();
  totalTransform->Concatenate(inputRASToIJK);
  resliceFilter->SetResliceTransform(totalTransform);

  //
  // resample the image
  resliceFilter->SetBackgroundLevel(backgroundLevel);
  resliceFilter->OptimizationOn();
  resliceFilter->SetInterpolationModeToLinear();
  resliceFilter->Update();
  outputImageData->ShallowCopy(resliceFilter->GetOutput());

  //
  // clean up
  outputIJKToRAS->Delete();
  inputRASToIJK->Delete();
  resliceFilter->Delete();
  totalTransform->Delete();
}
