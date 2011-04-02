#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>

#include "vtkObjectFactory.h"
#include "vtkImageChangeInformation.h"

#include "vtkEMSegmentLogic.h"
#include "vtkEMSegment.h"

#include "vtkMRMLScene.h"

#include "vtkMRMLEMSTemplateNode.h"
#include "vtkMRMLEMSTreeNode.h"
#include "vtkMRMLEMSTreeParametersLeafNode.h"
#include "vtkMRMLEMSTreeParametersParentNode.h"
#include "vtkMRMLEMSTreeParametersNode.h"
#include "vtkMRMLEMSWorkingDataNode.h"
#include "vtkImageEMLocalSegmenter.h"
#include "vtkImageEMLocalSuperClass.h"
#include "vtkMath.h"
#include "vtkImageReslice.h"
#include "vtkRigidRegistrator.h"
#include "vtkBSplineRegistrator.h"
#include "vtkTransformToGrid.h"
#include "vtkIdentityTransform.h"
#include "vtkSlicerApplication.h"
#include "../../Applications/GUI/Slicer3Helper.cxx"
#include "vtkKWTkUtilities.h"

#include "vtkMRMLEMSAtlasNode.h"
#include "vtkMRMLEMSGlobalParametersNode.h"

#include "vtkImageIslandFilter.h"
#include "vtkSlicerVolumesLogic.h"

// needed to translate between enums
#include "EMLocalInterface.h"

#include <math.h>
#include <exception>

#include <vtksys/SystemTools.hxx>
#include "vtkDirectory.h"
#include "vtkMatrix4x4.h"

#define ERROR_NODE_VTKID 0

// A helper class to compare two maps
template <class T>
class MapCompare
{
public:
  static bool 
  map_value_comparer(typename std::map<T, unsigned int>::value_type &i1, 
                     typename std::map<T, unsigned int>::value_type &i2)
  {
  return i1.second<i2.second;
  }
};

//----------------------------------------------------------------------------
vtkEMSegmentLogic* vtkEMSegmentLogic::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = 
    vtkObjectFactory::CreateInstance("vtkEMSegmentLogic");
  if(ret)
    {
    return (vtkEMSegmentLogic*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkEMSegmentLogic;
}


//----------------------------------------------------------------------------
vtkEMSegmentLogic::vtkEMSegmentLogic()
{
  this->ModuleName = NULL;

  this->ProgressCurrentAction = NULL;
  this->ProgressGlobalFractionCompleted = 0.0;
  this->ProgressCurrentFractionCompleted = 0.0;

  //this->DebugOn();

  this->MRMLManager = NULL; // NB: must be set before SetMRMLManager is called
  vtkEMSegmentMRMLManager* manager = vtkEMSegmentMRMLManager::New();
  this->SetMRMLManager(manager);
  manager->Delete();
}

//----------------------------------------------------------------------------
vtkEMSegmentLogic::~vtkEMSegmentLogic()
{
  this->SetMRMLManager(NULL);
  this->SetProgressCurrentAction(NULL);
  this->SetModuleName(NULL);
}

//----------------------------------------------------------------------------
void vtkEMSegmentLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  // !!! todo
}


vtkMRMLScalarVolumeNode* 
vtkEMSegmentLogic::AddArchetypeScalarVolume (const char* filename, const char* volname, vtkSlicerApplicationLogic* appLogic,  vtkMRMLScene* mrmlScene)
{
  vtkSlicerVolumesLogic* volLogic  = vtkSlicerVolumesLogic::New();
  volLogic->SetMRMLScene(mrmlScene);
  volLogic->SetApplicationLogic(appLogic);
  vtkMRMLScalarVolumeNode* volNode = volLogic->AddArchetypeScalarVolume(filename, volname,2);
  volLogic->Delete();
  return  volNode;
}


//----------------------------------------------------------------------------
bool
vtkEMSegmentLogic::SaveIntermediateResults(vtkSlicerApplication* app, vtkSlicerApplicationLogic *appLogic)
{
  //
  // get output directory
  std::string outputDirectory(this->MRMLManager->GetSaveWorkingDirectory());

  if (!vtksys::SystemTools::FileExists(outputDirectory.c_str()))
    {
       // try to create directory
       bool createdOK = true;
       createdOK = vtksys::SystemTools::MakeDirectory(outputDirectory.c_str());
       if (!createdOK) {
              std::string  msg = "SaveIntermediateResults: could not create " + outputDirectory  + "!" ;
              ErrorMsg += msg + "\n";
              vtkErrorMacro(<< msg);
              return false;
       }
    }

  // check again whether or not directory exists
  if (!vtksys::SystemTools::FileExists(outputDirectory.c_str()))
    {
      std::string  msg = "SaveIntermediateResults: Directory " + outputDirectory  + " does not exist !" ;
      ErrorMsg += msg + "\n"; 
      vtkErrorMacro(<< msg);
      return false;
    }  

  //
  // package EMSeg-related parameters together and write them to disk
  bool writeSuccessful = this->PackageAndWriteData(app,appLogic,outputDirectory.c_str());

  return writeSuccessful;
}

//----------------------------------------------------------------------------
// New Task Specific Pipeline
//----------------------------------------------------------------------------

int vtkEMSegmentLogic::SourceTclFile(vtkSlicerApplication*app,const char *tclFile)
{
  // Load Tcl File defining the setting
  if (!app->LoadScript(tclFile))
    {
      vtkErrorMacro("Could not load in data for task. The following file does not exist: " << tclFile);
      return 1;
    }
  return 0 ;
}

//----------------------------------------------------------------------------

int vtkEMSegmentLogic::SourceTaskFiles(vtkSlicerApplication* app) { 
  vtkstd::string generalFile = this->DefineTclTaskFullPathName(app, vtkMRMLEMSGlobalParametersNode::GetDefaultTaskTclFileName());
  vtkstd::string specificFile = this->DefineTclTaskFileFromMRML(app);
  cout << "Sourcing general Task file : " << generalFile.c_str() << endl;
  // Have to first source the default file to set up the basic structure"
  if (this->SourceTclFile(app,generalFile.c_str()))
    {
      return 1;
    }
  // Now we overwrite anything from the default
  if (specificFile.compare(generalFile))
    {
      cout << "Sourcing task specific file: " <<   specificFile << endl;
      return this->SourceTclFile(app,specificFile.c_str()); 
    }
  return 0;
}

//----------------------------------------------------------------------------  
int vtkEMSegmentLogic::SourcePreprocessingTclFiles(vtkSlicerApplication* app) 
{
  if (this->SourceTaskFiles(app))
    {
      return 1;
    }
   // Source all files here as we otherwise sometimes do not find the function as Tcl did not finish sourcing but our cxx file is already trying to call the function 
   vtkstd::string tclFile =  this->GetModuleShareDirectory();
#ifdef _WIN32
   tclFile.append("\\Tcl\\EMSegmentAutoSample.tcl");
#else
   tclFile.append("/Tcl/EMSegmentAutoSample.tcl");
#endif
   return this->SourceTclFile(app,tclFile.c_str());
}

//----------------------------------------------------------------------------
bool
vtkEMSegmentLogic::
StartPreprocessingInitializeInputData()
{
  this->MRMLManager->GetWorkingDataNode()->SetInputTargetNodeIsValid(1);
  this->MRMLManager->GetWorkingDataNode()->SetInputAtlasNodeIsValid(1);
  this->MRMLManager->GetWorkingDataNode()->SetAlignedTargetNodeIsValid(0);
  this->MRMLManager->GetWorkingDataNode()->SetAlignedAtlasNodeIsValid(0);

  return true;
}

void
vtkEMSegmentLogic::
PrintImageInfo(vtkMRMLVolumeNode* volumeNode)
{
  if (volumeNode == NULL || volumeNode->GetImageData() == NULL)
    {
    std::cout << "Volume node or image data is null" << std::endl;
    return;
    }

  // extent
  int extent[6];
  volumeNode->GetImageData()->GetExtent(extent);
  std::cout << "Extent: " << std::endl;
  std::cout  << extent[0] << " " << extent[1] << " " << extent[2] << " " << extent[3] << " " << extent[4] << " " << extent[5] << std::endl;

  // ijkToRAS
  vtkMatrix4x4* matrix = vtkMatrix4x4::New();
  volumeNode->GetIJKToRASMatrix(matrix);
  std::cout << "IJKtoRAS Matrix: " << std::endl;
  for (unsigned int r = 0; r < 4; ++r)
    {
    std::cout << "   ";
    for (unsigned int c = 0; c < 4; ++c)
      {
      std::cout 
        << matrix->GetElement(r,c)
        << "   ";
      }
    std::cout << std::endl;
    }  
  matrix->Delete();
}

// a utility to print out a vtk image origin, spacing, and extent
//----------------------------------------------------------------------------
void
vtkEMSegmentLogic::
PrintImageInfo(vtkImageData* image)
{
  double spacing[3];
  double origin[3];
  int extent[6];

  image->GetSpacing(spacing);
  image->GetOrigin(origin);
  image->GetExtent(extent);

  std::cout << "Spacing: " << spacing[0] << " " << spacing[1] << " " << spacing[2] << std::endl;
  std::cout << "Origin: " << origin[0] << " " << origin[1] << " " << origin[2] << std::endl;
  std::cout << "Extent: " << extent[0] << " " << extent[1] << " " << extent[2] << " " << extent[3] << " " << extent[4] << " " << extent[5] << std::endl;
}

bool 
vtkEMSegmentLogic::
IsVolumeGeometryEqual(vtkMRMLVolumeNode* lhs,
                      vtkMRMLVolumeNode* rhs)
{
  if (lhs == NULL || rhs == NULL ||
      lhs->GetImageData() == NULL || rhs->GetImageData() == NULL)
    {
    return false;
    }

  // check extent
  int extentLHS[6];
  lhs->GetImageData()->GetExtent(extentLHS);
  int extentRHS[6];
  rhs->GetImageData()->GetExtent(extentRHS);
  bool equalExent = std::equal(extentLHS, extentLHS+6, extentRHS);
  
  // check ijkToRAS
  vtkMatrix4x4* matrixLHS = vtkMatrix4x4::New();
  lhs->GetIJKToRASMatrix(matrixLHS);
  vtkMatrix4x4* matrixRHS = vtkMatrix4x4::New();
  rhs->GetIJKToRASMatrix(matrixRHS);  
  bool equalMatrix = true;
  for (int r = 0; r < 4; ++r)
    {
    for (int c = 0; c < 4; ++c)
      {
        // Otherwise small errors will cause that they are not equal but should be ignored !
    if (double(int((*matrixLHS)[r][c]*100000)/100000.0) != double(int((*matrixRHS)[r][c]*100000)/100000.0))
        {
        equalMatrix = false;
    break;
        }
      }
    }

  matrixLHS->Delete();
  matrixRHS->Delete();
  return equalExent && equalMatrix;
}

// loops through the faces of the image bounding box and counts all the different image values and stores them in a map
// T represents the image data type
template <class T>
T
vtkEMSegmentLogic::
GuessRegistrationBackgroundLevel(vtkImageData* imageData)
{
  int borderWidth = 5;
  T inLevel;
  typedef std::map<T, unsigned int> MapType;
  MapType m;
  long totalVoxelsCounted = 0;

  T* inData = static_cast<T*>(imageData->GetScalarPointer());
  int dim[3];
  imageData->GetDimensions(dim);

  vtkIdType inc[3];
  vtkIdType iInc, jInc, kInc;
  imageData->GetIncrements(inc);

   // k first slice
  for (int k = 0; k < borderWidth; ++k)
    {
    kInc = k*inc[2];
    for (int j = 0; j < dim[1]; ++j)
      {
      jInc = j*inc[1];
      for (int i = 0; i < dim[0]; ++i)
        {
        iInc = i*inc[0];
        inLevel = inData[iInc+jInc+kInc];
        if (m.count(inLevel))
          {
          ++m[inLevel];
          }
        else
          {
          m[inLevel] = 1;
          }
        ++totalVoxelsCounted;
        }
      }
    }

  // k last slice
  for (int k = dim[2]-borderWidth; k < dim[2]; ++k)
    {
    kInc = k*inc[2];
    for (int j = 0; j < dim[1]; ++j)
      {
      jInc = j*inc[1];
      for (int i = 0; i < dim[0]; ++i)
        {
        iInc = i*inc[0];
        inLevel = inData[iInc+jInc+kInc];
        if (m.count(inLevel))
          {
          ++m[inLevel];
          }
        else
          {
          m[inLevel] = 1;
          }
        ++totalVoxelsCounted;
        }
      }
    }

  // j first slice
  for (int j = 0; j < borderWidth; ++j)
    {
    jInc = j*inc[1];
    for (int k = 0; k < dim[2]; ++k)
      {
      kInc = k*inc[2];
      for (int i = 0; i < dim[0]; ++i)
        {
        iInc = i*inc[0];
        inLevel = inData[iInc+jInc+kInc];
        if (m.count(inLevel))
          {
          ++m[inLevel];
          }
        else
          {
          m[inLevel] = 1;
          }
        ++totalVoxelsCounted;
        }
      }
    }

  // j last slice
  for (int j = dim[1]-borderWidth; j < dim[1]; ++j)
    {
    jInc = j*inc[1];
    for (int k = 0; k < dim[2]; ++k)
      {
      kInc = k*inc[2];
      for (int i = 0; i < dim[0]; ++i)
        {
        iInc = i*inc[0];
        inLevel = inData[iInc+jInc+kInc];
        if (m.count(inLevel))
          {
          ++m[inLevel];
          }
        else
          {
          m[inLevel] = 1;
          }
        ++totalVoxelsCounted;
        }
      }
    }

  // i first slice
  for (int i = 0; i < borderWidth; ++i)
    {
    iInc = i*inc[0];
    for (int k = 0; k < dim[2]; ++k)
      {
      kInc = k*inc[2];
      for (int j = 0; j < dim[1]; ++j)
        {
        jInc = j*inc[1];
        inLevel = inData[iInc+jInc+kInc];
        if (m.count(inLevel))
          {
          ++m[inLevel];
          }
        else
          {
          m[inLevel] = 1;
          }
        ++totalVoxelsCounted;
        }
      }
    }

  // i last slice
  for (int i = dim[0]-borderWidth; i < dim[0]; ++i)
    {
    iInc = i*inc[0];
    for (int k = 0; k < dim[2]; ++k)
      {
      kInc = k*inc[2];
      for (int j = 0; j < dim[1]; ++j)
        {
        jInc = j*inc[1];
        inLevel = inData[iInc+jInc+kInc];
        if (m.count(inLevel))
          {
          ++m[inLevel];
          }
        else
          {
          m[inLevel] = 1;
          }
        ++totalVoxelsCounted;
        }
      }
    }

  // all the information is stored in map m :  std::map<T, unsigned int>

  if (m.empty())
    {
    // no image data provided?
    return 0;
    }
  else if (m.size() == 1)
    {
      // Homogeneous background
      return m.begin()->first;
   }
  else
    {
    // search for the largest element
    typename MapType::iterator itor = 
      std::max_element(m.begin(), m.end(),
                       MapCompare<T>::map_value_comparer);

    // the iterator is pointing to the element with the largest value in the range [m.begin(), m.end()]
    T backgroundLevel = itor->first;

    // how many counts?
    double percentageOfVoxels = 
      100.0 * static_cast<double>(itor->second)/totalVoxelsCounted;

    std::cout << "   Background level guess : "<< std::endl
              << "   first place: "
              << static_cast<int>(backgroundLevel) << " (" << percentageOfVoxels << "%) "
              << std::endl;


    // erase largest element
    m.erase(itor);


    // again, search for the largest element (second place)
    typename MapType::iterator itor2 = 
      std::max_element(m.begin(), m.end(),
                       MapCompare<T>::map_value_comparer);

    T backgroundLevel_second_place = itor2->first;

    double percentageOfVoxels_secondplace =
      100.0 * static_cast<double>(itor2->second)/totalVoxelsCounted;

    std::cout << "   second place: "
              << static_cast<int>(backgroundLevel_second_place) << " (" << percentageOfVoxels_secondplace << "%)"
              << std::endl;

    return backgroundLevel;
    }
}

//
// A Slicer3 wrapper around vtkImageReslice.  Reslice the image data
// from inputVolumeNode into outputVolumeNode with the output image
// geometry specified by outputVolumeGeometryNode.  Optionally specify
// a transform.  The reslice transform will be:
//
// outputIJK->outputRAS->(outputRASToInputRASTransform)->inputRAS->inputIJK
//
//----------------------------------------------------------------------------
void
vtkEMSegmentLogic::
SlicerImageReslice(vtkMRMLVolumeNode* inputVolumeNode,
                   vtkMRMLVolumeNode* outputVolumeNode,
                   vtkMRMLVolumeNode* outputVolumeGeometryNode,
                   vtkTransform* outputRASToInputRASTransform,
                   int interpolationType,
                   double backgroundLevel)
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

  switch (interpolationType)
    {
    case vtkEMSegmentMRMLManager::InterpolationNearestNeighbor:
      resliceFilter->SetInterpolationModeToNearestNeighbor();
      break;
    case vtkEMSegmentMRMLManager::InterpolationCubic:
      resliceFilter->SetInterpolationModeToCubic();
      break;
    case vtkEMSegmentMRMLManager::InterpolationLinear:
    default:
      resliceFilter->SetInterpolationModeToLinear();
    }

  resliceFilter->Update();
  outputImageData->ShallowCopy(resliceFilter->GetOutput());

  //
  // clean up
  outputIJKToRAS->Delete();
  inputRASToIJK->Delete();
  resliceFilter->Delete();
  totalTransform->Delete();
}

// Assume geometry is already specified, create
// outGrid(p) = postMultiply \circ inGrid \circ preMultiply (p)
//
// right now simplicity over speed.  Optimize later?
//----------------------------------------------------------------------------
void
vtkEMSegmentLogic::
ComposeGridTransform(vtkGridTransform* inGrid,
                     vtkMatrix4x4*     preMultiply,
                     vtkMatrix4x4*     postMultiply,
                     vtkGridTransform* outGrid)
{
  // iterate over output grid
  double inPt[4] = {0, 0, 0, 1};
  double pt[4]   = {0, 0, 0, 1};
  double* outDataPtr = 
    static_cast<double*>(outGrid->GetDisplacementGrid()->GetScalarPointer());  
  vtkIdType numOutputVoxels = outGrid->GetDisplacementGrid()->
    GetNumberOfPoints();

  for (vtkIdType i = 0; i < numOutputVoxels; ++i)
    {
    outGrid->GetDisplacementGrid()->GetPoint(i, inPt);
    preMultiply->MultiplyPoint(inPt, pt);
    inGrid->TransformPoint(pt, pt);
    postMultiply->MultiplyPoint(pt, pt);
    
    *outDataPtr++ = pt[0] - inPt[0];
    *outDataPtr++ = pt[1] - inPt[1];
    *outDataPtr++ = pt[2] - inPt[2];
    }
}

//
// A Slicer3 wrapper around vtkImageReslice.  Reslice the image data
// from inputVolumeNode into outputVolumeNode with the output image
// geometry specified by outputVolumeGeometryNode.  Optionally specify
// a transform.  The reslice transorm will be:
//
// outputIJK->outputRAS->(outputRASToInputRASTransform)->inputRAS->inputIJK
//
//----------------------------------------------------------------------------
void
vtkEMSegmentLogic::
SlicerImageResliceWithGrid(vtkMRMLVolumeNode* inputVolumeNode,
                           vtkMRMLVolumeNode* outputVolumeNode,
                           vtkMRMLVolumeNode* outputVolumeGeometryNode,
                           vtkGridTransform* outputRASToInputRASTransform,
                           int interpolationType,
                           double backgroundLevel)
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
  // create total transform
  vtkTransformToGrid* gridSource = vtkTransformToGrid::New();
  vtkIdentityTransform* idTransform = vtkIdentityTransform::New();
  gridSource->SetInput(idTransform);
  //gridSource->SetGridScalarType(VTK_FLOAT);
  idTransform->Delete();

  //
  // set geometry
  if (outputGeometryData != NULL)
    {
    resliceFilter->SetInformationInput(outputGeometryData);
    outputVolumeNode->CopyOrientation(outputVolumeGeometryNode);

    gridSource->SetGridExtent(outputGeometryData->GetExtent());
    gridSource->SetGridSpacing(outputGeometryData->GetSpacing());
    gridSource->SetGridOrigin(outputGeometryData->GetOrigin());
    }
  else
    {
    gridSource->SetGridExtent(outputImageData->GetExtent());
    gridSource->SetGridSpacing(outputImageData->GetSpacing());
    gridSource->SetGridOrigin(outputImageData->GetOrigin());
    }
  gridSource->Update();
  vtkGridTransform* totalTransform = vtkGridTransform::New();
  totalTransform->SetDisplacementGrid(gridSource->GetOutput());
//  totalTransform->SetInterpolationModeToCubic();
  gridSource->Delete();
  
  //
  // fill in total transform
  // ijk of output -> RAS -> XFORM -> RAS -> ijk of input
  vtkMatrix4x4* outputIJKToRAS  = vtkMatrix4x4::New();
  outputVolumeNode->GetIJKToRASMatrix(outputIJKToRAS);
  vtkMatrix4x4* inputRASToIJK = vtkMatrix4x4::New();
  inputVolumeNode->GetRASToIJKMatrix(inputRASToIJK);
  vtkEMSegmentLogic::ComposeGridTransform(outputRASToInputRASTransform,
                                          outputIJKToRAS,
                                          inputRASToIJK,
                                          totalTransform);
  resliceFilter->SetResliceTransform(totalTransform);

  //
  // resample the image
  resliceFilter->SetBackgroundLevel(backgroundLevel);
  resliceFilter->OptimizationOn();

  switch (interpolationType)
    {
    case vtkEMSegmentMRMLManager::InterpolationNearestNeighbor:
      resliceFilter->SetInterpolationModeToNearestNeighbor();
      break;
    case vtkEMSegmentMRMLManager::InterpolationCubic:
      resliceFilter->SetInterpolationModeToCubic();
      break;
    case vtkEMSegmentMRMLManager::InterpolationLinear:
    default:
      resliceFilter->SetInterpolationModeToLinear();
    }

  resliceFilter->Update();
  outputImageData->ShallowCopy(resliceFilter->GetOutput());

  //
  // clean up
  outputIJKToRAS->Delete();
  inputRASToIJK->Delete();
  resliceFilter->Delete();
  totalTransform->Delete();
}


void vtkEMSegmentLogic::StartPreprocessingResampleAndCastToTarget(vtkMRMLVolumeNode* movingVolumeNode, vtkMRMLVolumeNode* fixedVolumeNode, vtkMRMLVolumeNode* outputVolumeNode)
{
  if (!vtkEMSegmentLogic::IsVolumeGeometryEqual(fixedVolumeNode, outputVolumeNode))
    {

      std::cout << "Warning: Target-to-target registration skipped but "
                << "target images have differenent geometries. "
                << std::endl
                << "Suggestion: If you are not positive that your images are "
                << "aligned, you should enable target-to-target registration."
                << std::endl;

      std::cout << "Fixed Volume Node: " << std::endl;
      PrintImageInfo(fixedVolumeNode);
      std::cout << "Output Volume Node: " << std::endl;
      PrintImageInfo(outputVolumeNode);

      // std::cout << "Resampling target image " << i << "...";
      double backgroundLevel = 0;
      switch (movingVolumeNode->GetImageData()->GetScalarType())
        {  
          vtkTemplateMacro(backgroundLevel = (GuessRegistrationBackgroundLevel<VTK_TT>(movingVolumeNode->GetImageData())););
        }
      std::cout << "   Guessed background level: " << backgroundLevel << std::endl;

      vtkEMSegmentLogic::SlicerImageReslice(movingVolumeNode, 
                                            outputVolumeNode, 
                                            fixedVolumeNode,
                                            NULL,
                                            vtkEMSegmentMRMLManager::InterpolationLinear,
                                            backgroundLevel);
    }

  if (fixedVolumeNode->GetImageData()->GetScalarType() != movingVolumeNode->GetImageData()->GetScalarType())
    {
      //cast
      vtkImageCast* cast = vtkImageCast::New();
      cast->SetInput(outputVolumeNode->GetImageData());
      cast->SetOutputScalarType(fixedVolumeNode->GetImageData()->GetScalarType());
      cast->Update();
      outputVolumeNode->GetImageData()->DeepCopy(cast->GetOutput());
      cast->Delete();
    }
  std::cout << "Resampling and casting output volume \"" << outputVolumeNode->GetName() << "\" to reference target \"" << fixedVolumeNode->GetName() <<  "\" DONE" << std::endl;
}

//----------------------------------------------------------------------------
double vtkEMSegmentLogic::GuessRegistrationBackgroundLevel(vtkMRMLVolumeNode* volumeNode)
{
  if (!volumeNode ||  !volumeNode->GetImageData())  
    {
      std::cerr << "double vtkEMSegmentLogic::GuessRegistrationBackgroundLevel(vtkMRMLVolumeNode* volumeNode) : volumeNode or volumeNode->GetImageData is null" << std::endl;
      return -1;
    }

  // guess background level    
  double backgroundLevel = 0;
  switch (volumeNode->GetImageData()->GetScalarType())
      {  
        vtkTemplateMacro(backgroundLevel = (GuessRegistrationBackgroundLevel<VTK_TT>(volumeNode->GetImageData())););
      }
  std::cout << "   Guessed background level: " << backgroundLevel << std::endl;
  return backgroundLevel;
}

//----------------------------------------------------------------------------
int vtkEMSegmentLogic::StartSegmentationWithoutPreprocessing(vtkSlicerApplication* app, vtkSlicerApplicationLogic *appLogic)
{
  //
  // make sure we're ready to start
  //
  ErrorMsg.clear();

  if (!this->MRMLManager->GetWorkingDataNode()->GetAlignedTargetNodeIsValid() ||
      !this->MRMLManager->GetWorkingDataNode()->GetAlignedAtlasNodeIsValid())
    {
    ErrorMsg = "Preprocessing pipeline not up to date!  Aborting Segmentation.";
    vtkErrorMacro( << ErrorMsg );
    return EXIT_FAILURE;
    }


  // find output volume
  if (!this->MRMLManager->GetNode())
    {
    ErrorMsg     = "Template node is null---aborting segmentation.";
    vtkErrorMacro( << ErrorMsg );
    return EXIT_FAILURE;
    }
  vtkMRMLScalarVolumeNode *outVolume = this->MRMLManager->GetOutputVolumeNode();
  if (outVolume == NULL)
    {
    ErrorMsg     = "No output volume found---aborting segmentation.";
    vtkErrorMacro( << ErrorMsg );
    return EXIT_FAILURE;
    }

  //
  // Copy RASToIJK matrix, and other attributes from input to
  // output. Use first target volume as source for this data.
  //
  
  // get attributes from first target input volume
  const char* inMRLMID = 
    this->MRMLManager->GetTargetInputNode()->GetNthVolumeNodeID(0);
  vtkMRMLScalarVolumeNode *inVolume = vtkMRMLScalarVolumeNode::
    SafeDownCast(this->GetMRMLScene()->GetNodeByID(inMRLMID));
  if (inVolume == NULL)
    {
    ErrorMsg     = "Can't get first target image.";
    vtkErrorMacro( << ErrorMsg); 
    return EXIT_FAILURE;
    }

  outVolume->CopyOrientation(inVolume);
  outVolume->SetAndObserveTransformNodeID(inVolume->GetTransformNodeID());

  //
  // create segmenter class
  //
  vtkImageEMLocalSegmenter* segmenter = vtkImageEMLocalSegmenter::New();
  if (segmenter == NULL)
    {
    ErrorMsg = "Could not create vtkImageEMLocalSegmenter pointer";
    vtkErrorMacro( << ErrorMsg );
    return EXIT_FAILURE;
    }

  //
  // copy mrml data to segmenter class
  //
  vtkstd::cout << "EMSEG: Copying data to algorithm class...";
  this->CopyDataToSegmenter(segmenter);
  vtkstd::cout << "DONE" << vtkstd::endl;

  if (this->GetDebug())
  {
    vtkstd::cout << "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB" << vtkstd::endl;
    vtkstd::cout << "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB" << vtkstd::endl;
    vtkstd::cout << "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB" << vtkstd::endl;
    vtkIndent indent;
    segmenter->PrintSelf(vtkstd::cout, indent);
    vtkstd::cout << "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB" << vtkstd::endl;
    vtkstd::cout << "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB" << vtkstd::endl;
    vtkstd::cout << "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB" << vtkstd::endl;
  }

  //
  // start segmentation
  //
  try 
    {
    vtkstd::cout << "[Start] Segmentation algorithm..." << vtkstd::endl;
    segmenter->Update();
    vtkstd::cout << "[Done]  Segmentation algorithm." << vtkstd::endl;
    }
  catch (std::exception& e)
    {
    ErrorMsg = "Exception thrown during segmentation: "  + std::string(e.what()) + "\n";
    vtkErrorMacro( << ErrorMsg );
    return EXIT_FAILURE;
    } 

  if (this->GetDebug())
  {
    vtkstd::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" << vtkstd::endl;
    vtkstd::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" << vtkstd::endl;
    vtkstd::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" << vtkstd::endl;
    segmenter->PrintSelf(vtkstd::cout, static_cast<vtkIndent>(0));
    vtkstd::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" << vtkstd::endl;
    vtkstd::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" << vtkstd::endl;
    vtkstd::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" << vtkstd::endl;
  }

  // POST PROCESSING 
  vtkstd::cout << "[Start] Postprocessing ..." << vtkstd::endl;
  vtkImageData* postProcessing = vtkImageData::New(); 
  postProcessing->ShallowCopy(segmenter->GetOutput());

  // Subparcellation
  if (this->MRMLManager->GetEnableSubParcellation()) {
    vtkstd::cout << "=== Sub-Parcellation === " << vtkstd::endl;
    this->SubParcelateSegmentation(postProcessing,this->MRMLManager->GetTreeRootNodeID()); 
  }

  // Island Removal
  if (this->MRMLManager->GetMinimumIslandSize() > 1) {
     vtkstd::cout << "=== Island removal === " << vtkstd::endl;
     vtkImageData* input = vtkImageData::New();
     input->DeepCopy(postProcessing);
     vtkImageIslandFilter* islandFilter = vtkImageIslandFilter::New();
     islandFilter->SetInput(input);
     islandFilter->SetIslandMinSize(this->MRMLManager->GetMinimumIslandSize());
     islandFilter->SetNeighborhoodDim3D();
     islandFilter->SetPrintInformation(1);
     islandFilter->Update();
     postProcessing->DeepCopy(islandFilter->GetOutput());
     islandFilter->Delete();
     input->Delete();
  }
  vtkstd::cout << "[Done] Postprocessing" << vtkstd::endl;
  //
  // copy result to output volume
  //
  
  // set output of the filter to VolumeNode's ImageData

  outVolume->SetAndObserveImageData(postProcessing);
  postProcessing->Delete();
  // make sure the output volume is a labelmap
  if (!outVolume->GetLabelMap())
  {
    vtkWarningMacro("Changing output image to labelmap");
    outVolume->LabelMapOn();
  }

  vtkMRMLVolumeDisplayNode *outDisplayNode = vtkMRMLVolumeDisplayNode::SafeDownCast(outVolume->GetDisplayNode());
  if (!outDisplayNode) 
    {
       vtkWarningMacro("Did not define lookup table bc display node is not defined ");
    } 
  else 
    {
      const char* colorID = this->MRMLManager->GetColormap();
      if (colorID) 
    {
             outDisplayNode->SetAndObserveColorNodeID(colorID);
    }
    }

    

  outVolume->SetModifiedSinceRead(1);

  //
  // clean up
  //
  segmenter->Delete();

  //
  // save intermediate results
  if (this->MRMLManager->GetSaveIntermediateResults())
    {
    vtkstd::cout << "[Start] Saving intermediate results..." << vtkstd::endl;
    bool savedResults = this->SaveIntermediateResults(app,appLogic);
    vtkstd::cout << "[Done]  Saving intermediate results." << vtkstd::endl;
    if (!savedResults)
      {
    std::string msg = "Error writing intermediate results"; 
        ErrorMsg += msg + "\n";
        vtkErrorMacro( << msg);
        return EXIT_FAILURE;
      }
    }

  return EXIT_SUCCESS;
}

//-----------------------------------------------------------------------------
void
vtkEMSegmentLogic::
PopulateTestingData()
{
  vtkDebugMacro("Begin populating test data");

  //
  // add some nodes to the hierarchy
  //
  vtkDebugMacro("Setting parameters for root node");
  double color[3];
  vtkIdType rootNodeID         = this->MRMLManager->GetTreeRootNodeID();
  this->MRMLManager->SetTreeNodeName(rootNodeID, "Root");
  color[0] = 1.0; color[1] = 0.0; color[2] = 0.0;
  this->MRMLManager->SetTreeNodeColor(rootNodeID, color);
  this->MRMLManager->SetTreeNodeSpatialPriorWeight(rootNodeID, 0.5);
  this->MRMLManager->SetTreeNodeClassProbability(rootNodeID, 0.5);
  this->MRMLManager->SetTreeNodeAlpha(rootNodeID, 0.5);
  this->MRMLManager->SetTreeNodePrintWeight(rootNodeID, 1);
  this->MRMLManager->SetTreeNodeStoppingConditionEMType(rootNodeID, 1);
  this->MRMLManager->SetTreeNodeStoppingConditionEMIterations(rootNodeID, 15);
  this->MRMLManager->SetTreeNodeStoppingConditionEMValue(rootNodeID, 0.5);
  this->MRMLManager->SetTreeNodeStoppingConditionMFAType(rootNodeID, 2);
  this->MRMLManager->SetTreeNodeStoppingConditionMFAIterations(rootNodeID, 16);
  this->MRMLManager->SetTreeNodeStoppingConditionMFAValue(rootNodeID, 0.6);

  vtkDebugMacro("Setting parameters for background node");
  vtkIdType backgroundNodeID   = this->MRMLManager->AddTreeNode(rootNodeID);
  this->MRMLManager->SetTreeNodeName(backgroundNodeID, "Background");
  color[0] = 0.0; color[1] = 0.0; color[2] = 0.0;
  this->MRMLManager->SetTreeNodeColor(backgroundNodeID, color);
  this->MRMLManager->SetTreeNodeSpatialPriorWeight(backgroundNodeID, 0.4);
  this->MRMLManager->SetTreeNodeClassProbability(backgroundNodeID, 0.4);
  this->MRMLManager->SetTreeNodePrintWeight(backgroundNodeID, 1);

  vtkDebugMacro("Setting parameters for icc node");
  vtkIdType iccNodeID          = this->MRMLManager->AddTreeNode(rootNodeID);
  this->MRMLManager->SetTreeNodeName(iccNodeID, "ICC");
  color[0] = 0.0; color[1] = 1.0; color[2] = 0.0;
  this->MRMLManager->SetTreeNodeColor(iccNodeID, color);
  this->MRMLManager->SetTreeNodeSpatialPriorWeight(iccNodeID, 0.3);
  this->MRMLManager->SetTreeNodeClassProbability(iccNodeID, 0.3);
  this->MRMLManager->SetTreeNodeAlpha(iccNodeID, 0.3);
  this->MRMLManager->SetTreeNodePrintWeight(iccNodeID, 1);
  this->MRMLManager->SetTreeNodeStoppingConditionEMType(iccNodeID, 0);
  this->MRMLManager->SetTreeNodeStoppingConditionEMIterations(iccNodeID, 13);
  this->MRMLManager->SetTreeNodeStoppingConditionEMValue(iccNodeID, 0.3);
  this->MRMLManager->SetTreeNodeStoppingConditionMFAType(iccNodeID, 1);
  this->MRMLManager->SetTreeNodeStoppingConditionMFAIterations(iccNodeID, 14);
  this->MRMLManager->SetTreeNodeStoppingConditionMFAValue(iccNodeID, 0.4);

  vtkDebugMacro("Setting parameters for grey matter node");
  vtkIdType greyMatterNodeID   = this->MRMLManager->AddTreeNode(iccNodeID);
  this->MRMLManager->SetTreeNodeName(greyMatterNodeID, "Grey Matter");
  color[0] = 0.0; color[1] = 1.0; color[2] = 1.0;
  this->MRMLManager->SetTreeNodeColor(greyMatterNodeID, color);
  this->MRMLManager->SetTreeNodeSpatialPriorWeight(greyMatterNodeID, 0.2);
  this->MRMLManager->SetTreeNodeClassProbability(greyMatterNodeID, 0.2);
  this->MRMLManager->SetTreeNodePrintWeight(greyMatterNodeID, 1);

  vtkDebugMacro("Setting parameters for white matter node");
  vtkIdType whiteMatterNodeID  = this->MRMLManager->AddTreeNode(iccNodeID);
  this->MRMLManager->SetTreeNodeName(whiteMatterNodeID, "White Matter");
  color[0] = 1.0; color[1] = 1.0; color[2] = 0.0;
  this->MRMLManager->SetTreeNodeColor(whiteMatterNodeID, color);
  this->MRMLManager->SetTreeNodeSpatialPriorWeight(whiteMatterNodeID, 0.1);
  this->MRMLManager->SetTreeNodeClassProbability(whiteMatterNodeID, 0.1);
  this->MRMLManager->SetTreeNodePrintWeight(whiteMatterNodeID, 1);

  vtkDebugMacro("Setting parameters for csf node");
  vtkIdType csfNodeID  = this->MRMLManager->AddTreeNode(iccNodeID);
  this->MRMLManager->SetTreeNodeName(csfNodeID, "CSF");

  //
  // set registration parameters
  //
  vtkDebugMacro("Setting registration parameters");
  this->MRMLManager->SetRegistrationAffineType(0);
  this->MRMLManager->SetRegistrationDeformableType(0);
  this->MRMLManager->SetRegistrationInterpolationType(1);

  //
  // set save parameters
  //
  vtkDebugMacro("Setting save parameters");
  this->MRMLManager->SetSaveWorkingDirectory("/tmp");
  this->MRMLManager->SetSaveTemplateFilename("/tmp/EMSTemplate.mrml");
  this->MRMLManager->SetSaveTemplateAfterSegmentation(1);
  this->MRMLManager->SetSaveIntermediateResults(1);
  this->MRMLManager->SetSaveSurfaceModels(1);
  
  this->MRMLManager->SetEnableMultithreading(1);
  this->SetProgressGlobalFractionCompleted(0.9);

  vtkDebugMacro("Done populating test data");
}

//-----------------------------------------------------------------------------
void
vtkEMSegmentLogic::
SpecialTestingFunction()
{
}

//-----------------------------------------------------------------------------
vtkIntArray*
vtkEMSegmentLogic::
NewObservableEvents()
{
  vtkIntArray *events = vtkIntArray::New();
  events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
  events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);

  return events;
}

//-----------------------------------------------------------------------------
void
vtkEMSegmentLogic::
CopyDataToSegmenter(vtkImageEMLocalSegmenter* segmenter)
{
  //
  // copy atlas related parameters to algorithm
  //
  vtkstd::cout << "atlas data...";
  this->CopyAtlasDataToSegmenter(segmenter);

  //
  // copy target related parameters to algorithm
  //
  vtkstd::cout << "target data...";
  this->CopyTargetDataToSegmenter(segmenter);

  //
  // copy global parameters to algorithm 
  //
  vtkstd::cout << "global data...";
  this->CopyGlobalDataToSegmenter(segmenter);

  //
  // copy tree base parameters to algorithm
  //
  vtkstd::cout << "tree data...";
  vtkImageEMLocalSuperClass* rootNode = vtkImageEMLocalSuperClass::New();
  this->CopyTreeDataToSegmenter(rootNode, 
                                this->MRMLManager->GetTreeRootNodeID());
  segmenter->SetHeadClass(rootNode);
  rootNode->Delete();
}

//-----------------------------------------------------------------------------
void
vtkEMSegmentLogic::
CopyAtlasDataToSegmenter(vtkImageEMLocalSegmenter* segmenter)
{
  segmenter->
    SetNumberOfTrainingSamples(this->MRMLManager->
                               GetAtlasNumberOfTrainingSamples());
}

//-----------------------------------------------------------------------------
void
vtkEMSegmentLogic::
CopyTargetDataToSegmenter(vtkImageEMLocalSegmenter* segmenter)
{
  // !!! todo: TESTING HERE!!!
  vtkMRMLEMSVolumeCollectionNode* workingTarget = 
    this->MRMLManager->GetWorkingDataNode()->GetAlignedTargetNode();
  unsigned int numTargetImages = workingTarget->GetNumberOfVolumes();
  std::cout << "Setting number of target images: " << numTargetImages 
            << std::endl;
  segmenter->SetNumInputImages(numTargetImages);

  for (unsigned int i = 0; i < numTargetImages; ++i)
    {
    std::string mrmlID = workingTarget->GetNthVolumeNodeID(i);
    vtkDebugMacro("Setting target image " << i << " mrmlID=" 
                  << mrmlID.c_str());

    vtkImageData* imageData = 
      workingTarget->GetNthVolumeNode(i)->GetImageData();

    std::cout << "AddingTargetImage..." << std::endl;
    this->PrintImageInfo(imageData);
    imageData->Update();
    this->PrintImageInfo(imageData);

    segmenter->SetImageInput(i, imageData);
    }
}

//-----------------------------------------------------------------------------
void
vtkEMSegmentLogic::
CopyGlobalDataToSegmenter(vtkImageEMLocalSegmenter* segmenter)
{
  if (this->MRMLManager->GetEnableMultithreading())
    {
    segmenter->
      SetDisableMultiThreading(0);
    }
  else
    {
    segmenter->
      SetDisableMultiThreading(1);
    }
  segmenter->SetPrintDir(this->MRMLManager->GetSaveWorkingDirectory());
  
  //
  // NB: In the algorithm code smoothing widht and sigma parameters
  // are defined globally.  In this logic, they are defined for each
  // parent node.  For now copy parameters from the root tree
  // node. !!!todo!!!
  //
  vtkIdType rootNodeID = this->MRMLManager->GetTreeRootNodeID();
  segmenter->
    SetSmoothingWidth(this->MRMLManager->
                      GetTreeNodeSmoothingKernelWidth(rootNodeID));

  // type mismatch between logic and algorithm !!!todo!!!
  int intSigma = 
    vtkMath::Round(this->MRMLManager->
                   GetTreeNodeSmoothingKernelSigma(rootNodeID));
  segmenter->SetSmoothingSigma(intSigma);

  //
  // registration parameters
  //
  int algType = this->ConvertGUIEnumToAlgorithmEnumInterpolationType
    (this->MRMLManager->GetRegistrationInterpolationType());
  segmenter->SetRegistrationInterpolationType(algType);
}

//-----------------------------------------------------------------------------
void
vtkEMSegmentLogic::
CopyTreeDataToSegmenter(vtkImageEMLocalSuperClass* node, vtkIdType nodeID)
{
  // need this here because the vtkImageEM* classes don't use
  // virtual functions and so failed initializations lead to
  // memory errors
  node->SetNumInputImages(this->MRMLManager->
                          GetTargetNumberOfSelectedVolumes());

  // copy generic tree node data to segmenter
  this->CopyTreeGenericDataToSegmenter(node, nodeID);
  
  // copy parent specific tree node data to segmenter
  this->CopyTreeParentDataToSegmenter(node, nodeID);

  // add children
  unsigned int numChildren = 
    this->MRMLManager->GetTreeNodeNumberOfChildren(nodeID);
  double totalProbability = 0.0;
  for (unsigned int i = 0; i < numChildren; ++i)
    {
    vtkIdType childID = this->MRMLManager->GetTreeNodeChildNodeID(nodeID, i);
    bool isLeaf = this->MRMLManager->GetTreeNodeIsLeaf(childID);

    if (isLeaf)
      {
      vtkImageEMLocalClass* childNode = vtkImageEMLocalClass::New();
      // need this here because the vtkImageEM* classes don't use
      // virtual functions and so failed initializations lead to
      // memory errors
      childNode->SetNumInputImages(this->MRMLManager->
                                   GetTargetNumberOfSelectedVolumes());
      this->CopyTreeGenericDataToSegmenter(childNode, childID);
      this->CopyTreeLeafDataToSegmenter(childNode, childID);
      node->AddSubClass(childNode, i);
      childNode->Delete();
      }
    else
      {
      vtkImageEMLocalSuperClass* childNode = vtkImageEMLocalSuperClass::New();
      this->CopyTreeDataToSegmenter(childNode, childID);
      node->AddSubClass(childNode, i);
      childNode->Delete();
      }

    totalProbability += 
      this->MRMLManager->GetTreeNodeClassProbability(childID);
    }

  if (totalProbability != 1.0)
    {
    vtkWarningMacro("Warning: child probabilities don't sum to unity for node "
                    << this->MRMLManager->GetTreeNodeName(nodeID)
                    << " they sum to " << totalProbability);
    }

  // Set Markov matrices
  const unsigned int numDirections = 6;
  for (unsigned int d = 0; d < numDirections; ++d)
    {
    for (unsigned int r = 0; r < numChildren; ++r)
      {
      for (unsigned int c = 0; c < numChildren; ++c)
        {
          double val = (r == c ? 1.0 : 0.0);
          node->SetMarkovMatrix(val, d, c, r);
        }
      }
    }
  node->Update();
}

//-----------------------------------------------------------------------------
void vtkEMSegmentLogic::DefineValidSegmentationBoundary() 
{
 //
  // Setup ROI.  If if looks bogus then use the default (entire image)
  bool useDefaultBoundary = false;
  int boundMin[3];
  int boundMax[3];

  // get dimensions of target image
  int targetImageDimensions[3];
  this->MRMLManager->GetTargetInputNode()->GetNthVolumeNode(0)->
    GetImageData()->GetDimensions(targetImageDimensions);

  this->MRMLManager->GetSegmentationBoundaryMin(boundMin);
  this->MRMLManager->GetSegmentationBoundaryMax(boundMax);
  // Specify boundary in 1-based, NOT 0-based as you might expect
  for (unsigned int i = 0; i < 3; ++i)
    {
    if (boundMin[i] <  1 || 
        boundMin[i] >  targetImageDimensions[i]   ||
        boundMax[i] <  1                   ||
        boundMax[i] >  targetImageDimensions[i]   ||
        boundMax[i] <  boundMin[i])
      {
      useDefaultBoundary = true;
      break;
      }
    }
  if (useDefaultBoundary)
    {
    std::cout 
      << std::endl
      << "====================================================================" << std::endl
      << "Warning: the segmentation ROI was bogus, setting ROI to entire image"  << std::endl
      << "Axis 0 -  Image Min: 1 <= RoiMin(" << boundMin[0] << ") <= ROIMax(" << boundMax[0] <<") <=  Image Max:" << targetImageDimensions[0] <<  std::endl
      << "Axis 1 -  Image Min: 1 <= RoiMin(" << boundMin[1] << ") <= ROIMax(" << boundMax[1] << ") <=  Image Max:" << targetImageDimensions[1] <<  std::endl
      << "Axis 2 -  Image Min: 1 <= RoiMin(" << boundMin[2] << ") <= ROIMax(" << boundMax[2] << ") <=  Image Max:" << targetImageDimensions[2] <<  std::endl
      << "NOTE: The above warning about ROI should not lead to poor segmentation results;  the entire image should be segmented.  It only indicates an error if you intended to segment a subregion of the image."
      << std::endl
      << "Define Boundary as: ";
      for (unsigned int i = 0; i < 3; ++i)
        {
          boundMin[i] = 1;
          boundMax[i] = targetImageDimensions[i];
          std::cout << boundMin[i] << ", " << boundMax[i] << ",   ";
        }
      std::cout << std::endl << "====================================================================" << std::endl;

      this->MRMLManager->SetSegmentationBoundaryMin(boundMin);
      this->MRMLManager->SetSegmentationBoundaryMax(boundMax); 
    }
}

void
vtkEMSegmentLogic::
CopyTreeGenericDataToSegmenter(vtkImageEMLocalGenericClass* node, 
                               vtkIdType nodeID)
{
  unsigned int numTargetImages = 
  this->MRMLManager->GetTargetNumberOfSelectedVolumes();

 
  this->DefineValidSegmentationBoundary();
  int boundMin[3];
  int boundMax[3];
  this->MRMLManager->GetSegmentationBoundaryMin(boundMin);
  this->MRMLManager->GetSegmentationBoundaryMax(boundMax);
  node->SetSegmentationBoundaryMin(boundMin[0], boundMin[1], boundMin[2]);
  node->SetSegmentationBoundaryMax(boundMax[0], boundMax[1], boundMax[2]);

  node->SetProbDataWeight(this->MRMLManager->
                          GetTreeNodeSpatialPriorWeight(nodeID));

  node->SetTissueProbability(this->MRMLManager->
                             GetTreeNodeClassProbability(nodeID));

  node->SetPrintWeights(this->MRMLManager->GetTreeNodePrintWeight(nodeID));

  // set target input channel weights
  for (unsigned int i = 0; i < numTargetImages; ++i)
    {
    node->SetInputChannelWeights(this->MRMLManager->
                                 GetTreeNodeInputChannelWeight(nodeID, 
                                                               i), i);
    }

  //
  // registration related data
  //
  //!!!bcd!!!

  //
  // set probability data
  //

  // get working atlas
  // !!! error checking!
  vtkMRMLVolumeNode*  atlasNode = this->MRMLManager->GetAlignedSpatialPriorFromTreeNodeID(nodeID);
  if (atlasNode)
    {
    vtkDebugMacro("Setting spatial prior: node=" 
                  << this->MRMLManager->GetTreeNodeName(nodeID));
    vtkImageData* imageData = atlasNode->GetImageData();
    node->SetProbDataPtr(imageData);
    }

  int exclude =  this->MRMLManager->GetTreeNodeExcludeFromIncompleteEStep(nodeID);
  node->SetExcludeFromIncompleteEStepFlag(exclude);
}


//-----------------------------------------------------------------------------
void
vtkEMSegmentLogic::
CopyTreeParentDataToSegmenter(vtkImageEMLocalSuperClass* node, 
                              vtkIdType nodeID)
{
  node->SetPrintFrequency (this->MRMLManager->
                           GetTreeNodePrintFrequency(nodeID));
  node->SetPrintBias      (this->MRMLManager->
                           GetTreeNodePrintBias(nodeID));
  node->SetPrintLabelMap  (this->MRMLManager->
                           GetTreeNodePrintLabelMap(nodeID));

  node->SetPrintEMLabelMapConvergence
    (this->MRMLManager->GetTreeNodePrintEMLabelMapConvergence(nodeID));
  node->SetPrintEMWeightsConvergence
    (this->MRMLManager->GetTreeNodePrintEMWeightsConvergence(nodeID));
  node->SetStopEMType(this->ConvertGUIEnumToAlgorithmEnumStoppingConditionType
                      (this->MRMLManager->
                      GetTreeNodeStoppingConditionEMType(nodeID)));
  node->SetStopEMValue(this->MRMLManager->
                       GetTreeNodeStoppingConditionEMValue(nodeID));
  node->SetStopEMMaxIter
    (this->MRMLManager->GetTreeNodeStoppingConditionEMIterations(nodeID));

  node->SetPrintMFALabelMapConvergence
    (this->MRMLManager->GetTreeNodePrintMFALabelMapConvergence(nodeID));
  node->SetPrintMFAWeightsConvergence
    (this->MRMLManager->GetTreeNodePrintMFAWeightsConvergence(nodeID));
  node->SetStopMFAType(this->ConvertGUIEnumToAlgorithmEnumStoppingConditionType
                       (this->MRMLManager->
                       GetTreeNodeStoppingConditionMFAType(nodeID)));
  node->SetStopMFAValue(this->MRMLManager->
                        GetTreeNodeStoppingConditionMFAValue(nodeID));
  node->SetStopMFAMaxIter
    (this->MRMLManager->GetTreeNodeStoppingConditionMFAIterations(nodeID));

  node->SetStopBiasCalculation
    (this->MRMLManager->GetTreeNodeBiasCalculationMaxIterations(nodeID));

  node->SetPrintShapeSimularityMeasure(0);         // !!!bcd!!!

  node->SetPCAShapeModelType(0);                   // !!!bcd!!!

  node->SetRegistrationIndependentSubClassFlag(0); // !!!bcd!!!
  node->SetRegistrationType(0);                    // !!!bcd!!!

  node->SetGenerateBackgroundProbability
    (this->MRMLManager->GetTreeNodeGenerateBackgroundProbability(nodeID));

  // New in 3.6. : Alpha now reflects user interface and is now correctly set for each parent node
  // cout << "Alpha setting for " << this->MRMLManager->GetTreeNodeName(nodeID) << " " << this->MRMLManager->GetTreeNodeAlpha(nodeID) << endl;
  node->SetAlpha(this->MRMLManager->GetTreeNodeAlpha(nodeID)); 
                      
}

//-----------------------------------------------------------------------------
void
vtkEMSegmentLogic::
CopyTreeLeafDataToSegmenter(vtkImageEMLocalClass* node, 
                            vtkIdType nodeID)
{
  unsigned int numTargetImages = 
    this->MRMLManager->GetTargetNumberOfSelectedVolumes();

  // this label describes the output intensity value for this class in
  // the segmentation result
  node->SetLabel(this->MRMLManager->GetTreeNodeIntensityLabel(nodeID));

  // set log mean and log covariance
  for (unsigned int r = 0; r < numTargetImages; ++r)
    {
    node->SetLogMu(this->MRMLManager->
                   GetTreeNodeDistributionLogMeanWithCorrection(nodeID, r), r);

    for (unsigned int c = 0; c < numTargetImages; ++c)
      {
      node->SetLogCovariance(this->MRMLManager->
                             GetTreeNodeDistributionLogCovarianceWithCorrection(nodeID,
                                                                  r, c), 
                             r, c);
      }
    }

  node->SetPrintQuality(this->MRMLManager->GetTreeNodePrintQuality(nodeID));
}

//-----------------------------------------------------------------------------
int
vtkEMSegmentLogic::
ConvertGUIEnumToAlgorithmEnumStoppingConditionType(int guiEnumValue)
{
  switch (guiEnumValue)
    {
    case (vtkEMSegmentMRMLManager::StoppingConditionIterations):
      return EMSEGMENT_STOP_FIXED;
    case (vtkEMSegmentMRMLManager::StoppingConditionLabelMapMeasure):
      return EMSEGMENT_STOP_LABELMAP;
    case (vtkEMSegmentMRMLManager::StoppingConditionWeightsMeasure):
      return EMSEGMENT_STOP_WEIGHTS;
    default:
      vtkErrorMacro("Unknown stopping condition type: " << guiEnumValue);
      return -1;
    }
}

//-----------------------------------------------------------------------------
int
vtkEMSegmentLogic::
ConvertGUIEnumToAlgorithmEnumInterpolationType(int guiEnumValue)
{
  switch (guiEnumValue)
    {
    case (vtkEMSegmentMRMLManager::InterpolationLinear):
      return EMSEGMENT_REGISTRATION_INTERPOLATION_LINEAR;
    case (vtkEMSegmentMRMLManager::InterpolationNearestNeighbor):
      return EMSEGMENT_REGISTRATION_INTERPOLATION_NEIGHBOUR;
    case (vtkEMSegmentMRMLManager::InterpolationCubic):
      // !!! not implemented
      vtkErrorMacro("Cubic interpolation not implemented: " << guiEnumValue);
      return -1;
    default:
      vtkErrorMacro("Unknown interpolation type: " << guiEnumValue);
      return -1;
    }
}

//----------------------------------------------------------------------------
vtkstd::string  vtkEMSegmentLogic::GetTclTaskDirectory(vtkSlicerApplication* app)
{
  //workaround for the mrml library, we need to have write access to this folder
  const char* tmp_dir = app->GetTemporaryDirectory();
  if (tmp_dir)
    {
      vtkstd::string copied_task_dir(std::string(tmp_dir) + std::string("/EMSegmentTaskCopy"));

      /**
        * Copy content directory to another directory with all files and
        * sub-directories.  If the "always" argument is true all files are
        * always copied.  If it is false, only files that have changed or
        * are new are copied.
        */
       // copy not always, only new files
       // Later do automatically
      vtkstd::string orig_task_dir = this->GetModuleShareDirectory() + vtkstd::string("/Tasks");
      
      if ( !vtksys::SystemTools::CopyADirectory(orig_task_dir.c_str(), copied_task_dir.c_str(), false, true) )
      {
          cout << "GetTclTaskDirectory:: Couldn't copy task directory " << orig_task_dir.c_str() << " to " << copied_task_dir.c_str() << endl;
          vtkErrorMacro("GetTclTaskDirectory:: Couldn't copy task directory " << orig_task_dir.c_str() << " to " << copied_task_dir.c_str());
          return vtksys::SystemTools::ConvertToOutputPath("");
      }
      return copied_task_dir;
    }
  else
    {
      // FIXME, make sure there is always a valid temporary directory
      vtkErrorMacro("GetTclTaskDirectory:: Tcl Task Directory was not found, set temporary directory first");
    }

  // return empty string if not found
  return vtksys::SystemTools::ConvertToOutputPath("");

}

//----------------------------------------------------------------------------
vtkstd::string  vtkEMSegmentLogic::GetTclGeneralDirectory()
{
  // Later do automatically
  vtkstd::string file_path = this->GetModuleShareDirectory() +  vtkstd::string("/Tcl");
  return vtksys::SystemTools::ConvertToOutputPath(file_path.c_str());
}

//----------------------------------------------------------------------------
std::string vtkEMSegmentLogic::DefineTclTaskFileFromMRML(vtkSlicerApplication *app)
{
  std::string tclFile("");
  tclFile = this->DefineTclTaskFullPathName(app, this->MRMLManager->GetTclTaskFilename());

  if (vtksys::SystemTools::FileExists(tclFile.c_str()) && (!vtksys::SystemTools::FileIsDirectory(tclFile.c_str())) )
    {
      return tclFile;
    }

  cout << "vtkEMSegmentLogic::DefineTclTaskFileFromMRML: " << tclFile.c_str() << " does not exist - using default file" << endl;

  tclFile = this->DefineTclTaskFullPathName(app, vtkMRMLEMSGlobalParametersNode::GetDefaultTaskTclFileName()); 
  return tclFile;  
}

void vtkEMSegmentLogic::TransferIJKToRAS(vtkMRMLVolumeNode* volumeNode, int ijk[3], double ras[3])
{
  vtkMatrix4x4* matrix = vtkMatrix4x4::New();
  volumeNode->GetIJKToRASMatrix(matrix);
  float input[4] = {ijk[0],ijk[1],ijk[2],1};
  float output[4];
  matrix->MultiplyPoint(input, output);
  ras[0]= output[0];
  ras[1]= output[1];
  ras[2]= output[2];
}

void vtkEMSegmentLogic::TransferRASToIJK(vtkMRMLVolumeNode* volumeNode, double ras[3], int ijk[3])
{
  vtkMatrix4x4* matrix = vtkMatrix4x4::New();
  volumeNode->GetRASToIJKMatrix(matrix);
  double input[4] = {ras[0],ras[1],ras[2],1};
  double output[4];
  matrix->MultiplyPoint(input, output);
  ijk[0]= int(output[0]);
  ijk[1]= int(output[1]);
  ijk[2]= int(output[2]);
}

// works for running stuff in TCL so that you do not need to look in two windows 
void vtkEMSegmentLogic::PrintText(char *TEXT) {
  cout << TEXT << endl;
} 

void vtkEMSegmentLogic::PrintTextNoNewLine(char *TEXT) {
  cout << TEXT;
  cout.flush();
} 

//-----------------------------------------------------------------------------
// Make sure you source EMSegmentAutoSample.tcl

int vtkEMSegmentLogic::ComputeIntensityDistributionsFromSpatialPrior(vtkKWApplication* app)
{
  // iterate over tree nodes
  typedef vtkstd::vector<vtkIdType>  NodeIDList;
  typedef NodeIDList::const_iterator NodeIDListIterator;
  NodeIDList nodeIDList;

  this->MRMLManager->GetListOfTreeNodeIDs(this->MRMLManager->GetTreeRootNodeID(), nodeIDList);
  for (NodeIDListIterator i = nodeIDList.begin(); i != nodeIDList.end(); ++i)
    {
      if (this->MRMLManager->GetTreeNodeIsLeaf(*i)) 
        {      
      this->UpdateIntensityDistributionAuto(app,*i);
        }
    }
  return 0;
}

//-----------------------------------------------------------------------------
void vtkEMSegmentLogic::UpdateIntensityDistributionAuto(vtkKWApplication* app, vtkIdType nodeID)
{

  if (!this->MRMLManager->GetTreeNodeSpatialPriorVolumeID(nodeID)) {
    vtkWarningMacro("Nothing to update for " << nodeID << " as atlas is not defined");
    return ;
  }

  vtkMRMLVolumeNode*  atlasNode = this->MRMLManager->GetAlignedSpatialPriorFromTreeNodeID(nodeID);
  if (!this->MRMLManager->GetTreeNodeSpatialPriorVolumeID(nodeID)) 
  {
    vtkErrorMacro("Atlas not yet aligned for " << nodeID << " ! ");
    return ;
  }

  // get working node 
  vtkMRMLEMSVolumeCollectionNode* workingTarget = NULL;
  if (this->MRMLManager->GetWorkingDataNode()->GetAlignedTargetNode() &&
      this->MRMLManager->GetWorkingDataNode()->GetAlignedTargetNodeIsValid())
    {
    workingTarget = this->MRMLManager->GetWorkingDataNode()->GetAlignedTargetNode();
    }
  else 
    {
       vtkErrorMacro("Cannot update intensity distribution bc Aligned Target is not correctly defined for node " << nodeID);
       return ;
    }

  int numTargetImages = workingTarget->GetNumberOfVolumes();
  
   // Sample
  {
    vtkstd::stringstream CMD ;
    CMD <<  "::EMSegmenterAutoSampleTcl::EMSegmentGaussCurveCalculationFromID " << vtkKWTkUtilities::GetTclNameFromPointer(app->GetMainInterp(), this) << " " << vtkKWTkUtilities::GetTclNameFromPointer(app->GetMainInterp(), this->MRMLManager) << " 0.95 1 { " ;
    for (int i = 0 ; i < numTargetImages; i++) {
      CMD << workingTarget->GetNthVolumeNodeID(i) << " " ;
    }
    CMD << " } ";
    CMD << atlasNode->GetID() << " {" <<  this->MRMLManager->GetTreeNodeName(nodeID) << "} \n";
    // cout << CMD.str().c_str() << endl;
    if (atoi(app->Script(CMD.str().c_str()))) { return; }
  }
  

  //
  // propagate data to mrml node
  //

  vtkMRMLEMSTreeParametersLeafNode* leafNode = this->MRMLManager->GetTreeParametersLeafNode(nodeID);  
  for (int r = 0; r < numTargetImages; ++r)
    {
      {
        double value = atof(app->Script("expr $::EMSegment(GaussCurveCalc,Mean,%d)",r));
        leafNode->SetLogMean(r, value);
      }
      for (int c = 0; c < numTargetImages; ++c)
      {
        double value = atof(app->Script("expr $::EMSegment(GaussCurveCalc,Covariance,%d,%d)",r,c));
        leafNode->SetLogCovariance(r, c, value);
      }
    }
}

//----------------------------------------------------------------------------
void  vtkEMSegmentLogic::AutoCorrectSpatialPriorWeight(vtkIdType nodeID)
{ 
   unsigned int numChildren = this->MRMLManager->GetTreeNodeNumberOfChildren(nodeID);
   for (unsigned int i = 0; i < numChildren; ++i)
    {
    vtkIdType childID = this->MRMLManager->GetTreeNodeChildNodeID(nodeID, i);
    bool isLeaf = this->MRMLManager->GetTreeNodeIsLeaf(childID);
    if (isLeaf)
      {
    if ((this->MRMLManager->GetTreeNodeSpatialPriorWeight(childID) > 0.0) && (!this->MRMLManager->GetAlignedSpatialPriorFromTreeNodeID(childID)))
      {
        vtkWarningMacro("Class with ID " <<  childID << " is set to 0 bc no atlas assigned to class!" );
        this->MRMLManager->SetTreeNodeSpatialPriorWeight(childID,0.0);
      }
      }
    else
      {
    this->AutoCorrectSpatialPriorWeight(childID);
      }
   }
}


//----------------------------------------------------------------------------
// cannot be moved to vtkEMSEgmentGUI bc of command line interface !
// This function is used for the UpdateButton in vtkEMSegmentParametersSetStep
vtkstd::string vtkEMSegmentLogic::GetTemporaryTaskDirectory(vtkSlicerApplication* app)
{
  // FIXME, what happens if user has no write permission to this directory
  std::string taskDir("");
  if (!app)
    {
      return taskDir;
    }

  const char* tmpDir = app->GetTemporaryDirectory();
  if (tmpDir)
    {
      std::string tmpTaskDir( std::string(tmpDir) + "/" + std::string(app->GetSvnRevision()) + std::string("/EMSegmentTask") );
      taskDir = vtksys::SystemTools::ConvertToOutputPath(tmpTaskDir.c_str());
    }
  else
    {
      // FIXME, make sure there is always a valid temporary directory
      vtkErrorMacro("GetTemporaryTaskDirectory:: Temporary Directory was not defined");
    }
  return taskDir;
} 

//----------------------------------------------------------------------------
// cannot be moved to vtkEMSEgmentGUI bc of command line interface !
std::string vtkEMSegmentLogic::DefineTclTaskFullPathName(vtkSlicerApplication* app, const char* TclFileName)
{

//  std::string task_dir = this->GetTclTaskDirectory(app);
//  cout << "TEST 1" << task_dir << " " << vtksys::SystemTools::FileExists(task_dir.c_str()) << endl;

  vtkstd::string tmp_full_file_path = this->GetTclTaskDirectory(app) + vtkstd::string("/") + vtkstd::string(TclFileName);
//  vtkstd::string full_file_path = vtksys::SystemTools::ConvertToOutputPath(tmp_full_file_path.c_str());
  if (vtksys::SystemTools::FileExists(tmp_full_file_path.c_str()))
    {
      return tmp_full_file_path;
    }

  tmp_full_file_path = this->GetTemporaryTaskDirectory(app) + vtkstd::string("/") + vtkstd::string(TclFileName);
//  full_file_path = vtksys::SystemTools::ConvertToOutputPath(tmp_full_file_path.c_str());
  if (vtksys::SystemTools::FileExists(tmp_full_file_path.c_str()))
    {
       return tmp_full_file_path;
    }

  vtkErrorMacro("DefineTclTaskFullPathName : could not find tcl file with name  " << TclFileName ); 
  tmp_full_file_path = vtkstd::string("");
  return  tmp_full_file_path;
}

bool vtkEMSegmentLogic::PackageAndWriteData(vtkSlicerApplication* app, vtkSlicerApplicationLogic* appLogic, const char* packageDirectory)
{
  //
  // create a scene and copy the EMSeg related nodes to it
  //
  if (!this->GetMRMLManager())
    {
      return false;
    }

  std::string outputDirectory(packageDirectory);
  std::string mrmlURL(outputDirectory + "/_EMSegmenterScene.mrml");

  vtkMRMLScene* newScene = vtkMRMLScene::New();
  newScene->SetRootDirectory(packageDirectory);
  newScene->SetURL(mrmlURL.c_str());

  vtkDataIOManagerLogic* dataIOManagerLogic = vtkDataIOManagerLogic::New();
  Slicer3Helper::AddDataIOToScene(newScene,app,appLogic,dataIOManagerLogic);

  // newScene->SetRootDirectory(outputDirectory.c_str());

  //std::cout << std::endl;
  this->GetMRMLManager()->CopyEMRelatedNodesToMRMLScene(newScene);

  // update filenames to match standardized package structure
  this->CreatePackageFilenames(newScene, packageDirectory);

  //
  // create directory structure on disk
  bool errorFlag = !this->CreatePackageDirectories(packageDirectory);

  if (errorFlag)
    {
    vtkErrorMacro("PackageAndWriteData: failed to create directories");
    }
  else 
    {
      //
      // write the scene out to disk
      errorFlag = !this->WritePackagedScene(newScene);
      if (errorFlag)
    {
      vtkErrorMacro("PackageAndWrite: failed to write scene");
    }
    }

    Slicer3Helper::RemoveDataIOFromScene(newScene,dataIOManagerLogic);
    dataIOManagerLogic->Delete();
    dataIOManagerLogic = NULL;
    newScene->Delete();

    return !errorFlag;
}


//-----------------------------------------------------------------------------
void
vtkEMSegmentLogic::
CreatePackageFilenames(vtkMRMLScene* scene, 
                       const char* vtkNotUsed(packageDirectoryName))
{
  //
  // set up mrml manager for this new scene
  vtkEMSegmentMRMLManager* newSceneManager = vtkEMSegmentMRMLManager::New();
  newSceneManager->SetMRMLScene(scene);
  vtkMRMLEMSTemplateNode* newEMSTemplateNode = dynamic_cast<vtkMRMLEMSTemplateNode*>(scene->GetNthNodeByClass(0, "vtkMRMLEMSTemplateNode"));
  if (newEMSTemplateNode == NULL)
    {
      vtkWarningMacro("CreatePackageFilenames: no EMSSegmenter node!");
      newSceneManager->Delete();
      return;
    }
  if (newSceneManager->SetNodeWithCheck(newEMSTemplateNode))
    {
       vtkWarningMacro("CreatePackageFilenames: not a valid template node!");
       newSceneManager->Delete();
       return;
    }
   
  vtkMRMLEMSWorkingDataNode* workingDataNode = 
    newSceneManager->GetWorkingDataNode();

  //
  // We might be creating volume storage nodes.  We must decide if the
  // images should be automatically centered when they are read.  Look
  // at the original input target node zero to decide if we will use
  // centering.
  bool centerImages = false;
  if (workingDataNode && workingDataNode->GetInputTargetNode())
    {
    if (workingDataNode->GetInputTargetNode()->GetNumberOfVolumes() > 0)
      {
    if (!workingDataNode->GetInputTargetNode()->GetNthVolumeNode(0)) 
      {
              vtkErrorMacro("CreatePackageFilenames: the first InputTagetNode is not defined!");
              vtkIndent ind;
          workingDataNode->GetInputTargetNode()->PrintSelf(cerr,ind);
              cout << endl;
      } 
        else 
          {
            vtkMRMLStorageNode* firstTargetStorageNode = workingDataNode->GetInputTargetNode()->GetNthVolumeNode(0)->GetStorageNode();
            vtkMRMLVolumeArchetypeStorageNode* firstTargetVolumeStorageNode = dynamic_cast<vtkMRMLVolumeArchetypeStorageNode*> (firstTargetStorageNode);
            if (firstTargetVolumeStorageNode != NULL)
            { 
             centerImages = firstTargetVolumeStorageNode->GetCenterImage();
            }
      }
       }
    }

   // get the full path to the scene
  std::vector<std::string> scenePathComponents;
  vtkstd::string rootDir = newSceneManager->GetMRMLScene()->GetRootDirectory();
  if (rootDir.find_last_of("/") == rootDir.length() - 1)
    {
      vtkDebugMacro("em seg: found trailing slash in : " << rootDir);
      rootDir = rootDir.substr(0, rootDir.length()-1);
    }
  vtkDebugMacro("em seg scene manager root dir = " << rootDir);
  vtksys::SystemTools::SplitPath(rootDir.c_str(), scenePathComponents);

  // change the storage file for the segmentation result
    {
    vtkMRMLVolumeNode* volumeNode = newSceneManager->GetOutputVolumeNode();
    if (volumeNode != NULL)
      {
      vtkMRMLStorageNode* storageNode = volumeNode->GetStorageNode();
      vtkMRMLVolumeArchetypeStorageNode* volumeStorageNode = 
        dynamic_cast<vtkMRMLVolumeArchetypeStorageNode*>(storageNode);
      if (volumeStorageNode == NULL)
      {
      // create a new storage node for this volume
      volumeStorageNode = vtkMRMLVolumeArchetypeStorageNode::New();
      scene->AddNodeNoNotify(volumeStorageNode);
      volumeNode->SetAndObserveStorageNodeID(volumeStorageNode->GetID());
      std::cout << "Added storage node : " << volumeStorageNode->GetID() 
                << std::endl;
      volumeStorageNode->Delete();
      storageNode = volumeStorageNode;
      }
      volumeStorageNode->SetCenterImage(centerImages);
    
      // create new filename
      std::string oldFilename       = 
        (storageNode->GetFileName() ? storageNode->GetFileName() :
         "SegmentationResult.mhd");
      std::string oldFilenameNoPath = 
        vtksys::SystemTools::GetFilenameName(oldFilename);

      scenePathComponents.push_back("Segmentation");
      scenePathComponents.push_back(oldFilenameNoPath);

      std::string newFilename = 
        vtksys::SystemTools::JoinPath(scenePathComponents);
      storageNode->SetFileName(newFilename.c_str());
      scenePathComponents.pop_back();
      scenePathComponents.pop_back();

      }
    }

  //
  // change the storage file for the targets
  int numTargets = newSceneManager->GetTargetNumberOfSelectedVolumes();

  // input target volumes
  if (workingDataNode->GetInputTargetNode())
    {
    for (int i = 0; i < numTargets; ++i)
      {
      vtkMRMLVolumeNode* volumeNode =
        workingDataNode->GetInputTargetNode()->GetNthVolumeNode(i);
      if (volumeNode != NULL)
        {
        vtkMRMLStorageNode* storageNode = volumeNode->GetStorageNode();
        vtkMRMLVolumeArchetypeStorageNode* volumeStorageNode = 
          dynamic_cast<vtkMRMLVolumeArchetypeStorageNode*>(storageNode);
        if (volumeStorageNode == NULL)
          {
          // create a new storage node for this volume
          volumeStorageNode = vtkMRMLVolumeArchetypeStorageNode::New();
          scene->AddNodeNoNotify(volumeStorageNode);
          volumeNode->SetAndObserveStorageNodeID(volumeStorageNode->GetID());
          std::cout << "Added storage node : " << volumeStorageNode->GetID() 
                    << std::endl;
          volumeStorageNode->Delete();
          storageNode = volumeStorageNode;
          }
        volumeStorageNode->SetCenterImage(centerImages);
      
        // create new filename
        vtkstd::stringstream defaultFilename;
        defaultFilename << "Target" << i << "_Input.mhd";
        std::string oldFilename       = 
          (storageNode->GetFileName() ? storageNode->GetFileName() :
           defaultFilename.str().c_str());
        std::string oldFilenameNoPath = 
          vtksys::SystemTools::GetFilenameName(oldFilename);
        scenePathComponents.push_back("Target");
        scenePathComponents.push_back("Input");
        scenePathComponents.push_back(oldFilenameNoPath);
        std::string newFilename = 
          vtksys::SystemTools::JoinPath(scenePathComponents);
        
        storageNode->SetFileName(newFilename.c_str());
        scenePathComponents.pop_back();
        scenePathComponents.pop_back();
        scenePathComponents.pop_back();
        }
      }  
    }

  // aligned target volumes
  if (workingDataNode->GetAlignedTargetNode())
    {
    for (int i = 0; i < numTargets; ++i)
      {
      vtkMRMLVolumeNode* volumeNode =
        workingDataNode->GetAlignedTargetNode()->GetNthVolumeNode(i);
      if (volumeNode != NULL)
        {
        vtkMRMLStorageNode* storageNode = volumeNode->GetStorageNode();
        vtkMRMLVolumeArchetypeStorageNode* volumeStorageNode = 
          dynamic_cast<vtkMRMLVolumeArchetypeStorageNode*>(storageNode);
        if (volumeStorageNode == NULL)
          {
          // create a new storage node for this volume
          volumeStorageNode = vtkMRMLVolumeArchetypeStorageNode::New();
          scene->AddNodeNoNotify(volumeStorageNode);
          volumeNode->SetAndObserveStorageNodeID(volumeStorageNode->GetID());
          std::cout << "Added storage node : " << volumeStorageNode->GetID() 
                    << std::endl;
          volumeStorageNode->Delete();
          storageNode = volumeStorageNode;
          }
        volumeStorageNode->SetCenterImage(centerImages);
          
        // create new filename
        vtkstd::stringstream defaultFilename;
        defaultFilename << "Target" << i << "_Aligned.mhd";
        std::string oldFilename       = 
          (storageNode->GetFileName() ? storageNode->GetFileName() :
           defaultFilename.str().c_str());
        std::string oldFilenameNoPath = 
          vtksys::SystemTools::GetFilenameName(oldFilename);
        scenePathComponents.push_back("Target");
        scenePathComponents.push_back("Aligned");
        scenePathComponents.push_back(oldFilenameNoPath);
        std::string newFilename = 
          vtksys::SystemTools::JoinPath(scenePathComponents);
        
        storageNode->SetFileName(newFilename.c_str());
    scenePathComponents.pop_back();
    scenePathComponents.pop_back();
    scenePathComponents.pop_back();
        }
      }  
    }

  //
  // change the storage file for the atlas
  int numAtlasVolumes = newSceneManager->GetAtlasInputNode()->
    GetNumberOfVolumes();

  // input atlas volumes
  if (newSceneManager->GetAtlasInputNode())
    {
    for (int i = 0; i < numAtlasVolumes; ++i)
      {
      vtkMRMLVolumeNode* volumeNode =
         newSceneManager->GetAtlasInputNode()->GetNthVolumeNode(i);
      if (volumeNode != NULL)
        {
        vtkMRMLStorageNode* storageNode = volumeNode->GetStorageNode();
        vtkMRMLVolumeArchetypeStorageNode* volumeStorageNode = 
          dynamic_cast<vtkMRMLVolumeArchetypeStorageNode*>(storageNode);
        if (volumeStorageNode == NULL)
          {
          // create a new storage node for this volume
          volumeStorageNode = vtkMRMLVolumeArchetypeStorageNode::New();
          scene->AddNodeNoNotify(volumeStorageNode);
          volumeNode->SetAndObserveStorageNodeID(volumeStorageNode->GetID());
          std::cout << "Added storage node : " << volumeStorageNode->GetID() 
                    << std::endl;
          volumeStorageNode->Delete();
          storageNode = volumeStorageNode;
          }
        volumeStorageNode->SetCenterImage(centerImages);
      
        // create new filename
        vtkstd::stringstream defaultFilename;
        defaultFilename << "Atlas" << i << "_Input.mhd";
        std::string oldFilename       = 
          (storageNode->GetFileName() ? storageNode->GetFileName() :
           defaultFilename.str().c_str());
        std::string oldFilenameNoPath = 
          vtksys::SystemTools::GetFilenameName(oldFilename);
        scenePathComponents.push_back("Atlas");
        scenePathComponents.push_back("Input");
        scenePathComponents.push_back(oldFilenameNoPath);
        std::string newFilename = 
          vtksys::SystemTools::JoinPath(scenePathComponents);
        
        storageNode->SetFileName(newFilename.c_str());
    scenePathComponents.pop_back();
    scenePathComponents.pop_back();
    scenePathComponents.pop_back();
        }
      }  
    }

  // aligned atlas volumes
  if (workingDataNode->GetAlignedAtlasNode())
    {
    for (int i = 0; i < numAtlasVolumes; ++i)
      {
      vtkMRMLVolumeNode* volumeNode =
        workingDataNode->GetAlignedAtlasNode()->GetNthVolumeNode(i);
      if (volumeNode != NULL)
        {
        vtkMRMLStorageNode* storageNode = volumeNode->GetStorageNode();
        vtkMRMLVolumeArchetypeStorageNode* volumeStorageNode = 
          dynamic_cast<vtkMRMLVolumeArchetypeStorageNode*>(storageNode);
        if (volumeStorageNode == NULL)
          {
          // create a new storage node for this volume
          volumeStorageNode = vtkMRMLVolumeArchetypeStorageNode::New();
          scene->AddNodeNoNotify(volumeStorageNode);
          volumeNode->SetAndObserveStorageNodeID(volumeStorageNode->GetID());
          std::cout << "Added storage node : " << volumeStorageNode->GetID() 
                    << std::endl;
          volumeStorageNode->Delete();
          storageNode = volumeStorageNode;
          }
        volumeStorageNode->SetCenterImage(centerImages);
        
        // create new filename
        vtkstd::stringstream defaultFilename;
        defaultFilename << "Atlas" << i << "_Aligned.mhd";
        std::string oldFilename       = 
          (storageNode->GetFileName() ? storageNode->GetFileName() :
           defaultFilename.str().c_str());
        std::string oldFilenameNoPath = 
          vtksys::SystemTools::GetFilenameName(oldFilename);
        scenePathComponents.push_back("Atlas");
        scenePathComponents.push_back("Aligned");
        scenePathComponents.push_back(oldFilenameNoPath);
        std::string newFilename = 
          vtksys::SystemTools::JoinPath(scenePathComponents);
        
        storageNode->SetFileName(newFilename.c_str());
    scenePathComponents.pop_back();
    scenePathComponents.pop_back();
    scenePathComponents.pop_back();
        }
      }  
    }

  // clean up
  newSceneManager->Delete();
}

//-----------------------------------------------------------------------------
bool
vtkEMSegmentLogic::
CreatePackageDirectories(const char* packageDirectoryName)
{
  vtkstd::string packageDirectory(packageDirectoryName);
  
  // check that parent directory exists
  std::string parentDirectory = 
    vtksys::SystemTools::GetParentDirectory(packageDirectory.c_str());
  if (!vtksys::SystemTools::FileExists(parentDirectory.c_str()))
    {
    vtkWarningMacro
      ("CreatePackageDirectories: Parent directory does not exist!");
    return false;
    }
  
  // create package directories
  bool createdOK = true;
  std::string newDir = packageDirectory + "/Atlas/Input";
  createdOK = createdOK &&
    vtksys::SystemTools::MakeDirectory(newDir.c_str());  
  newDir = packageDirectory + "/Atlas/Aligned";
  createdOK = createdOK &&
    vtksys::SystemTools::MakeDirectory(newDir.c_str());  
  newDir = packageDirectory + "/Target/Input";
  createdOK = createdOK &&
    vtksys::SystemTools::MakeDirectory(newDir.c_str());  
  newDir = packageDirectory + "/Target/Normalized";
  createdOK = createdOK &&
    vtksys::SystemTools::MakeDirectory(newDir.c_str());  
  newDir = packageDirectory + "/Target/Aligned";
  createdOK = createdOK &&
    vtksys::SystemTools::MakeDirectory(newDir.c_str());  
  newDir = packageDirectory + "/Segmentation";
  createdOK = createdOK &&
    vtksys::SystemTools::MakeDirectory(newDir.c_str());  

  if (!createdOK)
    {
    vtkWarningMacro("CreatePackageDirectories: Could not create directories!");
    return false;
    }

  return true;
}

//-----------------------------------------------------------------------------
bool
vtkEMSegmentLogic::
WritePackagedScene(vtkMRMLScene* scene)
{
  //
  // write the volumes
  scene->InitTraversal();
  vtkMRMLNode* currentNode;
  bool allOK = true;
  while ((currentNode = scene->GetNextNodeByClass("vtkMRMLVolumeNode")) &&
         (currentNode != NULL))
    {
    vtkMRMLVolumeNode* volumeNode = 
      dynamic_cast<vtkMRMLVolumeNode*>(currentNode);

    if (volumeNode == NULL)
      {
      vtkWarningMacro("Volume node is null for node: " 
                    << currentNode->GetID());
      scene->RemoveNode(currentNode);
      allOK = false;
      continue;
      }
    if (volumeNode->GetImageData() == NULL)
      {
    vtkWarningMacro("Volume data is null for volume node: " << currentNode->GetID() << " Name : " <<  (currentNode->GetName() ? currentNode->GetName(): "(none)" ));
      scene->RemoveNode(currentNode);
      allOK = false;
      continue;
      }
    if (volumeNode->GetStorageNode() == NULL)
      {
      vtkWarningMacro("Volume storage node is null for volume node: " 
                    << currentNode->GetID());
      scene->RemoveNode(currentNode);
      allOK = false;
      continue;
      }

    try
      {
      std::cout << "Writing volume: " << volumeNode->GetName() 
                << ": " << volumeNode->GetStorageNode()->GetFileName() << "...";
      volumeNode->GetStorageNode()->SetUseCompression(0);
      volumeNode->GetStorageNode()->WriteData(volumeNode);
      std::cout << "DONE" << std::endl;
      }
    catch (...)
      {
      vtkErrorMacro("Problem writing volume: " << volumeNode->GetID());
      allOK = false;
      }
    }
  
  //
  // write the MRML scene file
  try 
    {
    scene->Commit();
    }
  catch (...)
    {
    vtkErrorMacro("Problem writing scene.");
    allOK = false;
    }  

  return allOK;
}

//-----------------------------------------------------------------------------
void vtkEMSegmentLogic::SubParcelateSegmentation(vtkImageData* segmentation, vtkIdType nodeID)
{
  unsigned int numChildren =  this->MRMLManager->GetTreeNodeNumberOfChildren(nodeID);
  for (unsigned int i = 0; i < numChildren; ++i)
    {
    vtkIdType childID = this->MRMLManager->GetTreeNodeChildNodeID(nodeID, i);
    if (this->MRMLManager->GetTreeNodeIsLeaf(childID))
      {
    vtkMRMLVolumeNode*  parcellationNode =  this->MRMLManager->GetAlignedSubParcellationFromTreeNodeID(childID);
        if ( ! parcellationNode || !parcellationNode->GetImageData() ) 
      {
            continue;
      }
        int childLabel =  this->MRMLManager->GetTreeNodeIntensityLabel(childID);
        cout << "==> Subparcellate " <<  childLabel << endl;
        vtkImageData* input = vtkImageData::New();
    input->DeepCopy(segmentation);
 
    vtkImageThreshold* roiMap =     vtkImageThreshold::New();
        roiMap->SetInput(input);
        roiMap->ThresholdBetween( childLabel,  childLabel);
        roiMap->ReplaceOutOn();
        roiMap->SetInValue(1);
        roiMap->SetOutValue(0);
        roiMap->Update();
 
        vtkImageCast* castParcellation = vtkImageCast::New();
        castParcellation->SetInput(parcellationNode->GetImageData());
    castParcellation->SetOutputScalarType(roiMap->GetOutput()->GetScalarType());
    castParcellation->Update();

        vtkImageMathematics* roiParcellation = vtkImageMathematics::New();
    roiParcellation->SetInput1(roiMap->GetOutput());
    roiParcellation->SetInput2(castParcellation->GetOutput());
        roiParcellation->SetOperationToMultiply();
        roiParcellation->Update();

    vtkImageThreshold* changedSeg =     vtkImageThreshold::New();
        changedSeg->SetInput(input);
        changedSeg->ThresholdBetween( childLabel,  childLabel);
        changedSeg->ReplaceOutOff();
        changedSeg->SetInValue(0);
        changedSeg->Update();

        vtkImageMathematics* parcellatedSeg = vtkImageMathematics::New();
    parcellatedSeg->SetInput1(changedSeg->GetOutput());
    parcellatedSeg->SetInput2(roiParcellation->GetOutput());
        parcellatedSeg->SetOperationToAdd();
        parcellatedSeg->Update();

        segmentation->DeepCopy(parcellatedSeg->GetOutput());
        parcellatedSeg->Delete();
        changedSeg->Delete();
        roiParcellation->Delete();
        castParcellation->Delete();
        roiMap->Delete();
    input->Delete();
      }
    else
      {
        this->SubParcelateSegmentation(segmentation, childID); 
      }

    }
}
