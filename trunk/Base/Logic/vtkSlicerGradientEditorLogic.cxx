#include "vtkObject.h"
#include "vtkObjectFactory.h"

#include "vtkSlicerGradientEditorLogic.h"

#include "vtkMRMLDiffusionWeightedVolumeNode.h"
#include "vtkMRMLNRRDStorageNode.h"
#include "vtkDoubleArray.h"
#include "vtkNRRDReader.h"

vtkCxxRevisionMacro(vtkSlicerGradientEditorLogic, "$Revision: 1.9.12.1 $");
vtkStandardNewMacro(vtkSlicerGradientEditorLogic);

//---------------------------------------------------------------------------
vtkSlicerGradientEditorLogic::vtkSlicerGradientEditorLogic(void)
  {
  }

//---------------------------------------------------------------------------
vtkSlicerGradientEditorLogic::~vtkSlicerGradientEditorLogic(void)
  {
  }

//---------------------------------------------------------------------------
void vtkSlicerGradientEditorLogic::PrintSelf ( ostream& os, vtkIndent indent )
  {
  this->vtkObject::PrintSelf ( os, indent );
  os << indent << "vtkSlicerGradientEditorLogic: " << this->GetClassName ( ) << "\n";
  }

//---------------------------------------------------------------------------
vtkMRMLDiffusionWeightedVolumeNode* vtkSlicerGradientEditorLogic::AddGradients (const char* filename)
  {
  // format the filename
  std::string fileString(filename);
  for (unsigned int i = 0; i < fileString.length(); i++)
    {
    if (fileString[i] == '\\')
      {
      fileString[i] = '/';
      }
    }

  vtkMRMLDiffusionWeightedVolumeNode *dwiNode = vtkMRMLDiffusionWeightedVolumeNode::New();

  // Instanciation of the I/O mechanism
  vtkMRMLNRRDStorageNode *storageNode = vtkMRMLNRRDStorageNode::New();
  storageNode->SetFileName(fileString.c_str());

  if (!storageNode->ReadData(dwiNode))
    {
    //TODO: txt File?
    vtkNRRDReader* reader = vtkNRRDReader::New();
    reader->CanReadFile(fileString.c_str());
    reader->GetClassName();
    }

  dwiNode->Delete();
  return dwiNode;
  }

