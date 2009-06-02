/*==========================================================================

  Portions (c) Copyright 2008 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See Doc/copyright/copyright.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $HeadURL: http://svn.slicer.org/Slicer3/branches/EndoTracking/Modules/EndoNavIF/vtkIGTLToViewerTransform.cxx $
  Date:      $Date: 2009-05-01 16:40:42 -0400 (Fri, 01 May 2009) $
  Version:   $Revision: 9387 $

==========================================================================*/

#include <vtksys/SystemTools.hxx>

#include "vtkObjectFactory.h"
#include "vtkIGTLToViewerTransform.h"
#include "vtkProp3D.h"

#include "vtkMRMLLinearTransformNode.h"
#include "igtlTransformMessage.h"

vtkStandardNewMacro(vtkIGTLToViewerTransform);
vtkCxxRevisionMacro(vtkIGTLToViewerTransform, "$Revision: 9387 $");


//---------------------------------------------------------------------------
vtkIGTLToViewerTransform::vtkIGTLToViewerTransform()
{
  this->Viewer = NULL;
  this->NodeCreated = 0;
  this->ImageViewerCT = NULL;

}


//---------------------------------------------------------------------------
vtkIGTLToViewerTransform::~vtkIGTLToViewerTransform()
{
}


//---------------------------------------------------------------------------
void vtkIGTLToViewerTransform::PrintSelf(ostream& os, vtkIndent indent)
{
}


//---------------------------------------------------------------------------
vtkMRMLNode* vtkIGTLToViewerTransform::CreateNewNode(vtkMRMLScene* scene, const char* name)
{
  vtkMRMLLinearTransformNode* transformNode;

  transformNode = vtkMRMLLinearTransformNode::New();
  transformNode->SetName(name);
  transformNode->SetDescription("Received by EndoNav");

  vtkMatrix4x4* transform = vtkMatrix4x4::New();
  transform->Identity();
  //transformNode->SetAndObserveImageData(transform);
  transformNode->ApplyTransform(transform);
  transform->Delete();

  scene->AddNode(transformNode);  

  this->NodeCreated = 1;

  return transformNode;
}


//---------------------------------------------------------------------------
vtkIntArray* vtkIGTLToViewerTransform::GetNodeEvents()
{
  return NULL;
}


//---------------------------------------------------------------------------
int vtkIGTLToViewerTransform::IGTLToMRML(igtl::MessageBase::Pointer buffer, vtkMRMLNode* node)
{
  // Create a message buffer to receive transform data
  igtl::TransformMessage::Pointer transMsg;
  transMsg = igtl::TransformMessage::New();
  transMsg->Copy(buffer);  // !! TODO: copy makes performance issue.

  // Deserialize the transform data
  // If you want to skip CRC check, call Unpack() without argument.
  int c = transMsg->Unpack(1);
  if (!(c & igtl::MessageHeader::UNPACK_BODY)) // if CRC check fails
    {
    // TODO: error handling
    return 0;
    }

  if (node == NULL)
    {
    return 0;
    }

  vtkMRMLLinearTransformNode* transformNode = 
    vtkMRMLLinearTransformNode::SafeDownCast(node);

  igtl::Matrix4x4 matrix;
  transMsg->GetMatrix(matrix);

  float tx = matrix[0][0];
  float ty = matrix[1][0];
  float tz = matrix[2][0];
  float sx = matrix[0][1];
  float sy = matrix[1][1];
  float sz = matrix[2][1];
  float nx = matrix[0][2];
  float ny = matrix[1][2];
  float nz = matrix[2][2];
  float px = matrix[0][3];
  float py = matrix[1][3];
  float pz = matrix[2][3];

  //std::cerr << "\n\nmatrix = " << std::endl;
  //std::cerr << tx << ", " << ty << ", " << tz << std::endl;
  //std::cerr << sx << ", " << sy << ", " << sz << std::endl;
  //std::cerr << nx << ", " << ny << ", " << nz << std::endl;
  //std::cerr << px << ", " << py << ", " << pz << std::endl;
  
  // set volume orientation
  vtkMatrix4x4* transform = vtkMatrix4x4::New();
  vtkMatrix4x4* transformToParent = transformNode->GetMatrixTransformToParent();

  transform->Identity();
  transform->SetElement(0, 0, tx);
  transform->SetElement(1, 0, ty);
  transform->SetElement(2, 0, tz);

  transform->SetElement(0, 1, sx);
  transform->SetElement(1, 1, sy);
  transform->SetElement(2, 1, sz);

  transform->SetElement(0, 2, nx);
  transform->SetElement(1, 2, ny);
  transform->SetElement(2, 2, nz);

  transform->SetElement(0, 3, px);
  transform->SetElement(1, 3, py);
  transform->SetElement(2, 3, pz);

  if (this->NodeCreated)
    {
    transformToParent->DeepCopy(transform);
    this->LocatorID = this->GetLocatorActorID(transformNode->GetScene());
    }
  else
    {
    if (this->LocatorID == "")
      {
      this->LocatorID = this->GetLocatorActorID(transformNode->GetScene());
      }
    if (std::string(transformNode->GetName()) == std::string("EnodNavSensor"))
      {
      transformNode->DisableModifiedEventOn();
      transformNode->SetAndObserveMatrixTransformToParent(NULL);

      if ( this->Viewer)
        {
        vtkProp3D *prop = this->Viewer->GetActorByID(this->LocatorID.c_str());
        if (prop)
          {
          prop->SetUserMatrix(transform);
          this->Viewer->GetMainViewer()->Render();
          }
        }

      transformNode->SetAndObserveMatrixTransformToParent(transform);

      transformNode->DisableModifiedEventOff();
      }

    }

  //std::cerr << "IGTL matrix = " << std::endl;
  //transform->Print(cerr);
  //std::cerr << "MRML matrix = " << std::endl;
  //transformToParent->Print(cerr);

  transform->Delete();

  this->NodeCreated = 0;

  return 1;
}


//---------------------------------------------------------------------------
int vtkIGTLToViewerTransform::MRMLToIGTL(unsigned long event, vtkMRMLNode* mrmlNode, int* size, void** igtlMsg)
{

  return 0;

}

//---------------------------------------------------------------------------
std::string vtkIGTLToViewerTransform::GetLocatorActorID(vtkMRMLScene*  scene)
{
  std::string locatorID = "";
  vtkMRMLModelNode*   locatorModel = NULL;
  vtkMRMLDisplayNode* locatorDisp = NULL;

  vtkCollection* collection = scene->GetNodesByName("IGTLLocator");

  if (collection != NULL && collection->GetNumberOfItems() == 1)
    {
    locatorModel = vtkMRMLModelNode::SafeDownCast(collection->GetItemAsObject(0));
    }

  if (locatorModel)
    {
    locatorDisp = locatorModel->GetDisplayNode();
    if (locatorDisp)
      {
      locatorID = locatorDisp->GetID();
      }
    }
  return locatorID;

}
