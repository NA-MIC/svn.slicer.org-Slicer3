#include <string>
#include <iostream>
#include <sstream>

#include "vtkObject.h"
#include "vtkObjectFactory.h"
#include "vtkCallbackCommand.h"

#include "vtkAbstractTransform.h"
#include "vtkMath.h"

#include "vtkMRMLROINode.h"
#include "vtkMRMLScene.h"

//------------------------------------------------------------------------------
vtkMRMLROINode* vtkMRMLROINode::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkMRMLROINode");
  if(ret)
    {
    return (vtkMRMLROINode*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkMRMLROINode;
}
 
//-----------------------------------------------------------------------------
vtkMRMLNode* vtkMRMLROINode::CreateNodeInstance()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkMRMLROINode");
  if(ret)
    {
    return (vtkMRMLROINode*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkMRMLROINode;
}

//----------------------------------------------------------------------------
void vtkMRMLROINode::PrintSelf(ostream& os, vtkIndent indent)
{  
  vtkObject::PrintSelf(os,indent);
  //Fill in ......
  return;
}

//----------------------------------------------------------------------------
vtkMRMLROINode::vtkMRMLROINode()
{
  this->XYZ[0] = this->XYZ[1] = this->XYZ[2] = 0.0;
  this->RadiusXYZ[0] = this->RadiusXYZ[1] = this->RadiusXYZ[2] = 10.0;

  this->IJK[0] = this->IJK[1] = this->IJK[2] = 0;
  this->RadiusIJK[0] = this->RadiusIJK[1] = this->RadiusIJK[2] = 0;
  // so that the SetLabelText macro won't try to free memory
  this->LabelText = NULL;
  this->SetLabelText(""); 
  this->ID = NULL;
  this->SetID("");
  this->Selected = 0;
  this->VolumeNodeID = NULL;
  this->Visibility = 1;
  this->InteractiveMode = 1;
  this->HideFromEditors = 0;
  return;
}

//----------------------------------------------------------------------------
vtkMRMLROINode::~vtkMRMLROINode()
{
  if (this->LabelText)
    {
    delete [] this->LabelText;
    this->LabelText = NULL;
    }
  if (this->ID)
    {
    delete [] this->ID;
    this->ID = NULL;
    }
  if (this->VolumeNodeID)
    {
    delete [] this->VolumeNodeID;
    this->VolumeNodeID = NULL;
    }
  return;
}

//----------------------------------------------------------------------------
void vtkMRMLROINode::WriteXML(ostream& of, int nIndent)
{
  // Write all attributes not equal to their defaults
  //Superclass::WriteXML(of, nIndent);

  //vtkIndent indent(nIndent);

  if (this->ID != NULL) 
    {
    of <<  " ROINodeID " << this->ID;
    }
  if (this->VolumeNodeID != NULL)
    {
    of << " VolumeNodeID " << this->VolumeNodeID;
    }
  if (this->LabelText != NULL)
    {
    of << " Labeltext " << this->LabelText;
    }

  of <<  " XYZ " 
    << this->XYZ[0] << " " << this->XYZ[1] << " " << this->XYZ[2];

  of <<  " RadiusXYZ " 
    << this->RadiusXYZ[0] << " " << this->RadiusXYZ[1] << " " << this->RadiusXYZ[2];
  of << " Selected " << this->Selected;

  of << " visibility=\"" << (this->Visibility ? "true" : "false") << "\"";

  of << " interactiveMode=\"" << (this->InteractiveMode ? "true" : "false") << "\"";


  return;
}

//----------------------------------------------------------------------------
void vtkMRMLROINode::ReadXMLAttributes( const char** atts)
{
  const char* attName;
  const char* attValue;

  while (*atts != NULL) 
    {
    attName = *(atts++);
    attValue = *(atts++);

    if (!strcmp(attName, "XYZ")) 
      {
      std::stringstream ss;
      double val;
      ss << attValue;
      for(int i=0; i<3; i++) 
        {
        ss >> val;
        this->XYZ[i] = val;
        }
      }
    if (!strcmp(attName, "RadiusXYZ")) 
      {
      std::stringstream ss;
      double val;
      ss << attValue;
      for(int i=0; i<3; i++) 
        {
        ss >> val;
        this->RadiusXYZ[i] = val;
        }
      }
    if (!strcmp(attName, "Selected")) 
      {
      std::stringstream ss;
      ss << attValue;
      ss >> this->Selected;
      }
    else if (!strcmp(attName, "ROINodeID")) 
      {
      this->SetID(attValue);
      }
    else if (!strcmp(attName, "VolumeNodeID")) 
      {
      this->SetVolumeNodeID(attValue);
      }
    else if (!strcmp(attName, "LabelText"))
      {
      this->SetLabelText(attValue);
      }
    else if (!strcmp(attName, "visibility")) 
      {
      if (!strcmp(attValue,"true")) 
        {
        this->Visibility = 1;
        }
      else
        {
        this->Visibility = 0;
        }
      }
    else if (!strcmp(attName, "interactiveMode")) 
      {
      if (!strcmp(attValue,"true")) 
        {
        this->InteractiveMode = 1;
        }
      else
        {
        this->InteractiveMode = 0;
        }
      }
 
    }

  return;
}

void vtkMRMLROINode::ReadXMLString(const char *keyValuePairs)
{
  // used because the ROI list gloms together the point's key and
  // values into one long string, VERY dependent on the order it's written
  // out in when WriteXML is used

  // insert the string into a stream
  std::stringstream ss;
  ss << keyValuePairs;

  char keyName[1024];

  // get out the id
  ss >> keyName;
  ss >> this->ID;
  vtkDebugMacro("ReadXMLString: got id " << this->ID);

  // get out the volume id
  ss >> keyName;
  if (!strcmp(keyName, "VolumeNodeID"))
    {
    char* IDValue = new char[1024];
    ss >> IDValue;
    this->SetVolumeNodeID(IDValue);
    delete [] IDValue;
    vtkDebugMacro("ReadXMLString: got VolumeNodeID " << this->VolumeNodeID);
    }
  else
    {
    // now get the label text value
    ss >> this->LabelText;
    vtkDebugMacro("ReadXMLString: got label text " << this->LabelText);
    }
  
  // get the xyz key
  ss >> keyName;
  // now get the x, y, z values
  ss >> this->XYZ[0];
  ss >> this->XYZ[1];
  ss >> this->XYZ[2];

  // get the Radiusxyz key
  ss >> keyName;
  // now get the x, y, z values
  ss >> this->RadiusXYZ[0];
  ss >> this->RadiusXYZ[1];
  ss >> this->RadiusXYZ[2];

  // get the selected flag
  ss >> keyName;
  ss >> this->Selected;
  return;
}

//----------------------------------------------------------------------------
// Copy the node's attributes to this object.
// Does NOT copy: ID, FilePrefix, Name, ID
void vtkMRMLROINode::Copy(vtkMRMLNode *anode)
{
  Superclass::Copy(anode);
  this->DisableModifiedEventOn();

  //  vtkObject::Copy(anode);
  vtkMRMLROINode *node = (vtkMRMLROINode *) anode;
  this->SetXYZ(node->XYZ);
  this->SetRadiusXYZ(node->RadiusXYZ);
  this->SetLabelText(node->GetLabelText());
  this->SetID(node->ID);
  this->SetSelected(node->GetSelected());
  this->SetVisibility(node->Visibility);
  this->SetInteractiveMode(node->InteractiveMode);

  this->DisableModifiedEventOff();
  this->InvokePendingModifiedEvent();

  return;
}

//-----------------------------------------------------------------------------
void vtkMRMLROINode::UpdateScene(vtkMRMLScene *scene)
{
  Superclass::UpdateScene(scene);
  return;
}

//-----------------------------------------------------------------------------
void vtkMRMLROINode::UpdateReferences()
{
  Superclass::UpdateReferences();
  return;
}

//-----------------------------------------------------------------------------
void vtkMRMLROINode::ProcessMRMLEvents ( vtkObject *caller,
                                        unsigned long event, 
                                        void *callData )
{
  Superclass::ProcessMRMLEvents(caller, event, callData);
  return;
}

//-----------------------------------------------------------------------------
void vtkMRMLROINode::SetXYZ(float X, float Y, float Z)
{
  this->XYZ[0] = X;
  this->XYZ[1] = Y;
  this->XYZ[2] = Z;

  this->Modified();
  return;
}

//-----------------------------------------------------------------------------
void vtkMRMLROINode::SetXYZ(float* xyz)
{
  this->SetXYZ(xyz[0], xyz[1], xyz[2]);
  return;
}

//-----------------------------------------------------------------------------
void vtkMRMLROINode::SetRadiusXYZ(float RadiusX, float RadiusY, float RadiusZ)
{ 
  this->RadiusXYZ[0] = RadiusX;
  this->RadiusXYZ[1] = RadiusY;
  this->RadiusXYZ[2] = RadiusZ;

  this->Modified();
  return;
}

//-----------------------------------------------------------------------------
void vtkMRMLROINode::SetRadiusXYZ(float* radiusXYZ)
{
  this->SetRadiusXYZ(radiusXYZ[0], radiusXYZ[1], radiusXYZ[2]);
  return;
}

//-----------------------------------------------------------------------------
void vtkMRMLROINode::SetIJK(float I, float J, float K)
{
  this->IJK[0] = I;
  this->IJK[1] = J;
  this->IJK[2] = K;

  //Update  

  return;
}

//-----------------------------------------------------------------------------
void vtkMRMLROINode::SetIJK(float* ijk)
{
  this->SetIJK(ijk[0], ijk[1], ijk[2]);
  this->Modified();
  return;
}

//-----------------------------------------------------------------------------
void vtkMRMLROINode::SetRadiusIJK(float RadiusI, float RadiusJ, float RadiusK)
{
  this->RadiusIJK[0] = RadiusI;
  this->RadiusIJK[1] = RadiusJ;
  this->RadiusIJK[2] = RadiusK;

  //Update

  this->Modified();
  return;
}

//-----------------------------------------------------------------------------
void vtkMRMLROINode::SetRadiusIJK(float* radiusIJK)
{
  this->SetRadiusIJK(radiusIJK[0], radiusIJK[1], radiusIJK[2]);
  return;
}

//---------------------------------------------------------------------------
void vtkMRMLROINode::ApplyTransform(vtkMatrix4x4* transformMatrix)
{
  double (*matrix)[4] = transformMatrix->Element;
  float *xyzIn  = this->GetXYZ();
  float xyzOut[3];

  xyzOut[0] = matrix[0][0]*xyzIn[0] + matrix[0][1]*xyzIn[1] + matrix[0][2]*xyzIn[2] + matrix[0][3];
  xyzOut[1] = matrix[1][0]*xyzIn[0] + matrix[1][1]*xyzIn[1] + matrix[1][2]*xyzIn[2] + matrix[1][3];
  xyzOut[2] = matrix[2][0]*xyzIn[0] + matrix[2][1]*xyzIn[1] + matrix[2][2]*xyzIn[2] + matrix[2][3];

  this->SetXYZ(xyzOut);

}

void vtkMRMLROINode::ApplyTransform(vtkAbstractTransform* transform)
{
  float *xyzIn  = this->GetXYZ();
  float xyzOut[3];
  transform->TransformPoint(xyzIn,xyzOut);
  this->SetXYZ(xyzOut);
}

