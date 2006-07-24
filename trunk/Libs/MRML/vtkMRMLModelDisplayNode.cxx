/*=auto=========================================================================

Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

See Doc/copyright/copyright.txt
or http://www.slicer.org/copyright/copyright.txt for details.

Program:   3D Slicer
Module:    $RCSfile: vtkMRMLModelDisplayNode.cxx,v $
Date:      $Date: 2006/03/03 22:26:39 $
Version:   $Revision: 1.3 $

=========================================================================auto=*/
#include <string>
#include <iostream>
#include <sstream>

#include "vtkObjectFactory.h"
#include "vtkCallbackCommand.h"

#include "vtkMRMLModelDisplayNode.h"
#include "vtkMRMLScene.h"

//------------------------------------------------------------------------------
vtkMRMLModelDisplayNode* vtkMRMLModelDisplayNode::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkMRMLModelDisplayNode");
  if(ret)
    {
    return (vtkMRMLModelDisplayNode*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkMRMLModelDisplayNode;
}

//-----------------------------------------------------------------------------

vtkMRMLNode* vtkMRMLModelDisplayNode::CreateNodeInstance()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkMRMLModelDisplayNode");
  if(ret)
    {
    return (vtkMRMLModelDisplayNode*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkMRMLModelDisplayNode;
}


//----------------------------------------------------------------------------
vtkMRMLModelDisplayNode::vtkMRMLModelDisplayNode()
{

  // Strings
  this->Color[0] = 0.5;
  this->Color[1] = 0.5;
  this->Color[2] = 0.5;

  // Numbers
  this->Opacity = 1.0;
  this->Visibility = 1;
  this->Clipping = 0;
  this->BackfaceCulling = 1;
  this->ScalarVisibility = 0;
  this->VectorVisibility = 0;
  this->TensorVisibility = 0;
  
  this->Ambient = 0;
  this->Diffuse = 100;
  this->Specular = 0;
  this->Power = 1;
  
  // Arrays
  this->ScalarRange[0] = 0;
  this->ScalarRange[1] = 100;

  // Scalars
  this->LUTName = -1;
  
  this->TextureImageData = NULL;

}

//----------------------------------------------------------------------------
vtkMRMLModelDisplayNode::~vtkMRMLModelDisplayNode()
{
  if (this->TextureImageData) 
    {
    this->TextureImageData->Delete();
    }
}

//----------------------------------------------------------------------------
void vtkMRMLModelDisplayNode::WriteXML(ostream& of, int nIndent)
{
  // Write all attributes not equal to their defaults
  
  Superclass::WriteXML(of, nIndent);

  vtkIndent indent(nIndent);

  if (this->Color)
    {
    of << indent << " color=\"" << this->Color[0] << " "
      << this->Color[1] << " "
      << this->Color[2] << "\"";
    }

  //if (this->LUTName && strcmp(this->LUTName,""))
  if (this->LUTName != -1)
    {
    of << indent << " lutName=\"" << this->LUTName << "\"";
    }
  
  of << indent << " ambient=\"" << this->Ambient << "\"";

  of << indent << " diffuse=\"" << this->Diffuse << "\"";

  of << indent << " specular=\"" << this->Specular << "\"";

  of << indent << " power=\"" << this->Power << "\"";

  of << indent << " opacity=\"" << this->Opacity << "\"";

  of << indent << " visibility=\"" << (this->Visibility ? "true" : "false") << "\"";

  of << indent << " clipping=\"" << (this->Clipping ? "true" : "false") << "\"";

  of << indent << " backfaceCulling=\"" << (this->BackfaceCulling ? "true" : "false") << "\"";

  of << indent << " scalarVisibility=\"" << (this->ScalarVisibility ? "true" : "false") << "\"";

  of << indent << " scalarRange=\"" << this->ScalarRange[0] << " "
     << this->ScalarRange[1] << "\"";
}

//----------------------------------------------------------------------------
void vtkMRMLModelDisplayNode::ReadXMLAttributes(const char** atts)
{

  vtkMRMLNode::ReadXMLAttributes(atts);

  const char* attName;
  const char* attValue;
  while (*atts != NULL) 
    {
    attName = *(atts++);
    attValue = *(atts++);
    if (!strcmp(attName, "color")) 
      {
      std::stringstream ss;
      ss << attValue;
      ss >> Color[0];
      ss >> Color[1];
      ss >> Color[2];
      }
    else if (!strcmp(attName, "scalarRange")) 
      {
      std::stringstream ss;
      ss << attValue;
      ss >> ScalarRange[0];
      ss >> ScalarRange[1];
      }
    else if (!strcmp(attName, "LUTName")) 
      {
      std::stringstream ss;
      ss << attValue;
      ss >> LUTName;
      }
    else if (!strcmp(attName, "ambient")) 
      {
      std::stringstream ss;
      ss << attValue;
      ss >> Ambient;
      }
    else if (!strcmp(attName, "diffuse")) 
      {
      std::stringstream ss;
      ss << attValue;
      ss >> Diffuse;
      }
    else if (!strcmp(attName, "specular")) 
      {
      std::stringstream ss;
      ss << attValue;
      ss >> Specular;
      }
    else if (!strcmp(attName, "power")) 
      {
      std::stringstream ss;
      ss << attValue;
      ss >> Power;
      }
    else if (!strcmp(attName, "opacity")) 
      {
      std::stringstream ss;
      ss << attValue;
      ss >> Opacity;
      }
    else if (!strcmp(attName, "visibility")) 
      {
      std::stringstream ss;
      ss << attValue;
      ss >> Visibility;
      }
    else if (!strcmp(attName, "backfaceCulling")) 
      {
      std::stringstream ss;
      ss << attValue;
      ss >> BackfaceCulling;
      }
    else if (!strcmp(attName, "scalarVisibility")) 
      {
      std::stringstream ss;
      ss << attValue;
      ss >> ScalarVisibility;
      }
    else if (!strcmp(attName, "vectorVisibility")) 
      {
      std::stringstream ss;
      ss << attValue;
      ss >> VectorVisibility;
      }
    else if (!strcmp(attName, "tensorVisibility")) 
      {
      std::stringstream ss;
      ss << attValue;
      ss >> TensorVisibility;
      }
    }  
}


//----------------------------------------------------------------------------
// Copy the node's attributes to this object.
// Does NOT copy: ID, FilePrefix, Name, ID
void vtkMRMLModelDisplayNode::Copy(vtkMRMLNode *anode)
{
  vtkMRMLNode::Copy(anode);
  vtkMRMLModelDisplayNode *node = (vtkMRMLModelDisplayNode *) anode;

  // Strings

  this->SetColor(node->Color);

  // Vectors
  this->SetScalarRange(node->ScalarRange);
  
  // Numbers
  this->SetOpacity(node->Opacity);
  this->SetAmbient(node->Ambient);
  this->SetDiffuse(node->Diffuse);
  this->SetSpecular(node->Specular);
  this->SetPower(node->Power);
  this->SetVisibility(node->Visibility);
  this->SetScalarVisibility(node->ScalarVisibility);
  this->SetBackfaceCulling(node->BackfaceCulling);
  this->SetClipping(node->Clipping);
  this->SetAndObserveTextureImageData(node->TextureImageData);

}

//----------------------------------------------------------------------------
void vtkMRMLModelDisplayNode::PrintSelf(ostream& os, vtkIndent indent)
{
  int idx;
  
  vtkMRMLNode::PrintSelf(os,indent);

  os << indent << "Color:             " << this->Color << "\n";
  os << indent << "Opacity:           " << this->Opacity << "\n";
  os << indent << "Ambient:           " << this->Ambient << "\n";
  os << indent << "Diffuse:           " << this->Diffuse << "\n";
  os << indent << "Specular:          " << this->Specular << "\n";
  os << indent << "Power:             " << this->Power << "\n";
  os << indent << "Visibility:        " << this->Visibility << "\n";
  os << indent << "ScalarVisibility:  " << this->ScalarVisibility << "\n";
  os << indent << "BackfaceCulling:   " << this->BackfaceCulling << "\n";
  os << indent << "Clipping:          " << this->Clipping << "\n";

  os << "ScalarRange:\n";
  for (idx = 0; idx < 2; ++idx)
    {
    os << indent << ", " << this->ScalarRange[idx];
    }

}

//----------------------------------------------------------------------------
void vtkMRMLModelDisplayNode::SetAndObserveTextureImageData(vtkImageData *ImageData)
{
  if (this->TextureImageData != NULL)
    {
    this->TextureImageData->RemoveObservers ( vtkCommand::ModifiedEvent, this->MRMLCallbackCommand );
    }

  this->SetTextureImageData(ImageData);
  if (this->TextureImageData != NULL)
    {
    this->TextureImageData->AddObserver ( vtkCommand::ModifiedEvent, this->MRMLCallbackCommand );
    }
}

//---------------------------------------------------------------------------
void vtkMRMLModelDisplayNode::ProcessMRMLEvents ( vtkObject *caller,
                                           unsigned long event, 
                                           void *callData )
{
  Superclass::ProcessMRMLEvents(caller, event, callData);

  if (this->TextureImageData == vtkImageData::SafeDownCast(caller) &&
    event ==  vtkCommand::ModifiedEvent)
    {
    this->InvokeEvent(vtkCommand::ModifiedEvent, NULL);
    }
  return;
}
