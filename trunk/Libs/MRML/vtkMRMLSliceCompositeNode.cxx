/*=auto=========================================================================

Portions (c) Copyright 2005 Brigham and Women\"s Hospital (BWH) All Rights Reserved.

See Doc/copyright/copyright.txt
or http://www.slicer.org/copyright/copyright.txt for details.

Program:   3D Slicer
Module:    $RCSfile: vtkMRMLSliceCompositeNode.cxx,v $
Date:      $Date: 2006/03/17 15:10:10 $
Version:   $Revision: 1.2 $

=========================================================================auto=*/

#include <string>
#include <iostream>
#include <sstream>

#include "vtkObjectFactory.h"
#include "vtkMRMLSliceCompositeNode.h"
#include "vtkMRMLScene.h"

#include "vtkMatrix4x4.h"

//------------------------------------------------------------------------------
vtkMRMLSliceCompositeNode* vtkMRMLSliceCompositeNode::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkMRMLSliceCompositeNode");
  if(ret)
    {
    return (vtkMRMLSliceCompositeNode*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkMRMLSliceCompositeNode;
}

//----------------------------------------------------------------------------

vtkMRMLNode* vtkMRMLSliceCompositeNode::CreateNodeInstance()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkMRMLSliceCompositeNode");
  if(ret)
    {
    return (vtkMRMLSliceCompositeNode*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkMRMLSliceCompositeNode;
}

//----------------------------------------------------------------------------
vtkMRMLSliceCompositeNode::vtkMRMLSliceCompositeNode()
{
  this->BackgroundVolumeID = NULL;
  this->ForegroundVolumeID = NULL;
  this->LabelVolumeID = NULL;
  this->ForegroundOpacity = 0.0; // start by showing only the background volume
  this->LabelOpacity = 1.0; // Show the label if there is one
  this->LayoutName = NULL;
}

//----------------------------------------------------------------------------
vtkMRMLSliceCompositeNode::~vtkMRMLSliceCompositeNode()
{
}

//----------------------------------------------------------------------------
void vtkMRMLSliceCompositeNode::WriteXML(ostream& of, int nIndent)
{
  Superclass::WriteXML(of, nIndent);

  vtkIndent indent(nIndent);

  of << indent << "backgroundVolumeID=\"" << 
   (this->BackgroundVolumeID ? this->BackgroundVolumeID : "") << "\" ";
  of << indent << "foregroundVolumeID=\"" << 
   (this->ForegroundVolumeID ? this->ForegroundVolumeID : "") << "\" ";
  of << indent << "labelVolumeID=\"" << 
   (this->LabelVolumeID ? this->LabelVolumeID : "") << "\" ";

  of << indent << "labelOpacity=\"" << this->LabelOpacity << "\" ";
  of << indent << "layoutName=\"" << this->LayoutName << "\" ";
}

//----------------------------------------------------------------------------
void vtkMRMLSliceCompositeNode::ReadXMLAttributes(const char** atts)
{

  Superclass::ReadXMLAttributes(atts);

  const char* attName;
  const char* attValue;
  while (*atts != NULL) 
    {
    attName = *(atts++);
    attValue = *(atts++);
    if (!strcmp(attName, "backgroundVolumeID")) 
      {
      if (attValue && *attValue == '\0')
        {
        this->SetBackgroundVolumeID(NULL);
        }
      else
        {
        this->SetBackgroundVolumeID(attValue);
        }
      }
    else if (!strcmp(attName, "foregroundVolumeID")) 
      {
      if (attValue && *attValue == '\0')
        {
        this->SetForegroundVolumeID(NULL);
        }
      else
        {
        this->SetForegroundVolumeID(attValue);
        }
      }
    else if (!strcmp(attName, "labelVolumeID")) 
      {
      if (attValue && *attValue == '\0')
        {
        this->SetLabelVolumeID(NULL);
        }
      else
        {
        this->SetLabelVolumeID(attValue);
        }
      }
    else if (!strcmp(attName, "foregroundOpacity")) 
      {
      this->SetForegroundOpacity( atof(attValue) );
      }
    else if (!strcmp(attName, "labelOpacity")) 
      {
      this->SetLabelOpacity( atof(attValue) );
      }
   else if (!strcmp(attName, "layoutName")) 
      {
      this->SetLayoutName( attValue );
      }
    }
}

//----------------------------------------------------------------------------
// Copy the node\"s attributes to this object.
// Does NOT copy: ID, FilePrefix, Name, SliceID
void vtkMRMLSliceCompositeNode::Copy(vtkMRMLNode *anode)
{
  Superclass::Copy(anode);
  vtkMRMLSliceCompositeNode *node = vtkMRMLSliceCompositeNode::SafeDownCast(anode);

  this->SetBackgroundVolumeID(node->GetBackgroundVolumeID());
  this->SetForegroundVolumeID(node->GetForegroundVolumeID());
  this->SetForegroundOpacity(node->GetForegroundOpacity());
  this->SetLabelOpacity(node->GetLabelOpacity());
}

//----------------------------------------------------------------------------
void vtkMRMLSliceCompositeNode::PrintSelf(ostream& os, vtkIndent indent)
{
  Superclass::PrintSelf(os,indent);

  os << indent << "BackgroundVolumeID: " << 
   (this->BackgroundVolumeID ? this->BackgroundVolumeID : "(none)") << "\n";
  os << indent << "ForegroundVolumeID: " << 
   (this->ForegroundVolumeID ? this->ForegroundVolumeID : "(none)") << "\n";
  os << indent << "LabelVolumeID: " << 
   (this->LabelVolumeID ? this->LabelVolumeID : "(none)") << "\n";
  os << indent << "ForegroundOpacity: " << this->ForegroundOpacity << "\n";
  os << indent << "LabelOpacity: " << this->LabelOpacity << "\n";
}

//----------------------------------------------------------------------------
void vtkMRMLSliceCompositeNode::UpdateScene(vtkMRMLScene* scene)
{
  vtkMRMLSliceCompositeNode *node= NULL;
  int nnodes = scene->GetNumberOfNodesByClass("vtkMRMLSliceCompositeNode");
  for (int n=0; n<nnodes; n++)
    {
    node = vtkMRMLSliceCompositeNode::SafeDownCast (
          scene->GetNthNodeByClass(n, "vtkMRMLSliceCompositeNode"));
    if (node != this && !strcmp(node->GetLayoutName(), this->GetLayoutName()))
      {
      break;
      }
    node = NULL;
    }
  if (node != NULL)
    {
    scene->RemoveNode(node);
    }
}
// End
