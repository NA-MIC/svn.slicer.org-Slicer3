/*=auto=========================================================================

Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

See Doc/copyright/copyright.txt
or http://www.slicer.org/copyright/copyright.txt for details.

Program:   3D Slicer
Module:    $RCSfile: vtkMRMLTumorGrowthNode.cxx,v $
Date:      $Date: 2006/03/17 15:10:10 $
Version:   $Revision: 1.2 $

=========================================================================auto=*/

#include <string>
#include <iostream>
#include <sstream>

#include "vtkObjectFactory.h"

#include "vtkMRMLTumorGrowthNode.h"
#include "vtkMRMLScene.h"


//------------------------------------------------------------------------------
vtkMRMLTumorGrowthNode* vtkMRMLTumorGrowthNode::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkMRMLTumorGrowthNode");
  if(ret)
    {
      return (vtkMRMLTumorGrowthNode*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkMRMLTumorGrowthNode;
}

//----------------------------------------------------------------------------

vtkMRMLNode* vtkMRMLTumorGrowthNode::CreateNodeInstance()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkMRMLTumorGrowthNode");
  if(ret)
    {
      return (vtkMRMLTumorGrowthNode*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkMRMLTumorGrowthNode;
}

//----------------------------------------------------------------------------
vtkMRMLTumorGrowthNode::vtkMRMLTumorGrowthNode()
{
   this->Conductance = 1.0;
   this->TimeStep = 0.1;
   this->HideFromEditors = true;

   this->Scan1_Ref = NULL;
   this->Scan2_Ref = NULL;
   this->Scan1_SuperSampleRef = NULL;
   this->Scan1_SegmentRef = NULL;

   // this->ROIMin[0] = this->ROIMin[1] = this->ROIMin[2] = this->ROIMax[0] = this->ROIMax[1] = this->ROIMax[2] = -1;
   this->ROIMin.resize(3,-1); 
   this->ROIMax.resize(3,-1); 

   this->SuperSampled_Spacing = -1;
   this->SuperSampled_VoxelVolume = -1;
   this->SuperSampled_RatioNewOldSpacing = -1;

   this->SegmentThreshold=-1;
}

//----------------------------------------------------------------------------
vtkMRMLTumorGrowthNode::~vtkMRMLTumorGrowthNode()
{
   this->SetScan1_Ref( NULL );
   this->SetScan2_Ref( NULL );
   this->SetScan1_SuperSampleRef( NULL);
   this->SetScan1_SegmentRef(NULL);
}

//----------------------------------------------------------------------------
void vtkMRMLTumorGrowthNode::WriteXML(ostream& of, int nIndent)
{
  Superclass::WriteXML(of, nIndent);

  // Write all MRML node attributes into output stream
  cout << "vtkMRMLTumorGrowthNode::WriteXML" << endl;
  vtkIndent indent(nIndent);

  {
    std::stringstream ss;
    ss << this->Conductance;
    of << indent << " Conductance=\"" << ss.str() << "\"";
  }
  {
    std::stringstream ss;
    ss << this->TimeStep;
    of << indent << " TimeStep=\"" << ss.str() << "\"";
  }
  {
    std::stringstream ss;
    if ( this->Scan1_Ref )
      {
      ss << this->Scan1_Ref;
      of << indent << " Scan1_Ref=\"" << ss.str() << "\"";
     }
  }
  {
    std::stringstream ss;
    if ( this->Scan2_Ref )
      {
      ss << this->Scan2_Ref;
      of << indent << " Scan2_Ref=\"" << ss.str() << "\"";
      }
  }

  of << indent << " ROIMin=\""<< this->ROIMin[0] << " "<< this->ROIMin[1] << " "<< this->ROIMin[2] <<"\"";
  of << indent << " ROIMax=\""<< this->ROIMax[0] << " "<< this->ROIMax[1] << " "<< this->ROIMax[2] <<"\"";

  // Do not write out the following parameters bc they are defined by rest
  // of << indent << " SuperSampled_Spacing=\""<< this->SuperSampled_Spacing  << "\"";
  // of << indent << " SuperSampled_VoxelVolume=\""<< this->SuperSampled_VoxelVolume  << "\"";
  // of << indent << " SuperSampled_RatioNewOldSpacing=\""<< this->SuperSampled_RatioNewOldSpacing  << "\"";

  of << indent << " SegmentThreshold=\""<< this->SegmentThreshold  << "\"";
}

//----------------------------------------------------------------------------
void vtkMRMLTumorGrowthNode::ReadXMLAttributes(const char** atts)
{
  vtkMRMLNode::ReadXMLAttributes(atts);

  // Read all MRML node attributes from two arrays of names and values
  const char* attName;
  const char* attValue;
  while (*atts != NULL) 
    {
    attName = *(atts++);
    attValue = *(atts++);
    if (!strcmp(attName, "Conductance")) 
      {
      std::stringstream ss;
      ss << attValue;
      ss >> this->Conductance;
      }
    else if (!strcmp(attName, "TimeStep")) 
      {
      std::stringstream ss;
      ss << attValue;
      ss >> this->TimeStep;
      }
    else if (!strcmp(attName, "Scan1_Ref"))
      {
      this->SetScan1_Ref(attValue);
      this->Scene->AddReferencedNodeID(this->Scan1_Ref, this);
      }
    else if (!strcmp(attName, "Scan2_Ref"))
      {
      this->SetScan2_Ref(attValue);
      this->Scene->AddReferencedNodeID(this->Scan2_Ref, this);
      }
    else if (!strcmp(attName, "ROIMin"))
      {
      // read data into a temporary vector
      vtksys_stl::stringstream ss;
      ss << attValue;
      ss >> this->ROIMin[0] >> this->ROIMin[1] >> this->ROIMin[2];
      }
    else if (!strcmp(attName, "ROIMax"))
      {
      // read data into a temporary vector
      vtksys_stl::stringstream ss;
      ss << attValue;
      ss >> this->ROIMax[0] >> this->ROIMax[1] >> this->ROIMax[2];
      }
    else if (!strcmp(attName, "SegmentThreshold"))
      {
    vtksys_stl::stringstream ss;
    ss << attValue;
    ss >>  this->SegmentThreshold; 
      }
     }
}


//----------------------------------------------------------------------------
// Copy the node's attributes to this object.
// Does NOT copy: ID, FilePrefix, Name, VolumeID
void vtkMRMLTumorGrowthNode::Copy(vtkMRMLNode *anode)
{
  Superclass::Copy(anode);
  vtkMRMLTumorGrowthNode *node = (vtkMRMLTumorGrowthNode *) anode;

  this->SetConductance(node->Conductance);
  this->SetTimeStep(node->TimeStep);
  this->SetScan1_Ref(node->Scan1_Ref);
  this->SetScan2_Ref(node->Scan2_Ref);
  this->ROIMin = node->ROIMin; 
  this->ROIMax = node->ROIMax; 
  this->SegmentThreshold = node->SegmentThreshold; 

}

//----------------------------------------------------------------------------
void vtkMRMLTumorGrowthNode::PrintSelf(ostream& os, vtkIndent indent)
{
  
  vtkMRMLNode::PrintSelf(os,indent);

  os << indent << "Conductance:       " << this->Conductance << "\n";
  os << indent << "TimeStep:          " << this->TimeStep << "\n";
  os << indent << "Scan1_Ref:      " << 
   (this->Scan1_Ref ? this->Scan1_Ref : "(none)") << "\n";
  os << indent << "OutputVolumeRef:   " << 
   (this->Scan2_Ref ? this->Scan2_Ref : "(none)") << "\n";
  os << indent << "ROIMin:            "<< this->ROIMin[0] << " "<< this->ROIMin[1] << " "<< this->ROIMin[2] <<"\n";
  os << indent << "ROIMax:            "<< this->ROIMax[0] << " "<< this->ROIMax[1] << " "<< this->ROIMax[2] <<"\n";
  os << indent << "SegmentThreshold:  "<< this->SegmentThreshold << "\n";

}

