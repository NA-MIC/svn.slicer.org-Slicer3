/*=auto=========================================================================

Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

See Doc/copyright/copyright.txt
or http://www.slicer.org/copyright/copyright.txt for details.

Program:   3D Slicer
Module:    $RCSfile: vtkGradientAnisotropicDiffusionFilterLogic.cxx,v $
Date:      $Date: 2006/03/17 15:10:10 $
Version:   $Revision: 1.2 $

=========================================================================auto=*/

#include <string>
#include <iostream>
#include <sstream>

#include "vtkObjectFactory.h"

#include "vtkGradientAnisotropicDiffusionFilterLogic.h"
#include "vtkITKGradientAnisotropicDiffusionImageFilter.h"

#include "vtkMRMLScene.h"
#include "vtkMRMLScalarVolumeNode.h"

vtkGradientAnisotropicDiffusionFilterLogic* vtkGradientAnisotropicDiffusionFilterLogic::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkGradientAnisotropicDiffusionFilterLogic");
  if(ret)
    {
      return (vtkGradientAnisotropicDiffusionFilterLogic*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkGradientAnisotropicDiffusionFilterLogic;
}


//----------------------------------------------------------------------------
vtkGradientAnisotropicDiffusionFilterLogic::vtkGradientAnisotropicDiffusionFilterLogic()
{
  GradientAnisotropicDiffusionFilterNode = vtkMRMLGradientAnisotropicDiffusionFilterNode::New();
  GradientAnisotropicDiffusionImageFilter = vtkITKGradientAnisotropicDiffusionImageFilter::New();
}

//----------------------------------------------------------------------------
vtkGradientAnisotropicDiffusionFilterLogic::~vtkGradientAnisotropicDiffusionFilterLogic()
{
  this->GradientAnisotropicDiffusionImageFilter->Delete();
  this->GradientAnisotropicDiffusionFilterNode->Delete();
}

//----------------------------------------------------------------------------
void vtkGradientAnisotropicDiffusionFilterLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  
}

void vtkGradientAnisotropicDiffusionFilterLogic::Apply()
{
  // chack if MRML node is present 
  if (this->GradientAnisotropicDiffusionFilterNode == NULL)
    {
    vtkErrorMacro("No input GradientAnisotropicDiffusionFilterNode found");
    return;
    }
  
  // find input volume
  vtkMRMLNode* inNode = this->GetMRMLScene()->GetNodeByID(this->GradientAnisotropicDiffusionFilterNode->GetInputVolumeRef());
  vtkMRMLScalarVolumeNode *inVolume =  dynamic_cast<vtkMRMLScalarVolumeNode *> (inNode);
  if (inVolume == NULL)
    {
    vtkErrorMacro("No input volume found with id= " << this->GradientAnisotropicDiffusionFilterNode->GetInputVolumeRef());
    return;
    }
  
  this->GradientAnisotropicDiffusionImageFilter->SetInput(inVolume->GetImageData());
  
  
  // set filter parameters
  this->GradientAnisotropicDiffusionImageFilter->SetConductanceParameter(this->GradientAnisotropicDiffusionFilterNode->GetConductance());
  this->GradientAnisotropicDiffusionImageFilter->SetNumberOfIterations(this->GradientAnisotropicDiffusionFilterNode->GetNumberOfIterations());
  this->GradientAnisotropicDiffusionImageFilter->SetTimeStep(this->GradientAnisotropicDiffusionFilterNode->GetTimeStep());
  
  // find output volume
  vtkMRMLScalarVolumeNode *outVolume = NULL;
  if (this->GradientAnisotropicDiffusionFilterNode->GetOutputVolumeRef() != NULL)
    {
    vtkMRMLNode* outNode = this->GetMRMLScene()->GetNodeByID(this->GradientAnisotropicDiffusionFilterNode->GetOutputVolumeRef());
    outVolume =  dynamic_cast<vtkMRMLScalarVolumeNode *> (outNode);
    if (outVolume == NULL)
      {
      vtkErrorMacro("No output volume found with id= " << this->GradientAnisotropicDiffusionFilterNode->GetOutputVolumeRef());
      return;
      }
    }
  else 
    {
    // create new volume Node and add it to mrml scene
    this->GetMRMLScene()->SaveStateForUndo();
    outVolume = vtkMRMLScalarVolumeNode::New();
    this->GetMRMLScene()->AddNode(outVolume);  
    outVolume->Delete();
    }

  outVolume->SetImageData(this->GradientAnisotropicDiffusionImageFilter->GetOutput());
  this->GradientAnisotropicDiffusionImageFilter->Update();
}
