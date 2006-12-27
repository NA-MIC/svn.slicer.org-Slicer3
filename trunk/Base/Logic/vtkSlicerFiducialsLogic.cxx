/*=auto=========================================================================

  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See Doc/copyright/copyright.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: vtkSlicerFiducialsLogic.cxx,v $
  Date:      $Date: 2006/01/06 17:56:48 $
  Version:   $Revision: 1.58 $

=========================================================================auto=*/

#include "vtkObjectFactory.h"
#include "vtkCallbackCommand.h"
#include <vtksys/SystemTools.hxx> 

#include "vtkSlicerFiducialsLogic.h"

#include "vtkMRMLFiducial.h"
#include "vtkMRMLFiducialListNode.h"
#include "vtkMRMLSelectionNode.h"

vtkCxxRevisionMacro(vtkSlicerFiducialsLogic, "$Revision: 1.9.12.1 $");
vtkStandardNewMacro(vtkSlicerFiducialsLogic);

//----------------------------------------------------------------------------
vtkSlicerFiducialsLogic::vtkSlicerFiducialsLogic()
{
 
}

//----------------------------------------------------------------------------
vtkSlicerFiducialsLogic::~vtkSlicerFiducialsLogic()
{
    
}

//----------------------------------------------------------------------------
void vtkSlicerFiducialsLogic::ProcessMRMLEvents()
{
  // TODO: implement if needed
}

//----------------------------------------------------------------------------
void vtkSlicerFiducialsLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->vtkObject::PrintSelf(os, indent);

  os << indent << "vtkSlicerFiducialsLogic:             " << this->GetClassName() << "\n";

}

//----------------------------------------------------------------------------
int vtkSlicerFiducialsLogic::AddFiducial(float x, float y, float z)
{
  // get the selection node
  vtkMRMLSelectionNode *selnode;
  selnode = vtkMRMLSelectionNode::SafeDownCast (
            this->MRMLScene->GetNthNodeByClass(0, "vtkMRMLSelectionNode"));
  
    if (selnode->GetActiveFiducialListID() == NULL)
      {
      vtkErrorMacro("FiducialsLogic: selection node doesn't have an active fiducial list right now, make one first before adding a fiducial");
      return -1;
      }
    // get the selected fiducial list
    vtkMRMLFiducialListNode *flist = vtkMRMLFiducialListNode::SafeDownCast(this->MRMLScene->GetNodeByID(selnode->GetActiveFiducialListID()));
    if (flist == NULL)
      {
      vtkErrorMacro("FiducialsLogic: selected fiducial list is null");
      return -1;
      }

    // add a fiducial
    int index = flist->AddFiducial();
    flist->SetNthFiducialXYZ(index, x, y, z);
    
    return index;
}
