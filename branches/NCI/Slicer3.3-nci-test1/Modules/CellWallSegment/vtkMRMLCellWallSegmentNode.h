/*=auto=========================================================================

  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See Doc/copyright/copyright.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: vtkMRMLCellWallSegmentNode.h,v $
  Date:      $Date: 2006/03/19 17:12:29 $
  Version:   $Revision: 1.3 $

=========================================================================auto=*/
#ifndef __vtkMRMLCellWallSegmentNode_h
#define __vtkMRMLCellWallSegmentNode_h

#include "vtkMRML.h"
#include "vtkMRMLNode.h"
#include "vtkMRMLStorageNode.h"
#include "vtkMRMLFiducialListNode.h"

#include "vtkMatrix4x4.h"
#include "vtkTransform.h"
#include "vtkImageData.h"

#include "vtkCellWallSegment.h"

class vtkImageData;

class VTK_CellWallSegment_EXPORT vtkMRMLCellWallSegmentNode : public vtkMRMLNode
{
  public:
  static vtkMRMLCellWallSegmentNode *New();
  vtkTypeMacro(vtkMRMLCellWallSegmentNode,vtkMRMLNode);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Create instance of a GAD node.
  virtual vtkMRMLNode* CreateNodeInstance();

  // Description:
  // Set node attributes from name/value pairs
  virtual void ReadXMLAttributes( const char** atts);

  // Description:
  // Write this node's information to a MRML file in XML format.
  virtual void WriteXML(ostream& of, int indent);

  // Description:
  // Copy the node's attributes to this object
  virtual void Copy(vtkMRMLNode *node);

  // Description:
  // Get unique node XML tag name (like Volume, Model)
  virtual const char* GetNodeTagName() {return "GADParameters";};

  // Description:
  // Get/Set Number of iterations (module parameter)
  vtkGetMacro(NumberOfIterations, int);
  vtkSetMacro(NumberOfIterations, int);

  // Description:
  // Get/Set Conductance (module parameter)
  vtkGetMacro(Conductance, double);
  vtkSetMacro(Conductance, double);

  // Description:
  // Get/Set time step (module parameter)
  vtkGetMacro(TimeStep, double);
  vtkSetMacro(TimeStep, double);
 
  // Description:
  // Get/Set input volume MRML Id
  vtkGetStringMacro(InputVolumeRef);
  vtkSetStringMacro(InputVolumeRef);
  
  // Description:
  // Get/Set output volume MRML Id
  vtkGetStringMacro(OutputVolumeRef);
  vtkSetStringMacro(OutputVolumeRef);

  // Description:
    // Get/Set output volume MRML Id
    vtkGetStringMacro(SegmentationVolumeRef);
    vtkSetStringMacro(SegmentationVolumeRef);
    
  // Description:
  // Get/Set output volume MRML Id
  vtkGetStringMacro(FiducialListRef);
  vtkSetStringMacro(FiducialListRef);
  
  // Description:
  // Update the stored reference to another node in the scene
  virtual void UpdateReferenceID(const char *oldID, const char *newID);

 
protected:
  vtkMRMLCellWallSegmentNode();
  ~vtkMRMLCellWallSegmentNode();
  vtkMRMLCellWallSegmentNode(const vtkMRMLCellWallSegmentNode&);
  void operator=(const vtkMRMLCellWallSegmentNode&);

  double Conductance;
  double TimeStep;
  int NumberOfIterations;
  
  // reference a MRML fiducial list node, because the user will access the list and browse fiducials
  char* FiducialListRef;
  
  char* InputVolumeRef;
  char* OutputVolumeRef;
  char* SegmentationVolumeRef;

};

#endif

