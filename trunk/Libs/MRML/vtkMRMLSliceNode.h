/*=auto=========================================================================

  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See Doc/copyright/copyright.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: vtkMRMLSliceNode.h,v $
  Date:      $Date: 2006/03/19 17:12:29 $
  Version:   $Revision: 1.3 $

=========================================================================auto=*/
// .NAME vtkMRMLSliceNode - MRML node for storing a slice through RAS space
// .SECTION Description
// This node stores the information about how to map from RAS space to 
// the desired slice plane.
// -- RASToSlice is the matrix that rotates and translates the slice plane
// -- FieldOfView tells the size of  slice plane
//

#ifndef __vtkMRMLSliceNode_h
#define __vtkMRMLSliceNode_h

#include "vtkMRML.h"
#include "vtkMRMLNode.h"

#include "vtkMatrix4x4.h"

class VTK_MRML_EXPORT vtkMRMLSliceNode : public vtkMRMLNode
{
  public:
  static vtkMRMLSliceNode *New();
  vtkTypeMacro(vtkMRMLSliceNode,vtkMRMLNode);
  void PrintSelf(ostream& os, vtkIndent indent);

  virtual vtkMRMLNode* CreateNodeInstance();

  // Description:
  // Set node attributes
  virtual void ReadXMLAttributes( const char** atts);

  // Description:
  // Write this node's information to a MRML file in XML format.
  virtual void WriteXML(ostream& of, int indent);

  // Description:
  // Copy the node's attributes to this object
  virtual void Copy(vtkMRMLNode *node);

  // Description:
  // Get node XML tag name (like Volume, Model)
  virtual const char* GetNodeTagName() {return "Slice";};

  // Description:
  // Mapping from RAS space onto the slice plane
  // TODO: maybe this should be a quaternion and a translate to avoid shears/scales
  vtkGetObjectMacro (RASToSlice, vtkMatrix4x4);
  vtkSetObjectMacro (RASToSlice, vtkMatrix4x4);

  // Description:
  // Size of the slice plane in 
  vtkGetVector3Macro (FieldOfView, double);
  vtkSetVector3Macro (FieldOfView, double);

  
protected:
  vtkMRMLSliceNode();
  ~vtkMRMLSliceNode();
  vtkMRMLSliceNode(const vtkMRMLSliceNode&);
  void operator=(const vtkMRMLSliceNode&);

  vtkMatrix4x4 *RASToSlice;
  double FieldOfView[3];

};

#endif

