/*=auto=========================================================================

  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See Doc/copyright/copyright.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: vtkMRMLVolumeDisplayNode.h,v $
  Date:      $Date: 2006/03/19 17:12:29 $
  Version:   $Revision: 1.3 $

=========================================================================auto=*/
// .NAME vtkMRMLVolumeDisplayNode - MRML node for representing a volume display attributes
// .SECTION Description
// vtkMRMLVolumeDisplayNode nodes describe how volume is displayed.

#ifndef __vtkMRMLVolumeDisplayNode_h
#define __vtkMRMLVolumeDisplayNode_h

#include "vtkMRML.h"
#include "vtkMRMLScene.h"
#include "vtkMRMLNode.h"
#include "vtkMRMLStorageNode.h"
#include "vtkMRMLColorNode.h"
#include "vtkMRMLDisplayNode.h"

#include "vtkMatrix4x4.h"
#include "vtkTransform.h"
#include "vtkImageData.h"

class vtkImageData;

class VTK_MRML_EXPORT vtkMRMLVolumeDisplayNode : public vtkMRMLDisplayNode
{
  public:
  static vtkMRMLVolumeDisplayNode *New();
  vtkTypeMacro(vtkMRMLVolumeDisplayNode,vtkMRMLDisplayNode);
  void PrintSelf(ostream& os, vtkIndent indent);

  virtual vtkMRMLNode* CreateNodeInstance();

  // Description:
  // Read node attributes from XML file
  virtual void ReadXMLAttributes( const char** atts);

  // Description:
  // Write this node's information to a MRML file in XML format.
  virtual void WriteXML(ostream& of, int indent);

  // Description:
  // Copy the node's attributes to this object
  virtual void Copy(vtkMRMLNode *node);

  // Description:
  // Get node XML tag name (like Volume, Model)
  virtual const char* GetNodeTagName() = 0;

  // Description:
  // Update the stored reference to another node in the scene
  virtual void UpdateReferenceID(const char *oldID, const char *newID);

  // Description:
  // Updates this node if it depends on other nodes 
  // when the node is deleted in the scene
  virtual void UpdateReferences();

  // Description:
  // Finds the storage node and read the data
  virtual void UpdateScene(vtkMRMLScene *scene);
  
  // Description:
  // Gets ImageData converted from the real data in the node
  virtual vtkImageData* GetImageData() {return NULL;};

  // Description:
  // Update the pipeline based on this node attributes
  virtual void UpdateImageDataPipeline() {};

  // Description:
  // String ID of the color MRML node
  void SetAndObserveColorNodeID(const char *ColorNodeID);
  //BTX
  void SetAndObserveColorNodeID(std::string ColorNodeID);
  //ETX
  vtkGetStringMacro(ColorNodeID);

  // Description:
  // Get associated color MRML node
  vtkMRMLColorNode* GetColorNode();

  // Description:
  // alternative method to propagate events generated in Display nodes
  virtual void ProcessMRMLEvents ( vtkObject * /*caller*/, 
                                   unsigned long /*event*/, 
                                   void * /*callData*/ );
  // Description:
  // set gray colormap
  virtual void SetDefaultColorMap();

protected:
  vtkMRMLVolumeDisplayNode();
  ~vtkMRMLVolumeDisplayNode();
  vtkMRMLVolumeDisplayNode(const vtkMRMLVolumeDisplayNode&);
  void operator=(const vtkMRMLVolumeDisplayNode&);

  char *ColorNodeID;

  vtkSetReferenceStringMacro(ColorNodeID);

  vtkMRMLColorNode *ColorNode;

};

#endif

