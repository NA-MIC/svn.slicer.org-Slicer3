/*=auto=========================================================================

  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See Doc/copyright/copyright.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: vtkMRMLScalarVolumeDisplayNode.h,v $
  Date:      $Date: 2006/03/19 17:12:29 $
  Version:   $Revision: 1.3 $

=========================================================================auto=*/
// .NAME vtkMRMLScalarVolumeDisplayNode - MRML node for representing a volume display attributes
// .SECTION Description
// vtkMRMLScalarVolumeDisplayNode nodes describe how volume is displayed.

#ifndef __vtkMRMLScalarVolumeDisplayNode_h
#define __vtkMRMLScalarVolumeDisplayNode_h

#include "vtkMRML.h"
#include "vtkMRMLScene.h"
#include "vtkMRMLNode.h"
#include "vtkMRMLVolumeDisplayNode.h"

#include "vtkImageData.h"

class vtkImageData;

class VTK_MRML_EXPORT vtkMRMLScalarVolumeDisplayNode : public vtkMRMLVolumeDisplayNode
{
  public:
  static vtkMRMLScalarVolumeDisplayNode *New();
  vtkTypeMacro(vtkMRMLScalarVolumeDisplayNode,vtkMRMLVolumeDisplayNode);
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
  virtual const char* GetNodeTagName() {return "VolumeDisplay";};


  //--------------------------------------------------------------------------
  // Display Information
  //--------------------------------------------------------------------------
  
  // Description:
  // Specifies whether windowing and leveling are to be performed automatically
  vtkBooleanMacro(AutoWindowLevel, int);
  vtkGetMacro(AutoWindowLevel, int);
  vtkSetMacro(AutoWindowLevel, int);

  // Description:
  // The window value to use when autoWindowLevel is 'no'
  vtkGetMacro(Window, double);
  vtkSetMacro(Window, double);

  // Description:
  // The level value to use when autoWindowLevel is 'no'
  vtkGetMacro(Level, double);
  vtkSetMacro(Level, double);

  // Description:
  // Specifies whether to apply the threshold
  vtkBooleanMacro(ApplyThreshold, int);
  vtkGetMacro(ApplyThreshold, int);
  vtkSetMacro(ApplyThreshold, int);

  // Description:
  // Specifies whether the threshold should be set automatically
  vtkBooleanMacro(AutoThreshold, int);
  vtkGetMacro(AutoThreshold, int);
  vtkSetMacro(AutoThreshold, int);

  // Description:
  // The upper threshold value to use when autoThreshold is 'no'
  vtkGetMacro(UpperThreshold, double);
  vtkSetMacro(UpperThreshold, double);

  // Description:
  // The lower threshold value to use when autoThreshold is 'no'
  vtkGetMacro(LowerThreshold, double);
  vtkSetMacro(LowerThreshold, double);

  // Description:
  // Set/Get interpolate reformated slices
  vtkGetMacro(Interpolate, int);
  vtkSetMacro(Interpolate, int);
  vtkBooleanMacro(Interpolate, int);

  virtual void SetDefaultColorMap();

  // Description:
  // alternative method to propagate events generated in Display nodes
  virtual void ProcessMRMLEvents ( vtkObject * /*caller*/, 
                                   unsigned long /*event*/, 
                                   void * /*callData*/ );

protected:
  vtkMRMLScalarVolumeDisplayNode();
  ~vtkMRMLScalarVolumeDisplayNode();
  vtkMRMLScalarVolumeDisplayNode(const vtkMRMLScalarVolumeDisplayNode&);
  void operator=(const vtkMRMLScalarVolumeDisplayNode&);

  double Window;
  double Level;
  double UpperThreshold;
  double LowerThreshold;


  // Booleans
  int Interpolate;
  int AutoWindowLevel;
  int ApplyThreshold;
  int AutoThreshold;

};

#endif

