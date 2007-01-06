/*=auto=========================================================================

  Portions (c) Copyright 2006 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See Doc/copyright/copyright.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: vtkMRMLColorTableNode.h,v $
  Date:      $Date: 2006/03/19 17:12:28 $
  Version:   $Revision: 1.0 $

=========================================================================auto=*/
// .NAME vtkMRMLColorTableNode - MRML node to represent discrete color information.
// .SECTION Description
// Color nodes describe colour look up tables. The tables may be pre-generated by
// Slicer (the label map colours, a random one) or created by
// a user. More than one model or label volume or editor can access the prebuilt
// nodes.

#ifndef __vtkMRMLColorTableNode_h
#define __vtkMRMLColorTableNode_h

#include <string>
#include <vector>

#include "vtkMRML.h"
#include "vtkMRMLNode.h"
#include "vtkMRMLColorNode.h"
#include "vtkLookupTable.h"

class vtkLookupTable;
class VTK_MRML_EXPORT vtkMRMLColorTableNode : public vtkMRMLColorNode
{
public:
  static vtkMRMLColorTableNode *New();
  vtkTypeMacro(vtkMRMLColorTableNode,vtkMRMLColorNode);
  void PrintSelf(ostream& os, vtkIndent indent);
  
  //--------------------------------------------------------------------------
  // MRMLNode methods
  //--------------------------------------------------------------------------

  virtual vtkMRMLNode* CreateNodeInstance();

  // Description:
  // Set node attributes
  virtual void ReadXMLAttributes( const char** atts);

  // Description:
  // Write this node's information to a MRML file in XML format.
  virtual void WriteXML(ostream& of, int indent);


  // Description:
  // Read in a text file holding colour table elements
  // id name r g b a
  // comments start with a hash mark
  virtual void ReadFile ();
  
  // Description:
  // Copy the node's attributes to this object
  virtual void Copy(vtkMRMLNode *node);
  
  // Description:
  // Get node XML tag name (like Volume, Model)
  virtual const char* GetNodeTagName() {return "ColorTable";};

  // Description:
  // 
  virtual void UpdateScene(vtkMRMLScene *scene);

  vtkGetObjectMacro(LookupTable, vtkLookupTable);
  vtkSetObjectMacro(LookupTable, vtkLookupTable);

  // Description:
  // Get/Set for Type
  void SetType(int type);
  vtkGetMacro(Type,int);
  void SetTypeToGrey();
  void SetTypeToIron();
  void SetTypeToRainbow();
  void SetTypeToOcean();
  void SetTypeToDesert();
  void SetTypeToInvGrey();
  void SetTypeToReverseRainbow();
  void SetTypeToFMRI();
  void SetTypeToFMRIPA();
  void SetTypeToLabels();
  void SetTypeToSPLBrainAtlas();
  void SetTypeToRandom();
  void SetTypeToUser();
  void SetTypeToFile();


  void ProcessMRMLEvents ( vtkObject *caller, unsigned long event, void *callData );

  //BTX
  // Description:
  // The list of valid table types
 
  // Grey - greyscale ramp
  // Iron - neutral
  // Rainbow - red-orange-yellow-blue-indigo-violet
  // Ocean - bluish ramp
  // Desert - orange ramp
  // InvGrey - inverted greyscale ramp
  // ReverseRainbow - inverted Rainbow
  // FMRI - fMRI map
  // FMRIPA - fMRI Positive Activation map
  // Labels - the Slicer2 default editor labels
  // SPLBrainAtlas - the SPL Brain atlas labels
  // Random - 255 random colors
  // User - user defined in the GUI
  // File - read in from file

  enum
    {
      Grey = 1,
      Iron = 2,
      Rainbow = 3,
      Ocean = 4,
      Desert = 5,
      InvGrey = 6,
      ReverseRainbow = 7,
      FMRI = 8,
      FMRIPA = 9,
      Labels = 10,
      SPLBrainAtlas = 11,
      Random = 12,
      User = 13,
      File = 14,
    };
  //ETX

  // Description:
  // Return the lowest and highest integers, for use in looping
  int GetFirstType () { return this->Grey; };
  int GetLastType () { return this->File; };
  
  // Description:
  // return a text string describing the colour look up table type
  virtual const char * GetTypeAsString();

  // Description:
  // Set the size of the colour table if it's a User table
  void SetNumberOfColors(int n);
  // Description:
  // Get the number of colours in the table
  int GetNumberOfColors();

  // Description:
  // keep track of where we last added a colour 
  int LastAddedColor;

  // Description:
  // Add a colour to the User colour table, at the end
  void AddColor(const char* name, double r, double g, double b);
  // Description:
  // Set a colour into the User colour table
  void SetColor(int entry, const char* name, double r, double g, double b);

protected:
  vtkMRMLColorTableNode();
  ~vtkMRMLColorTableNode();
  vtkMRMLColorTableNode(const vtkMRMLColorTableNode&);
  void operator=(const vtkMRMLColorTableNode&);

  // Description:
  // Set values in the names vector from the colour rgba entries in the colour
  // table
  virtual void SetNamesFromColors();

  // Description: 
  // The look up table, constructed according to the Type
  vtkLookupTable *LookupTable;
};

#endif
