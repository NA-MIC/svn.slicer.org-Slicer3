/*=auto=========================================================================

  Portions (c) Copyright 2006 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See Doc/copyright/copyright.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: vtkMRMLColorNode.h,v $
  Date:      $Date: 2006/03/19 17:12:28 $
  Version:   $Revision: 1.0 $

=========================================================================auto=*/
// .NAME vtkMRMLColorNode - MRML node to represent color information.
// .SECTION Description
// Color nodes describe colour look up tables. The tables may be pre-generated by
// Slicer (the label map colours, some default ramps, a random one) or created by
// a user. More than one model or label volume or editor can access the prebuilt
// nodes. This is used as a superclass for table based, procedural based, and
// implicit function based color nodes

#ifndef __vtkMRMLColorNode_h
#define __vtkMRMLColorNode_h

#include <string>
#include <vector>

#include "vtkMRML.h"
#include "vtkMRMLNode.h"

#include "vtkLookupTable.h"

class VTK_MRML_EXPORT vtkMRMLColorNode : public vtkMRMLNode
{
public:
  static vtkMRMLColorNode *New();
  vtkTypeMacro(vtkMRMLColorNode,vtkMRMLNode);
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
  // Read in a text file holding colours
  // Return 1 on sucess, 0 on failure
  virtual int ReadFile ();
  
  // Description:
  // Copy the node's attributes to this object
  virtual void Copy(vtkMRMLNode *node);
  
  // Description:
  // Get node XML tag name (like Volume, Model)
  virtual const char* GetNodeTagName() {return "Color";};

  // Description:
  // Reset node attributes to the initilal state as defined in the constructor.
  // NOTE:   it preserves values several dynamic attributes that may be set by an application: type, name
  virtual void Reset();
  
  // Description:
  // 
  virtual void UpdateScene(vtkMRMLScene *scene);

  // Description:
  // Set Type to type, then build colours and set names
  virtual void SetType(int type);
  // Description:
  // Get for Type
  vtkGetMacro(Type,int);

  void ProcessMRMLEvents ( vtkObject *caller, unsigned long event, void *callData );

  // Description:
  // Return the lowest and the highest type integers (defined in enum in
  // subclass), for use in looping
  virtual int GetFirstType();
  virtual int GetLastType ();
  
  // Description:
  // return a text string describing the colour look up table type
  virtual const char * GetTypeAsString();

  //BTX
  // Description:
  // TypeModifiedEvent is generated when the type of the colour look up table changes
  enum
    {
      TypeModifiedEvent = 20002,
    };
//ETX

  // Description:
  // Get the 0th based nth name of this colour
  const char *GetColorName(int ind);
  // Description:
  // Get the 0th based nth name of this colour, replacing the spaces with
  // subst
  //BTX
  std::string GetColorNameWithoutSpaces(int ind, const char *subst);
  //ETX
  
  // Description:
  // Add a color name to the vector
  void AddColorName(const char *name);
  // Description:
  // Set the 0th based nth name of this colour
  void SetColorName(int ind, const char *name);
  // Description:
  // Set the 0th based nth name of this colour, replacing the subst character
  // with spaces
  void SetColorNameWithSpaces(int ind, const char *name, const char *subst);
  
  // Description:
  // Name of the file name from which to read color information
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  // Description:
  // Most color nodes will implement a look up table, so provide a top level
  // get method
  virtual vtkLookupTable * GetLookupTable();

  // Description:
  // get/set the string used for an unnamed colour
  vtkGetStringMacro(NoName);
  vtkSetStringMacro(NoName);

  // Description:
  // Get/Set for the flag on names array having been initalised
  vtkGetMacro(NamesInitialised, int);
  vtkSetMacro(NamesInitialised, int);
  vtkBooleanMacro(NamesInitialised, int);
  
protected:
  vtkMRMLColorNode();
  virtual ~vtkMRMLColorNode();
  vtkMRMLColorNode(const vtkMRMLColorNode&);
  void operator=(const vtkMRMLColorNode&);

  // Description:
  // Set values in the names vector from the colours in the node
  virtual void SetNamesFromColors();
  
  // Description:
  // Which type of look up table does this node hold? 
  // Valid values are in the enumerated list
  int Type;

  //BTX
  // Description:
  // A vector of names for the color table elements
  std::vector<std::string> Names;
  //ETX

  // Description:
  // A file name to read text attributes from
  char *FileName;

  // Description:
  // the string used for an unnamed colour
  char *NoName;

  // Description:
  // Have the colour names been set? Used to do lazy copy of the Names array.
  int NamesInitialised;
};

#endif
