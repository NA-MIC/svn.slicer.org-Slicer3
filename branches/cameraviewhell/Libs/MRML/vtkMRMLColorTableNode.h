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

#include "vtkMRML.h"
#include "vtkMRMLNode.h"
#include "vtkMRMLColorNode.h"
#include "vtkMRMLColorTableStorageNode.h"
#include "vtkLookupTable.h"
#include "vtkMRMLScene.h"

class vtkLookupTable;
class vtkMRMLStorageNode;
class vtkMRMLScene;
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
  void SetTypeToFullRainbow();
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
  void SetTypeToRed();
  void SetTypeToGreen();
  void SetTypeToBlue();
  void SetTypeToCyan();
  void SetTypeToMagenta();
  void SetTypeToYellow();
  void SetTypeToWarm1();
  void SetTypeToWarm2();
  void SetTypeToWarm3();
  void SetTypeToCool1();
  void SetTypeToCool2();
  void SetTypeToCool3();  
  void SetTypeToWarmShade1();
  void SetTypeToWarmShade2();
  void SetTypeToWarmShade3();
  void SetTypeToCoolShade1();
  void SetTypeToCoolShade2();
  void SetTypeToCoolShade3();  
  void SetTypeToWarmTint1();
  void SetTypeToWarmTint2();
  void SetTypeToWarmTint3();
  void SetTypeToCoolTint1();
  void SetTypeToCoolTint2();
  void SetTypeToCoolTint3();  


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
  // Labels - the Slicer2 FullRainbow editor labels
  // SPLBrainAtlas - the SPL Brain atlas labels
  // Random - 255 random colors
  // User - user defined in the GUI
  // File - read in from file
  // Red - red ramp (like greyscale but with red, meant for layering with cyan)
  // Green - green ramp (like greyscale but with green, layering with magenta)
  // Blue - blue ramp (like greyscale but with blue, layering with yellow)
  // Yellow - yellow ramp (complementary ramp to blue, layering yeilds gray)
  // Cyan - cyan ramp (complementary ramp to red, layering yeilds gray)
  // Magenta - magenta ramp (complementary ramp to green, layering yeilds gray)
  // Warm# - ramps of warm colors that are complimentary to Cool#
  // WarmShade# - ramps of warm colors with variation in value that are
  //       complimentary to CoolShade# 
  // WarmTint# - ramps of warm colors with variation in saturation that are
  //       complimentary to CoolTint# 
  // Cool# - ramps of cool colors that are complimentary to Warm#
  // CoolShade# - ramps of cool colors with variation in value that are
  //       complimentary to WarmShade# 
  // CoolTint# - ramps of cool colors with variation in saturation that are
  //       complimentary to WarmSTint# 
  enum
    {
      FullRainbow = 0,
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
      Red = 15,
      Green = 16,
      Blue = 17,
      Yellow = 18,
      Cyan = 19,
      Magenta = 20,
      Warm1 = 21,
      Warm2 = 22,
      Warm3 = 23,
      Cool1 = 24,
      Cool2 = 25,
      Cool3 = 26,
      WarmShade1 = 27,
      WarmShade2 = 28,
      WarmShade3 = 29,
      CoolShade1 = 30,
      CoolShade2 = 31,
      CoolShade3 = 32,
      WarmTint1 = 33,
      WarmTint2 = 34,
      WarmTint3 = 35,
      CoolTint1 = 36,
      CoolTint2 = 37,
      CoolTint3 = 38,
    };
  //ETX

  // Description:
  // Return the lowest and highest integers, for use in looping
  int GetFirstType () { return this->FullRainbow; };
  int GetLastType () { return this->CoolTint3; };
  
  // Description:
  // return a text string describing the colour look up table type
  virtual const char * GetTypeAsString();

  // Description:
  // return a text string describing the class name and type name
  virtual const char * GetTypeAsIDString();

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

  // Description:
  // clear out the names list
  void ClearNames();

  // Description:
  // reset when close the scene
  virtual void Reset();

  // Description:
  // return the index associated with this color name, which can then be used
  // to get the colour. Returns -1 on failure.
  int GetColorIndexByName(const char *name);
 
  // Description:
  // Create default storage node or NULL if does not have one
  virtual vtkMRMLStorageNode* CreateDefaultStorageNode()
    {
    return vtkMRMLColorTableStorageNode::New();
    };

protected:
  vtkMRMLColorTableNode();
  virtual ~vtkMRMLColorTableNode();
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
