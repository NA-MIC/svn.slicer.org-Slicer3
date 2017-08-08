#ifndef __vtkMRMLEMSWorkingDataNode_h
#define __vtkMRMLEMSWorkingDataNode_h

#include "vtkMRMLNode.h"
#include "vtkEMSegment.h"
#include "vtkMRMLScene.h"

class vtkMRMLEMSAtlasNode;
class vtkMRMLLabelMapVolumeNode;
class vtkMRMLEMSVolumeCollectionNode;

class VTK_EMSEGMENT_EXPORT vtkMRMLEMSWorkingDataNode :  public vtkMRMLNode 
{
public:
  static vtkMRMLEMSWorkingDataNode *New();
  vtkTypeMacro(vtkMRMLEMSWorkingDataNode,vtkMRMLNode);
  void PrintSelf(ostream& os, vtkIndent indent) VTK_OVERRIDE;

  virtual vtkMRMLNode* CreateNodeInstance() VTK_OVERRIDE;

  // Description:
  // Set node attributes
  virtual void ReadXMLAttributes(const char** atts) VTK_OVERRIDE;

  // Description:
  // Write this node's information to a MRML file in XML format.
  virtual void WriteXML(ostream& of, int indent) VTK_OVERRIDE;

  // Description:
  // Copy the node's attributes to this object
  virtual void Copy(vtkMRMLNode *node) VTK_OVERRIDE;

  // Description:
  // Get node XML tag name (like Volume, Model)
  virtual const char* GetNodeTagName() VTK_OVERRIDE {return "EMSWorkingData";}

  // Description:
  // Set the references of the node to the scene - only used in Slicer4 
  virtual void SetSceneReferences() VTK_OVERRIDE;

  // Description:
  // Updates this node if it depends on other nodes
  // when the node is deleted in the scene
  virtual void UpdateReferences() VTK_OVERRIDE;

  // Description:
  // Update the stored reference to another node in the scene
  virtual void UpdateReferenceID(const char *oldID, const char *newID) VTK_OVERRIDE;

  vtkGetStringMacro(InputTargetNodeID);
  //BTX
  vtkSetReferenceStringMacro(InputTargetNodeID);
  //ETX
  // For tcl Wrapping
  void SetReferenceInputTargetNodeID(const char *input)
  {
    this->SetInputTargetNodeID(input);
  }
  vtkMRMLEMSVolumeCollectionNode* GetInputTargetNode();

  vtkGetStringMacro(AlignedTargetNodeID);
  //BTX
  vtkSetReferenceStringMacro(AlignedTargetNodeID);
  //ETX
  void SetReferenceAlignedTargetNodeID(const char* name)
  {
    this->SetAlignedTargetNodeID(name);
  } 
  vtkMRMLEMSVolumeCollectionNode* GetAlignedTargetNode();
  
  vtkGetStringMacro(AlignedAtlasNodeID);
  //BTX
  vtkSetReferenceStringMacro(AlignedAtlasNodeID);
  //ETX
  void SetReferenceAlignedAtlasNodeID(const char* name)
  {
    this->SetAlignedAtlasNodeID(name);
  } 
  vtkMRMLEMSAtlasNode* GetAlignedAtlasNode();

  vtkGetStringMacro(AlignedSubParcellationNodeID);
  //BTX
  vtkSetReferenceStringMacro(AlignedSubParcellationNodeID);
  //ETX
  void SetReferenceAlignedSubParcellationNodeID(const char* name)
  {
    this->SetAlignedSubParcellationNodeID(name);
  } 
  vtkMRMLEMSVolumeCollectionNode* GetAlignedSubParcellationNode();

  vtkGetMacro(InputTargetNodeIsValid, int);
  vtkSetMacro(InputTargetNodeIsValid, int);

  vtkGetMacro(AlignedTargetNodeIsValid, int);
  vtkSetMacro(AlignedTargetNodeIsValid, int);

  vtkGetMacro(InputAtlasNodeIsValid, int);
  vtkSetMacro(InputAtlasNodeIsValid, int);

  vtkGetMacro(AlignedAtlasNodeIsValid, int);
  vtkSetMacro(AlignedAtlasNodeIsValid, int);

  vtkGetStringMacro(OutputSegmentationNodeID);
  vtkSetReferenceStringMacro(OutputSegmentationNodeID);
  vtkMRMLLabelMapVolumeNode* GetOutputSegmentationNode();

  // For legacy - to be complient with older mrml files  - please do not add functions for these variables 
  vtkGetStringMacro(InputAtlasNodeID);
  vtkGetStringMacro(InputSubParcellationNodeID);

  void RemoveInputSubParcellationNodeID()
    {
      this->SetInputSubParcellationNodeID(NULL);
    }

  void RemoveInputAtlasNodeID()
    {
      this->SetInputAtlasNodeID(NULL);
    }

protected:
  vtkMRMLEMSWorkingDataNode();
  ~vtkMRMLEMSWorkingDataNode();
  vtkMRMLEMSWorkingDataNode(const vtkMRMLEMSWorkingDataNode&);
  void operator=(const vtkMRMLEMSWorkingDataNode&);

  char*                InputTargetNodeID;
  char*                AlignedTargetNodeID;
  char*                AlignedAtlasNodeID;
  char*                AlignedSubParcellationNodeID;

  int                  InputTargetNodeIsValid;
  int                  AlignedTargetNodeIsValid;  
  int                  InputAtlasNodeIsValid;  
  int                  AlignedAtlasNodeIsValid;  

  char*                OutputSegmentationNodeID;

  // For legacy - to be complient with older mrml files 
  char*                InputAtlasNodeID;
  char*                InputSubParcellationNodeID;

  //BTX
  vtkSetReferenceStringMacro(InputAtlasNodeID);
  vtkSetReferenceStringMacro(InputSubParcellationNodeID);
  //ETX
  
};

#endif
