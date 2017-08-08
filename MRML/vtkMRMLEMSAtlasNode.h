#ifndef __vtkMRMLEMSAtlasNode_h
#define __vtkMRMLEMSAtlasNode_h

#include "vtkMRML.h"
#include "vtkMRMLNode.h"
#include "vtkEMSegment.h"
#include "vtkMRMLEMSVolumeCollectionNode.h"

class VTK_EMSEGMENT_EXPORT vtkMRMLEMSAtlasNode : 
  public vtkMRMLEMSVolumeCollectionNode
{
public:
  static vtkMRMLEMSAtlasNode *New();
  vtkTypeMacro(vtkMRMLEMSAtlasNode,vtkMRMLEMSVolumeCollectionNode);
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
  virtual const char* GetNodeTagName() VTK_OVERRIDE {return "EMSAtlas";}

  vtkGetMacro(NumberOfTrainingSamples, int);
  vtkSetMacro(NumberOfTrainingSamples, int);

protected:
  vtkMRMLEMSAtlasNode();
  ~vtkMRMLEMSAtlasNode();
  vtkMRMLEMSAtlasNode(const vtkMRMLEMSAtlasNode&);
  void operator=(const vtkMRMLEMSAtlasNode&);

private:
  int NumberOfTrainingSamples;
};

#endif
