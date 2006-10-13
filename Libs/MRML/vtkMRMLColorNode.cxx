/*=auto=========================================================================

Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

See Doc/copyright/copyright.txt
or http://www.slicer.org/copyright/copyright.txt for details.

Program:   3D Slicer
Module:    $RCSfile: vtkMRMLColorNode.cxx,v $
Date:      $Date: 2006/03/03 22:26:39 $
Version:   $Revision: 1.0 $

=========================================================================auto=*/
#include <string>
#include <iostream>
#include <sstream>

#include "vtkObjectFactory.h"
#include "vtkCallbackCommand.h"

#include "vtkMRMLColorNode.h"
#include "vtkMRMLScene.h"

//------------------------------------------------------------------------------
vtkMRMLColorNode* vtkMRMLColorNode::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkMRMLColorNode");
  if(ret)
    {
    return (vtkMRMLColorNode*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkMRMLColorNode;
}

//-----------------------------------------------------------------------------

vtkMRMLNode* vtkMRMLColorNode::CreateNodeInstance()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkMRMLColorNode");
  if(ret)
    {
    return (vtkMRMLColorNode*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkMRMLColorNode;
}


//----------------------------------------------------------------------------
vtkMRMLColorNode::vtkMRMLColorNode()
{

  this->LookupTable = vtkLookupTable::New();
  this->SetTypeToGrey();
}

//----------------------------------------------------------------------------
vtkMRMLColorNode::~vtkMRMLColorNode()
{
    if (this->LookupTable)
    {
        this->LookupTable->Delete();
    }
  if (this->Name) {

      delete [] this->Name;
      this->Name = NULL;
  }
}

//----------------------------------------------------------------------------
void vtkMRMLColorNode::WriteXML(ostream& of, int nIndent)
{
  // Write all attributes not equal to their defaults
  
  Superclass::WriteXML(of, nIndent);
  
  vtkIndent indent(nIndent);
  
  of << " type=\"" << this->GetTypeAsString() << "\"";
  
  of << " color=\"" <<  "\"";
}

//----------------------------------------------------------------------------
void vtkMRMLColorNode::ReadXMLAttributes(const char** atts)
{

  Superclass::ReadXMLAttributes(atts);

  const char* attName;
  const char* attValue;
  while (*atts != NULL) 
  {
      attName = *(atts++);
      attValue = *(atts++);
      if (!strcmp(attName, "name"))
      {
          this->SetName(attValue);
      }
      else if (!strcmp(attName, "id"))
      {
          // handled at the vtkMRMLNode level
      }
      else  if (!strcmp(attName, "color")) 
      {
          std::stringstream ss;
          ss << attValue;
          //ss >> this->Color[0];
          //ss >> this->Color[1];
          //ss >> this->Color[2];
      }
      else if (!strcmp(attName, "type")) 
      {
          std::stringstream ss;
          ss << attValue;
          ss >> this->Type;
      }
      else
      {
          std::cerr << "Unknown attribute name " << attName << endl;
      }
  }
  vtkDebugMacro("Finished reading in xml attributes, list id = " << this->GetID() << " and name = " << this->GetName() << endl);
}


//----------------------------------------------------------------------------
// Copy the node's attributes to this object.
// Does NOT copy: ID, FilePrefix, Name, ID
void vtkMRMLColorNode::Copy(vtkMRMLNode *anode)
{
  Superclass::Copy(anode);
  vtkMRMLColorNode *node = (vtkMRMLColorNode *) anode;

  this->SetName(node->Name);
  this->SetLookupTable(node->LookupTable);
  this->SetType(node->Type);

}

//----------------------------------------------------------------------------
void vtkMRMLColorNode::PrintSelf(ostream& os, vtkIndent indent)
{
  
  Superclass::PrintSelf(os,indent);

  os << indent << "Name: " <<
      (this->Name ? this->Name : "(none)") << "\n";
  

  os << indent << "Type: (" << this->GetTypeAsString() << ")\n";

  if (this->LookupTable != NULL)
    {
    os << indent << "Look up table:\n";
    this->LookupTable->PrintSelf(os, indent.GetNextIndent());
    }
  
  if (this->Names.size() > 0)
    {
    os << indent << "Color Names:\n";
    for (int i = 0; i < this->Names.size(); i++)
      {
      os << indent << indent << i << " " << this->GetColorName(i) << endl;
      }
    }
}

//-----------------------------------------------------------

void vtkMRMLColorNode::UpdateScene(vtkMRMLScene *scene)
{
    Superclass::UpdateScene(scene);
    /*
    if (this->GetStorageNodeID() == NULL) 
    {
        //vtkErrorMacro("No reference StorageNodeID found");
        return;
    }

    vtkMRMLNode* mnode = scene->GetNodeByID(this->StorageNodeID);
    if (mnode) 
    {
        vtkMRMLStorageNode *node  = dynamic_cast < vtkMRMLStorageNode *>(mnode);
        node->ReadData(this);
        //this->SetAndObservePolyData(this->GetPolyData());
    }
    */
}

//----------------------------------------------------------------------------
void vtkMRMLColorNode::SetTypeToLabels()
{
    this->SetType(this->Labels);
}

//----------------------------------------------------------------------------
void vtkMRMLColorNode::SetTypeToRandom()
{
    this->SetType(this->Random);
}

//----------------------------------------------------------------------------
void vtkMRMLColorNode::SetTypeToOcean()
{
    this->SetType(this->Ocean);
}

//----------------------------------------------------------------------------
void vtkMRMLColorNode::SetTypeToGrey()
{
    this->SetType(this->Grey);
}

//----------------------------------------------------------------------------
void vtkMRMLColorNode::SetTypeToFMRIPA()
{
    this->SetType(this->FMRIPA);
}

//----------------------------------------------------------------------------
const char* vtkMRMLColorNode::GetTypeAsString()
{
  if (this->Type == this->Labels)
    {
    return "Labels";
    }
  if (this->Type == this->Random)
    {
    return "Random";
    }
  if (this->Type == this->Grey)
    {
    return "Grey";
    }
  if (this->Type == this->Ocean)
    {
    return "Ocean";
    }
  if (this->Type == this->FMRIPA)
    {
    return "fMRIPA";
    }
  return "(unknown)";
}

//---------------------------------------------------------------------------
void vtkMRMLColorNode::ProcessMRMLEvents ( vtkObject *caller,
                                           unsigned long event, 
                                           void *callData )
{
  Superclass::ProcessMRMLEvents(caller, event, callData);
/*
  vtkMRMLColorDisplayNode *dnode = this->GetDisplayNode();
  if (dnode != NULL && dnode == vtkMRMLColorDisplayNode::SafeDownCast(caller) &&
      event ==  vtkCommand::ModifiedEvent)
    {
        this->InvokeEvent(vtkMRMLColorNode::DisplayModifiedEvent, NULL);
    }
*/
  return;
}

//---------------------------------------------------------------------------
void vtkMRMLColorNode::SetType(int type)
{
    if (this->Type == type)
    {
        return;
    }
    vtkDebugMacro(<< this->GetClassName() << " (" << this << "): setting Type to " << type);
    this->Type = type;
    this->SetName(this->GetTypeAsString());
    
    if (this->Type == this->Random)
      {

      int size = 255;
      
      this->LookupTable->SetTableValue(0, 0, 0, 0, 0);
      this->LookupTable->SetRange(0, size);
      this->LookupTable->SetNumberOfColors(size + 1);
      for (int i = 1; i <= size; i++)
        {
        // table values have to be 0-1
        double r = (rand()%255)/255.0;
        double g = (rand()%255)/255.0;
        double b = (rand()%255)/255.0;
       
        this->LookupTable->SetTableValue(i, r, g, b, 1.0);
        }
      this->SetNamesFromColors();
      }

    if (this->Type == this->FMRIPA)
      {
      int size = 20;
      this->LookupTable->SetNumberOfTableValues(size);
      this->LookupTable->SetHueRange(0, 0.16667);
      this->LookupTable->SetSaturationRange(1, 1);
      this->LookupTable->SetValueRange(1, 1);
      this->LookupTable->SetRampToLinear();
      this->LookupTable->Build();
      this->SetNamesFromColors();
      }

    if (this->Type == this->Grey)
      {
      // from vtkSlicerSliceLayerLogic.cxx
      this->LookupTable->SetRampToLinear();
      this->LookupTable->SetTableRange(0, 255);
      this->LookupTable->SetHueRange(0, 0);
      this->LookupTable->SetSaturationRange(0, 0);
      this->LookupTable->SetValueRange(0, 1);
      this->LookupTable->SetAlphaRange(1, 1); // not used
      this->LookupTable->Build();
      this->SetNamesFromColors();
      }
    
    // invoke a modified event
    this->Modified();
    
    // invoke a type  modified event
    this->InvokeEvent(vtkMRMLColorNode::TypeModifiedEvent);
}

//---------------------------------------------------------------------------
void vtkMRMLColorNode::SetNamesFromColors()
{
  int size = this->LookupTable->GetNumberOfColors();
  double *rgba;
  // reset the names
  this->Names.clear();
  this->Names.resize(size);
  for (int i = 0; i < size; i++)
    {
    rgba = this->LookupTable->GetTableValue(i);
    std::stringstream ss;
    ss << "R=";
    ss << rgba[0];
    ss << " G=";
    ss << rgba[1];
    ss << " B=";
    ss << rgba[2];
    ss << " A=";
    ss << rgba[3];
    vtkDebugMacro("SetNamesFromColors: " << i << " Name = " << ss.str().c_str());
    this->SetColorName(i, ss.str().c_str());
    }
}

//---------------------------------------------------------------------------
const char *vtkMRMLColorNode::GetColorName(int ind)
{
  if (ind < this->Names.size() && ind >= 0)
    {
    return this->Names[ind].c_str();
    }
  else
    {
    return "invalid";
    }
}

//---------------------------------------------------------------------------
void vtkMRMLColorNode::SetColorName(int ind, const char *name)
{
  if (ind < this->Names.size() && ind >= 0)
    {
    this->Names[ind] = std::string(name);
    }
  else
    {
    std::cerr << "ERROR: SetColorName, index was out of bounds: " << ind << ", current size is " << this->Names.size() << endl;
    }
}

//---------------------------------------------------------------------------
void vtkMRMLColorNode::AddColorName(const char *name)
{
  this->Names.push_back(std::string(name));
}
