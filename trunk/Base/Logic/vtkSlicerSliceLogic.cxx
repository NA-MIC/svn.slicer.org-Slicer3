/*=auto=========================================================================

  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See Doc/copyright/copyright.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: vtkSlicerSliceLogic.cxx,v $
  Date:      $Date: 2006/01/06 17:56:48 $
  Version:   $Revision: 1.58 $

=========================================================================auto=*/

#include "vtkObjectFactory.h"
#include "vtkCallbackCommand.h"
#include "vtkPlaneSource.h"
#include "vtkPoints.h"

#include "vtkMRMLModelDisplayNode.h"

#include "vtkSlicerSliceLogic.h"

#include <sstream>

vtkCxxRevisionMacro(vtkSlicerSliceLogic, "$Revision: 1.9.12.1 $");
vtkStandardNewMacro(vtkSlicerSliceLogic);

//----------------------------------------------------------------------------
vtkSlicerSliceLogic::vtkSlicerSliceLogic()
{
  this->BackgroundLayer = NULL;
  this->ForegroundLayer = NULL;
  this->LabelLayer = NULL;
  this->SliceNode = NULL;
  this->SliceCompositeNode = NULL;
  this->ForegroundOpacity = 0.5; // Start by blending fg/bg
  this->LabelOpacity = 1.0;
  this->Blend = vtkImageBlend::New();
  this->SetForegroundOpacity(this->ForegroundOpacity);
  this->SetLabelOpacity(this->LabelOpacity);
  this->SliceModelNode = NULL;
  this->Name = NULL;
  this->SliceModelDisplayNode = NULL;
  this->ImageData = vtkImageData::New();
  this->SliceSpacing[0] = this->SliceSpacing[1] = this->SliceSpacing[2] = 1;
}

//----------------------------------------------------------------------------
vtkSlicerSliceLogic::~vtkSlicerSliceLogic()
{
  this->SetSliceNode(NULL);

  if (this->ImageData)
    {
    this->ImageData->Delete();
    }

  if ( this->Blend ) 
    {
    this->Blend->Delete();
    this->Blend = NULL;
    }

  this->SetBackgroundLayer (NULL);
  this->SetForegroundLayer (NULL);
  this->SetLabelLayer (NULL);

  if ( this->SliceCompositeNode ) 
    {
    vtkSetAndObserveMRMLNodeMacro( this->SliceCompositeNode, NULL );
    }

  this->SetName(NULL);

  this->DeleteSliceModel();

}

//----------------------------------------------------------------------------
void vtkSlicerSliceLogic::UpdateSliceNode()
{
  // find SliceNode in the scene
  vtkMRMLSliceNode *node= NULL;
  int nnodes = this->MRMLScene->GetNumberOfNodesByClass("vtkMRMLSliceNode");
  for (int n=0; n<nnodes; n++)
    {
    node = vtkMRMLSliceNode::SafeDownCast (
          this->MRMLScene->GetNthNodeByClass(n, "vtkMRMLSliceNode"));
    if (node->GetLayoutName() && !strcmp(node->GetLayoutName(), this->GetName()))
      {
      break;
      }
    node = NULL;
    }

  if ( this->SliceNode != NULL && node != NULL && 
    strcmp(this->SliceNode->GetID(), node->GetID()) != 0 )
    {
    // local SliceNode is out of sync with the scene
    this->SetSliceNode (NULL);
    }

  if ( this->SliceNode == NULL )
    {
    if ( node == NULL )
      {
      node = vtkMRMLSliceNode::New();
      node->SetLayoutName(this->GetName());
      //node->SetSingletonTag(this->GetName());
      this->SetSliceNode (node);
      this->UpdateSliceNodeFromLayout();
      node->Delete();
      }
    else
      {

      this->SetSliceNode (node);
      }
    }

  if ( this->MRMLScene->GetNodeByID(this->SliceNode->GetID()) == NULL)
    {
    // local node not in the scene
    node = this->SliceNode;
    node->Register(this);
    this->SetSliceNode (NULL);
    this->MRMLScene->AddNodeNoNotify(node);
    this->SetSliceNode (node);
    node->UnRegister(this);
    }
}

//----------------------------------------------------------------------------

void vtkSlicerSliceLogic::UpdateSliceNodeFromLayout()
{
  if (this->SliceNode == NULL)
    {
    return;
    }

  if ( !strcmp( this->GetName(), "Red" ) )
    {
    this->SliceNode->SetOrientationToAxial();
    }
  if ( !strcmp( this->GetName(), "Yellow" ) )
    {
    this->SliceNode->SetOrientationToSagittal();
    }
  if ( !strcmp( this->GetName(), "Green" ) )
    {
    this->SliceNode->SetOrientationToCoronal();
    }
}

//----------------------------------------------------------------------------

void vtkSlicerSliceLogic::UpdateSliceCompositeNode()
{
  // find SliceCompositeNode in the scene
  vtkMRMLSliceCompositeNode *node= NULL;
  int nnodes = this->MRMLScene->GetNumberOfNodesByClass("vtkMRMLSliceCompositeNode");
  for (int n=0; n<nnodes; n++)
    {
    node = vtkMRMLSliceCompositeNode::SafeDownCast (
          this->MRMLScene->GetNthNodeByClass(n, "vtkMRMLSliceCompositeNode"));
    if (node->GetLayoutName() && !strcmp(node->GetLayoutName(), this->GetName()))
      {
      break;
      }
    node = NULL;
    }

  if ( this->SliceCompositeNode != NULL && node != NULL && 
        strcmp(this->SliceCompositeNode->GetID(), node->GetID()) != 0 )
    {
    // local SliceCompositeNode is out of sync with the scene
    this->SetSliceCompositeNode (NULL);
    }

  if ( this->SliceCompositeNode == NULL )
    {
    if ( node == NULL )
      {
      node = vtkMRMLSliceCompositeNode::New();
      node->SetLayoutName(this->GetName());
      //node->SetSingletonTag(this->GetName());
      this->SetSliceCompositeNode (node);
      node->Delete();
      }
    else
      {
      this->SetSliceCompositeNode (node);
      }
    }

  if ( this->MRMLScene->GetNodeByID(this->SliceCompositeNode->GetID()) == NULL)
    {
    // local node not in the scene
    node = this->SliceCompositeNode;
    node->Register(this);
    this->SetSliceCompositeNode (NULL);
    this->MRMLScene->AddNodeNoNotify(node);
    this->SetSliceCompositeNode (node);
    node->UnRegister(this);
    }
}

//----------------------------------------------------------------------------
void vtkSlicerSliceLogic::ProcessMRMLEvents(vtkObject * caller, 
                                            unsigned long event, 
                                            void * callData )
{
 //
  // if you don't have a node yet, look in the scene to see if 
  // one exists for you to use.  If not, create one and add it to the scene
  //  
  if ( vtkMRMLScene::SafeDownCast(caller) == this->MRMLScene 
    && (event == vtkMRMLScene::NodeAddedEvent || event == vtkMRMLScene::NodeRemovedEvent ) )
    {
    vtkMRMLNode *node =  reinterpret_cast<vtkMRMLNode*> (callData);
    if (node == NULL || !(node->IsA("vtkMRMLSliceCompositeNode") || node->IsA("vtkMRMLSliceNode") || node->IsA("vtkMRMLVolumeNode")) )
      {
      return;
      }
    }

  if (event == vtkMRMLScene::SceneCloseEvent) 
    {
    this->UpdateSliceNodeFromLayout();
    this->DeleteSliceModel();

    return;
    }

  this->CreateSliceModel();

  this->UpdateSliceNode();
 
  this->UpdateSliceCompositeNode();

  //
  // check that our referenced nodes exist, and if not set to None
  //

  /** PROBABLY DONT NEED TO DO THAT And it causes load scene to override ID's
  if ( this->MRMLScene->GetNodeByID( this->SliceCompositeNode->GetForegroundVolumeID() ) == NULL )
    {
    this->SliceCompositeNode->SetForegroundVolumeID(NULL);
    }

  if ( this->MRMLScene->GetNodeByID( this->SliceCompositeNode->GetLabelVolumeID() ) == NULL )
    {
    this->SliceCompositeNode->SetLabelVolumeID(NULL);
    }

  if ( this->MRMLScene->GetNodeByID( this->SliceCompositeNode->GetBackgroundVolumeID() ) == NULL )
    {
    this->SliceCompositeNode->SetBackgroundVolumeID(NULL);
    }
   **/

  if (event != vtkMRMLScene::NewSceneEvent) 
    {
    this->UpdatePipeline();
    }
}

//----------------------------------------------------------------------------
void vtkSlicerSliceLogic::ProcessLogicEvents()
{
  //
  // if we don't have layers yet, create them 
  //
  if ( this->BackgroundLayer == NULL )
    {
    vtkSlicerSliceLayerLogic *layer = vtkSlicerSliceLayerLogic::New();
    this->SetBackgroundLayer (layer);
    layer->Delete();
    }
  if ( this->ForegroundLayer == NULL )
    {
    vtkSlicerSliceLayerLogic *layer = vtkSlicerSliceLayerLogic::New();
    this->SetForegroundLayer (layer);
    layer->Delete();
    }
  if ( this->LabelLayer == NULL )
    {
    vtkSlicerSliceLayerLogic *layer = vtkSlicerSliceLayerLogic::New();
    this->SetLabelLayer (layer);
    layer->Delete();
    }

  // Update slice plane geometry
  if (this->SliceNode != NULL && this->GetSliceModelNode() != NULL 
      && this->MRMLScene->GetNodeByID( this->SliceModelNode->GetID() ) != NULL && 
        this->SliceModelNode->GetPolyData() != NULL )
    {
    vtkPoints *points = this->SliceModelNode->GetPolyData()->GetPoints();
    unsigned int *dims = this->SliceNode->GetDimensions();
    vtkMatrix4x4 *xyToRAS = this->SliceNode->GetXYToRAS();

    double inPt[4]={0,0,0,1};
    double outPt[4];
    double *outPt3 = outPt;

    xyToRAS->MultiplyPoint(inPt, outPt);
    points->SetPoint(0, outPt3);

    inPt[0] = dims[0];
    xyToRAS->MultiplyPoint(inPt, outPt);
    points->SetPoint(1, outPt3);

    inPt[0] = 0;
    inPt[1] = dims[1];
    xyToRAS->MultiplyPoint(inPt, outPt);
    points->SetPoint(2, outPt3);

    inPt[0] = dims[0];
    inPt[1] = dims[1];
    xyToRAS->MultiplyPoint(inPt, outPt);
    points->SetPoint(3, outPt3);

    this->UpdateImageData();
    this->SliceModelNode->GetPolyData()->Modified();
    vtkMRMLModelDisplayNode *modelDisplayNode = this->SliceModelNode->GetDisplayNode();
    if ( modelDisplayNode )
      {
      modelDisplayNode->SetVisibility( this->SliceNode->GetSliceVisible() );
      }
    }

  // This is called when a slice layer is modified, so pass it on
  // to anyone interested in changes to this sub-pipeline
  this->Modified();
}

//----------------------------------------------------------------------------
void vtkSlicerSliceLogic::SetSliceNode(vtkMRMLSliceNode *sliceNode)
{
    // Don't directly observe the slice node -- the layers will observe it and
    // will notify us when things have changed.
    // This class takes care of passing the one slice node to each of the layers
    // so that users of this class only need to set the node one place.
  vtkSetMRMLNodeMacro( this->SliceNode, sliceNode );

  if (this->BackgroundLayer)
    {
    this->BackgroundLayer->SetSliceNode(sliceNode);
    }
  if (this->ForegroundLayer)
    {
    this->ForegroundLayer->SetSliceNode(sliceNode);
    }
  if (this->LabelLayer)
    {
    this->LabelLayer->SetSliceNode(sliceNode);
    }

  this->Modified();

}

//----------------------------------------------------------------------------
void vtkSlicerSliceLogic::SetSliceCompositeNode(vtkMRMLSliceCompositeNode *sliceCompositeNode)
{
    // Observe the composite node, since this holds the parameters for
    // this pipeline
  vtkSetAndObserveMRMLNodeMacro( this->SliceCompositeNode, sliceCompositeNode );
  this->UpdatePipeline();

}

//----------------------------------------------------------------------------
void vtkSlicerSliceLogic::SetBackgroundLayer(vtkSlicerSliceLayerLogic *BackgroundLayer)
{
  if (this->BackgroundLayer)
    {
    this->BackgroundLayer->RemoveObserver( this->LogicCallbackCommand );
    this->BackgroundLayer->SetAndObserveMRMLScene( NULL );
    this->BackgroundLayer->Delete();
    }
  this->BackgroundLayer = BackgroundLayer;

  if (this->BackgroundLayer)
    {
    this->BackgroundLayer->Register(this);

    vtkIntArray *events = vtkIntArray::New();
    events->InsertNextValue(vtkMRMLScene::NewSceneEvent);
    events->InsertNextValue(vtkMRMLScene::SceneCloseEvent);
    events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
    events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
    this->BackgroundLayer->SetAndObserveMRMLSceneEvents ( this->MRMLScene, events );
    events->Delete();

    this->BackgroundLayer->SetSliceNode(SliceNode);
    this->BackgroundLayer->AddObserver( vtkCommand::ModifiedEvent, this->LogicCallbackCommand );
    }

  this->Modified();
}

//----------------------------------------------------------------------------
void vtkSlicerSliceLogic::SetForegroundLayer(vtkSlicerSliceLayerLogic *ForegroundLayer)
{
  if (this->ForegroundLayer)
    {
    this->ForegroundLayer->RemoveObserver( this->LogicCallbackCommand );
    this->ForegroundLayer->SetAndObserveMRMLScene( NULL );
    this->ForegroundLayer->Delete();
    }
  this->ForegroundLayer = ForegroundLayer;

  if (this->ForegroundLayer)
    {
    this->ForegroundLayer->Register(this);

    vtkIntArray *events = vtkIntArray::New();
    events->InsertNextValue(vtkMRMLScene::NewSceneEvent);
    events->InsertNextValue(vtkMRMLScene::SceneCloseEvent);
    events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
    events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
    this->ForegroundLayer->SetAndObserveMRMLSceneEvents ( this->MRMLScene, events );
    events->Delete();

    this->ForegroundLayer->SetSliceNode(SliceNode);
    this->ForegroundLayer->AddObserver( vtkCommand::ModifiedEvent, this->LogicCallbackCommand );
    }

  this->Modified();
}

//----------------------------------------------------------------------------
void vtkSlicerSliceLogic::SetLabelLayer(vtkSlicerSliceLayerLogic *LabelLayer)
{
  if (this->LabelLayer)
    {
    this->LabelLayer->RemoveObserver( this->LogicCallbackCommand );
    this->LabelLayer->SetAndObserveMRMLScene( NULL );
    this->LabelLayer->Delete();
    }
  this->LabelLayer = LabelLayer;

  if (this->LabelLayer)
    {
    this->LabelLayer->Register(this);

    vtkIntArray *events = vtkIntArray::New();
    events->InsertNextValue(vtkMRMLScene::NewSceneEvent);
    events->InsertNextValue(vtkMRMLScene::SceneCloseEvent);
    events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
    events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
    this->LabelLayer->SetAndObserveMRMLSceneEvents ( this->MRMLScene, events );
    events->Delete();

    this->LabelLayer->SetSliceNode(SliceNode);
    this->LabelLayer->AddObserver( vtkCommand::ModifiedEvent, this->LogicCallbackCommand );
    }

  this->Modified();
}


//----------------------------------------------------------------------------
void vtkSlicerSliceLogic::SetForegroundOpacity(double ForegroundOpacity)
{
  this->ForegroundOpacity = ForegroundOpacity;

  if ( this->Blend->GetOpacity(1) != this->ForegroundOpacity )
    {
    this->Blend->SetOpacity(1, this->ForegroundOpacity);
    this->Modified();
    }
}

//----------------------------------------------------------------------------
void vtkSlicerSliceLogic::SetLabelOpacity(double LabelOpacity)
{
  this->LabelOpacity = LabelOpacity;

  if ( this->Blend->GetOpacity(2) != this->LabelOpacity )
    {
    this->Blend->SetOpacity(2, this->LabelOpacity);
    this->Modified();
    }
}

//----------------------------------------------------------------------------
void vtkSlicerSliceLogic::UpdatePipeline()
{

  int modified = 0;

  if ( this->SliceCompositeNode )
    {
    // get the background and foreground image data from the layers
    // so we can use them as input to the image blend
    // TODO: change logic to use a volume node superclass rather than
    // a scalar volume node once the superclass is sorted out for vector/tensor Volumes

    const char *id;
    
    // Background
    id = this->SliceCompositeNode->GetBackgroundVolumeID();
    vtkMRMLVolumeNode *bgnode = NULL;
    if (id)
      {
      bgnode = vtkMRMLVolumeNode::SafeDownCast (this->MRMLScene->GetNodeByID(id));
      }
    
    if (this->BackgroundLayer)
      {
      if ( this->BackgroundLayer->GetVolumeNode() != bgnode ) 
        {
        this->BackgroundLayer->SetVolumeNode (bgnode);
        modified = 1;
        }
      }

    // Foreground
    id = this->SliceCompositeNode->GetForegroundVolumeID();
    vtkMRMLVolumeNode *fgnode = NULL;
    if (id)
      {
      fgnode = vtkMRMLVolumeNode::SafeDownCast (this->MRMLScene->GetNodeByID(id));
      }
    
    if (this->ForegroundLayer)
      {
      if ( this->ForegroundLayer->GetVolumeNode() != fgnode ) 
        {
        this->ForegroundLayer->SetVolumeNode (fgnode);
        modified = 1;
        }
      }

    // Label
    id = this->SliceCompositeNode->GetLabelVolumeID();
    vtkMRMLVolumeNode *lbnode = NULL;
    if (id)
      {
      lbnode = vtkMRMLVolumeNode::SafeDownCast (this->MRMLScene->GetNodeByID(id));
      }
    
    if (this->LabelLayer)
      {
      if ( this->LabelLayer->GetVolumeNode() != lbnode ) 
        {
        this->LabelLayer->SetVolumeNode (lbnode);
        modified = 1;
        }
      }

    // Now update the image blend with the background and foreground and label
    // -- layer 0 opacity is ignored, but since not all inputs may be non-null, 
    //    we keep track so that someone could, for example, have a NULL background
    //    with a non-null foreground and label and everything will work with the 
    //    label opacity
    //
    // -- first make a temp blend instance and set it up according to the current 
    //    parameters.  Then check if this is the same as the current 'real' blend,
    //    and if not send a modified event
    //
    vtkImageBlend *tempBlend = vtkImageBlend::New();
    
    tempBlend->RemoveAllInputs ( );
    int layerIndex = 0;
    if ( this->BackgroundLayer && this->BackgroundLayer->GetImageData() )
      {
      tempBlend->AddInput( this->BackgroundLayer->GetImageData() );
      tempBlend->SetOpacity( layerIndex++, 1.0 );
      }
    if ( this->ForegroundLayer && this->ForegroundLayer->GetImageData() )
      {
      tempBlend->AddInput( this->ForegroundLayer->GetImageData() );
      tempBlend->SetOpacity( layerIndex++, this->SliceCompositeNode->GetForegroundOpacity() );
      }
    if ( this->LabelLayer && this->LabelLayer->GetImageData() )
      {
      tempBlend->AddInput( this->LabelLayer->GetImageData() );
      tempBlend->SetOpacity( layerIndex++, this->SliceCompositeNode->GetLabelOpacity() );
      }

    for (layerIndex = 0; layerIndex < tempBlend->GetNumberOfInputs(); layerIndex++)
      {
      if ( tempBlend->GetInput(layerIndex) != this->Blend->GetInput(layerIndex) )
        {
        this->Blend->SetInput(layerIndex, tempBlend->GetInput(layerIndex));
        modified = 1;
        }
      if ( tempBlend->GetOpacity(layerIndex) != this->Blend->GetOpacity(layerIndex) )
        {
        this->Blend->SetOpacity(layerIndex, tempBlend->GetOpacity(layerIndex));
        modified = 1;
        }
      }

    tempBlend->Delete();


    if ( this->SliceModelNode && 
          this->SliceModelNode->GetDisplayNode() &&
            this->SliceNode ) 
      {
      if (this->SliceModelNode->GetDisplayNode()->GetVisibility() != this->SliceNode->GetSliceVisible() )
        {
        this->SliceModelNode->GetDisplayNode()->SetVisibility( this->SliceNode->GetSliceVisible() );
        }
      if (this->SliceModelNode->GetDisplayNode()->GetTextureImageData() != this->GetImageData())
        {
        this->SliceModelNode->GetDisplayNode()->SetAndObserveTextureImageData(this->GetImageData());
        }
      }

    this->UpdateImageData();

    if ( modified )
      {
      this->Modified();
      }
    }
}

//----------------------------------------------------------------------------
void vtkSlicerSliceLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->vtkObject::PrintSelf(os, indent);

  os << indent << "SlicerSliceLogic:             " << this->GetClassName() << "\n";

  os << indent << "SliceNode: " <<
    (this->SliceNode ? this->SliceNode->GetName() : "(none)") << "\n";
  os << indent << "SliceCompositeNode: " <<
    (this->SliceCompositeNode ? this->SliceCompositeNode->GetName() : "(none)") << "\n";
  // TODO: fix printing of vtk objects
  os << indent << "BackgroundLayer: " <<
    (this->BackgroundLayer ? "not null" : "(none)") << "\n";
  os << indent << "ForegroundLayer: " <<
    (this->ForegroundLayer ? "not null" : "(none)") << "\n";
  os << indent << "LabelLayer: " <<
    (this->LabelLayer ? "not null" : "(none)") << "\n";
  os << indent << "Blend: " <<
    (this->Blend ? "not null" : "(none)") << "\n";
  os << indent << "ForegroundOpacity: " << this->ForegroundOpacity << "\n";
  os << indent << "LabelOpacity: " << this->LabelOpacity << "\n";

}

//----------------------------------------------------------------------------
void vtkSlicerSliceLogic::DeleteSliceModel()
{
  if (this->SliceModelNode != NULL)
    {
    this->SliceModelNode->SetAndObserveDisplayNodeID(NULL);
    this->SliceModelNode->SetAndObservePolyData(NULL);
    this->SliceModelNode->Delete();
    this->SliceModelNode = NULL;
    }
  if (this->SliceModelDisplayNode != NULL)
    {
    this->SliceModelDisplayNode->SetAndObserveTextureImageData(NULL);
    this->SliceModelDisplayNode->Delete();
    this->SliceModelDisplayNode = NULL;
    }
}

//----------------------------------------------------------------------------
void vtkSlicerSliceLogic::CreateSliceModel()
{

  if (this->SliceModelNode != NULL && this->MRMLScene->GetNodeByID( this->GetSliceModelNode()->GetID() ) == NULL )
    {
    this->DeleteSliceModel();
    }

  if ( this->SliceModelNode == NULL) 
  {
    this->SliceModelNode = vtkMRMLModelNode::New();
    this->SliceModelNode->SetScene(this->GetMRMLScene());
    this->SliceModelNode->SetHideFromEditors(1);
    this->SliceModelNode->SetSelectable(0);
    this->SliceModelNode->SetSaveWithScene(0);

    // create plane slice
    vtkPlaneSource *planeSource;
    planeSource = vtkPlaneSource::New();
    planeSource->GetOutput()->Update();
    this->SliceModelNode->SetAndObservePolyData(planeSource->GetOutput());
    planeSource->Delete();

    // create display node and set texture
    this->SliceModelDisplayNode = vtkMRMLModelDisplayNode::New();
    this->SliceModelDisplayNode->SetScene(this->GetMRMLScene());
    this->SliceModelDisplayNode->SetVisibility(0);
    this->SliceModelDisplayNode->SetOpacity(1);
    this->SliceModelDisplayNode->SetColor(1,1,1);
    this->SliceModelDisplayNode->SetAmbient(1);
    this->SliceModelDisplayNode->SetDiffuse(0);
    this->SliceModelDisplayNode->SetAndObserveTextureImageData(this->GetImageData());
    this->SliceModelDisplayNode->SetSaveWithScene(0);

    std::stringstream ss;
    char name[256];
    ss << this->Name << " Volume Slice";
    ss.getline(name,256);
    this->SliceModelNode->SetName(name);


    this->SliceModelNode->SetAndObserveDisplayNodeID(this->SliceModelDisplayNode->GetID());
  }

  if (this->SliceModelNode != NULL && this->MRMLScene->GetNodeByID( this->GetSliceModelNode()->GetID() ) == NULL )
    {
    vtkIntArray *events = vtkIntArray::New();
    //this->GetMRMLEvents(this->GetMRMLScene(), events);
    //this->SetMRMLScene(this->GetMRMLScene());
    this->MRMLScene->AddNodeNoNotify(this->SliceModelDisplayNode);
    this->MRMLScene->AddNode(this->SliceModelNode);
    this->SliceModelNode->SetAndObserveDisplayNodeID(this->SliceModelDisplayNode->GetID());
    this->SliceModelDisplayNode->SetAndObserveTextureImageData(this->GetImageData());
    //this->SetAndObserveMRMLSceneEvents(this->GetMRMLScene(), events);
    events->Delete();
    //this->SliceModelNode->Delete();
    //this->SliceModelDisplayNode->Delete();
    }

  // update the description to refer back to the slice and composite nodes
  // TODO: this doesn't need to be done unless the ID change, but it needs
  // to happen after they have been set, so do it every event for now
  if ( this->SliceModelNode != NULL ) {
    char description[256];
    std::stringstream ssD;
    vtkMRMLSliceNode *sliceNode = this->GetSliceNode();
    if ( sliceNode )
      {
      ssD << " SliceID " << sliceNode->GetID();
      }
    vtkMRMLSliceCompositeNode *compositeNode = this->GetSliceCompositeNode();
    if ( compositeNode )
      {
      ssD << " CompositeID " << compositeNode->GetID();
      }

    ssD.getline(description,256);
    this->SliceModelNode->SetDescription(description);
  }
}

// Get the size of the volume, transformed to RAS space
void vtkSlicerSliceLogic::GetBackgroundRASDimensions(double rasDimensions[3], double rasCenter[3])
{
  rasCenter[0] = rasDimensions[0] = 0.0;
  rasCenter[1] = rasDimensions[1] = 0.0;
  rasCenter[2] = rasDimensions[2] = 0.0;

  vtkMRMLSliceNode *sliceNode = this->GetSliceNode();
  vtkMRMLSliceCompositeNode *compositeNode = this->GetSliceCompositeNode();
  vtkMRMLVolumeNode *backgroundNode = NULL;

  if ( !sliceNode || !compositeNode )
    {
    return;
    }
  
  backgroundNode = vtkMRMLVolumeNode::SafeDownCast (
      this->MRMLScene->GetNodeByID( compositeNode->GetBackgroundVolumeID() ));

  vtkImageData *backgroundImage;
  if ( !backgroundNode || ! (backgroundImage = backgroundNode->GetImageData()) )
    {
    return;
    }

  //
  // Get the size of the volume in RAS space
  // - map the size of the volume in IJK into RAS
  // - map the middle of the volume to RAS for the center
  //   (IJK space always has origin at first pixel)
  //
  vtkMatrix4x4 *ijkToRAS = vtkMatrix4x4::New();
  int dimensions[3];
  backgroundImage->GetDimensions(dimensions);
  double doubleDimensions[4], rasHDimensions[4], rasHCenter[4];
  doubleDimensions[0] = dimensions[0];
  doubleDimensions[1] = dimensions[1];
  doubleDimensions[2] = dimensions[2];
  doubleDimensions[3] = 0;
  backgroundNode->GetIJKToRASMatrix (ijkToRAS);
  ijkToRAS->MultiplyPoint( doubleDimensions, rasHDimensions );
  doubleDimensions[0] = dimensions[0]/2.;
  doubleDimensions[1] = dimensions[1]/2.;
  doubleDimensions[2] = dimensions[2]/2.;
  doubleDimensions[3] = 1.;
  ijkToRAS->MultiplyPoint( doubleDimensions, rasHCenter );
  ijkToRAS->Delete();

  // ignore homogeneous coordinate
  rasDimensions[0] = rasHDimensions[0];
  rasDimensions[1] = rasHDimensions[1];
  rasDimensions[2] = rasHDimensions[2];
  rasCenter[0] = rasHCenter[0];
  rasCenter[1] = rasHCenter[1];
  rasCenter[2] = rasHCenter[2];

}

// Get the size of the volume, transformed to RAS space
void vtkSlicerSliceLogic::GetBackgroundSliceDimensions(double sliceDimensions[3], double sliceCenter[3])
{
  sliceCenter[0] = sliceDimensions[0] = 0.0;
  sliceCenter[1] = sliceDimensions[1] = 0.0;
  sliceCenter[2] = sliceDimensions[2] = 0.0;

  vtkMRMLSliceNode *sliceNode = this->GetSliceNode();

  if ( !sliceNode )
    {
    return;
    }
  
  double rasDimensions[3], rasHDimensions[4], sliceHDimensions[4];
  double rasCenter[3], rasHCenter[4], sliceHCenter[4];
  this->GetBackgroundRASDimensions(rasDimensions, rasCenter);
  rasHDimensions[0] = rasDimensions[0];
  rasHDimensions[1] = rasDimensions[1];
  rasHDimensions[2] = rasDimensions[2];
  rasHDimensions[3] = 0.0;
  rasHCenter[0] = rasCenter[0];
  rasHCenter[1] = rasCenter[1];
  rasHCenter[2] = rasCenter[2];
  rasHCenter[3] = 1.0;

  //
  // figure out how big that volume is on this particular slice plane
  //
  vtkMatrix4x4 *rasToSlice = vtkMatrix4x4::New();
  rasToSlice->DeepCopy(sliceNode->GetSliceToRAS());
  rasToSlice->MultiplyPoint( rasHCenter, sliceHCenter );
  rasToSlice->SetElement(0, 3, 0.0);
  rasToSlice->SetElement(1, 3, 0.0);
  rasToSlice->SetElement(2, 3, 0.0);
  rasToSlice->Invert();
  rasToSlice->MultiplyPoint( rasHDimensions, sliceHDimensions );
  rasToSlice->Delete();

  // ignore homogeneous coordinate
  sliceDimensions[0] = sliceHDimensions[0];
  sliceDimensions[1] = sliceHDimensions[1];
  sliceDimensions[2] = sliceHDimensions[2];
  sliceCenter[0] = sliceHCenter[0];
  sliceCenter[1] = sliceHCenter[1];
  sliceCenter[2] = sliceHCenter[2];
}

// Get the spacing of the volume, transformed to slice space
double *vtkSlicerSliceLogic::GetBackgroundSliceSpacing()
{
  vtkMatrix4x4 *ijkToRAS = vtkMatrix4x4::New();
  vtkMatrix4x4 *rasToSlice = vtkMatrix4x4::New();
  vtkMatrix4x4 *ijkToSlice = vtkMatrix4x4::New();

  vtkMRMLSliceNode *sliceNode = this->GetSliceNode();
  vtkMRMLSliceCompositeNode *compositeNode = this->GetSliceCompositeNode();
  vtkMRMLVolumeNode *backgroundNode = NULL;

  if ( !sliceNode || !compositeNode )
    {
    return (this->SliceSpacing);
    }
  
  backgroundNode = vtkMRMLVolumeNode::SafeDownCast (
      this->MRMLScene->GetNodeByID( compositeNode->GetBackgroundVolumeID() ));

  if ( !backgroundNode )
    {
    return (this->SliceSpacing);
    }

  backgroundNode->GetIJKToRASMatrix(ijkToRAS);
  rasToSlice->DeepCopy(sliceNode->GetSliceToRAS());
  rasToSlice->Invert();

  ijkToSlice->Multiply4x4(rasToSlice, ijkToRAS, ijkToSlice);

  double invector[4];
  invector[0] = invector[1] = invector[2] = 1.0;
  invector[3] = 0.0;
  double spacing[4];
  ijkToSlice->MultiplyPoint(invector, spacing);
  int i;
  for (i = 0; i < 3; i++)
    {
    this->SliceSpacing[i] = fabs(spacing[i]);
    }

  ijkToRAS->Delete();
  rasToSlice->Delete();
  ijkToSlice->Delete();

  return (this->SliceSpacing);
}

void vtkSlicerSliceLogic::GetBackgroundSliceBounds(double sliceBounds[6])
{
  sliceBounds[0] = sliceBounds[1] = 0.0;
  sliceBounds[2] = sliceBounds[3] = 0.0;
  sliceBounds[4] = sliceBounds[5] = 0.0;

  vtkMRMLSliceNode *sliceNode = this->GetSliceNode();

  if ( !sliceNode )
    {
    return;
    }
  
  double rasDimensions[3], rasCenter[3];
  double rasHMin[4], rasHMax[4]; 
  double sliceHMin[4], sliceHMax[4]; 
  this->GetBackgroundRASDimensions(rasDimensions, rasCenter);
  rasHMin[0] = rasCenter[0] - rasDimensions[0] / 2.;
  rasHMin[1] = rasCenter[1] - rasDimensions[1] / 2.;
  rasHMin[2] = rasCenter[2] - rasDimensions[2] / 2.;
  rasHMin[3] = 1.;
  rasHMax[0] = rasCenter[0] + rasDimensions[0] / 2.;
  rasHMax[1] = rasCenter[1] + rasDimensions[1] / 2.;
  rasHMax[2] = rasCenter[2] + rasDimensions[2] / 2.;
  rasHMax[3] = 1.;

  //
  // figure out how big that volume is on this particular slice plane
  //
  vtkMatrix4x4 *rasToSlice = vtkMatrix4x4::New();
  rasToSlice->DeepCopy(sliceNode->GetSliceToRAS());
  rasToSlice->SetElement(0, 3, 0.0);
  rasToSlice->SetElement(1, 3, 0.0);
  rasToSlice->SetElement(2, 3, 0.0);
  rasToSlice->Invert();
  rasToSlice->MultiplyPoint( rasHMin, sliceHMin );
  rasToSlice->MultiplyPoint( rasHMax, sliceHMax );
  rasToSlice->Delete();

  // ignore homogeneous coordinate
  sliceBounds[0] = sliceHMin[0];
  sliceBounds[1] = sliceHMax[0];
  sliceBounds[2] = sliceHMin[1];
  sliceBounds[3] = sliceHMax[1];
  sliceBounds[4] = sliceHMin[2];
  sliceBounds[5] = sliceHMax[2];
}

// adjust the node's field of view to match the extent of current background volume
void vtkSlicerSliceLogic::FitSliceToBackground(int width, int height)
{
  vtkMRMLSliceNode *sliceNode = this->GetSliceNode();
  vtkMRMLSliceCompositeNode *compositeNode = this->GetSliceCompositeNode();
  vtkMRMLVolumeNode *backgroundNode = NULL;

  if ( !sliceNode || !compositeNode )
    {
    return;
    }
  
  backgroundNode = vtkMRMLVolumeNode::SafeDownCast (
      this->MRMLScene->GetNodeByID( compositeNode->GetBackgroundVolumeID() ));

  vtkImageData *backgroundImage;
  if ( !backgroundNode || ! (backgroundImage = backgroundNode->GetImageData()) )
    {
    return;
    }

  double rasDimensions[3], rasCenter[3];
  this->GetBackgroundRASDimensions (rasDimensions, rasCenter);
  double sliceDimensions[3], sliceCenter[3];
  this->GetBackgroundSliceDimensions (sliceDimensions, sliceCenter);

  double absSliceDimensions[4];
  double fitX, fitY, fitZ, displayX, displayY;
  displayX = fitX = absSliceDimensions[0] = fabs(sliceDimensions[0]);
  displayY = fitY = absSliceDimensions[1] = fabs(sliceDimensions[1]);
  fitZ = absSliceDimensions[2] = fabs(sliceDimensions[2]);


  // fit fov to min dimension of window
  double pixelSize;
  if ( height > width )
    {
    pixelSize = fitX / (1.0 * width);
    fitY = pixelSize * height;
    }
  else
    {
    pixelSize = fitY / (1.0 * height);
    fitX = pixelSize * width;
    }

  // if volume is still too big, shrink some more
  if ( displayX > fitX )
    {
    fitY = fitY / ( fitX / (displayX * 1.0) );
    fitX = displayX;
    }
  if ( displayY > fitY )
    {
    fitX = fitX / ( fitY / (displayY * 1.0) );
    fitY = displayY;
    }

  sliceNode->SetFieldOfView(fitX, fitY, fitZ);

  //
  // set the origin to be the center of the volume in RAS
  //
  vtkMatrix4x4 *sliceToRAS = vtkMatrix4x4::New();
  sliceToRAS->DeepCopy(sliceNode->GetSliceToRAS());
  sliceToRAS->SetElement(0, 3, rasCenter[0]);
  sliceToRAS->SetElement(1, 3, rasCenter[1]);
  sliceToRAS->SetElement(2, 3, rasCenter[2]);
  sliceNode->GetSliceToRAS()->DeepCopy(sliceToRAS);
  sliceToRAS->Delete();

}

// adjust the node's field of view to match the extent of all volume layers
void vtkSlicerSliceLogic::FitSliceToAll(int width, int height)
{
  // TODO...
  vtkErrorMacro ("not yet implemented");
}

// Get/Set the current distance from the origin to the slice plane

double vtkSlicerSliceLogic::GetSliceOffset()
{
  //
  // - get the current translation in RAS space and convert it to Slice space
  //   by transforming it by the inverse of the upper 3x3 of SliceToRAS
  // - pull out the Z translation part
  //

  vtkMRMLSliceNode *sliceNode = this->GetSliceNode();

  if ( !sliceNode )
    {
    return 0.0;
    }

  vtkMatrix4x4 *sliceToRAS = vtkMatrix4x4::New();
  sliceToRAS->DeepCopy( this->SliceNode->GetSliceToRAS() );
  for (int i = 0; i < 3; i++)
    {
    sliceToRAS->SetElement( i, 3, 0.0 );  // Zero out the tranlation portion
    }
  sliceToRAS->Invert();
  double v1[4], v2[4];
  for (int i = 0; i < 4; i++)
    { // get the translation back as a vector
    v1[i] = sliceNode->GetSliceToRAS()->GetElement( i, 3 );
    }
  // bring the translation into slice space
  // and overwrite the z part
  sliceToRAS->MultiplyPoint(v1, v2);
  sliceToRAS->Delete();

  return ( v2[2] );

}

void vtkSlicerSliceLogic::SetSliceOffset(double offset)
{
  //
  // Set the Offset
  // - get the current translation in RAS space and convert it to Slice space
  //   by transforming it by the invers of the upper 3x3 of SliceToRAS
  // - replace the z value of the translation with the new value given by the slider
  // - this preserves whatever translation was already in place
  //

  vtkMRMLSliceNode *sliceNode = this->GetSliceNode();

  if ( !sliceNode )
    {
    return;
    }

  vtkMatrix4x4 *sliceToRAS = vtkMatrix4x4::New();
  sliceToRAS->DeepCopy( this->SliceNode->GetSliceToRAS() );
  for (int i = 0; i < 3; i++)
    {
    sliceToRAS->SetElement( i, 3, 0.0 );  // Zero out the tranlation portion
    }
  sliceToRAS->Invert();
  double v1[4], v2[4];
  for (int i = 0; i < 4; i++)
    { // get the translation back as a vector
    v1[i] = sliceNode->GetSliceToRAS()->GetElement( i, 3 );
    }
  // bring the translation into slice space
  // and overwrite the z part
  sliceToRAS->MultiplyPoint(v1, v2);

  v2[2] = offset;

  // Now bring the new translation vector back into RAS space
  sliceToRAS->Invert();
  sliceToRAS->MultiplyPoint(v2, v1);
  for (int i = 0; i < 4; i++)
    {
    sliceToRAS->SetElement( i, 3, v1[i] );
    }
 
  // if the translation has changed, update the rest of the matrices
  if ( sliceToRAS->GetElement( 0, 3 ) != sliceNode->GetSliceToRAS()->GetElement( 0, 3 ) ||
       sliceToRAS->GetElement( 1, 3 ) != sliceNode->GetSliceToRAS()->GetElement( 1, 3 ) ||
       sliceToRAS->GetElement( 2, 3 ) != sliceNode->GetSliceToRAS()->GetElement( 2, 3 ) )
    {
    sliceNode->GetSliceToRAS()->DeepCopy( sliceToRAS );
    sliceNode->UpdateMatrices();
    }
  sliceToRAS->Delete();
}


