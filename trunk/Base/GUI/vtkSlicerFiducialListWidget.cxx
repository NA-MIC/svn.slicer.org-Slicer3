#include <string>
#include <sstream>

#include "vtkObject.h"
#include "vtkObjectFactory.h"

#include "vtkSlicerFiducialListWidget.h"

#include "vtkSlicerApplicationGUI.h"
#include "vtkSlicerApplication.h"

#include "vtkActor.h"
#include "vtkFollower.h"
#include "vtkProperty.h"
#include "vtkTexture.h"
#include "vtkTransform.h"
#include "vtkPolyData.h"
#include "vtkLookupTable.h"
#include "vtkFloatArray.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkPolyDataMapper.h"
#include "vtkMapper.h"
#include "vtkVectorText.h"
#include "vtkRenderer.h"

#include "vtkMRMLTransformNode.h"
#include "vtkMRMLLinearTransformNode.h"

#include "vtkGlyphSource2D.h"

#include "vtkMRMLFiducialListNode.h"


#include "vtkKWWidget.h"
#include "vtkKWRenderWidget.h"

#include "vtkCollection.h"

#include "vtkSphereSource.h"

#include "vtkTransformPolyDataFilter.h"
#include "vtkGlyph3D.h"

//---------------------------------------------------------------------------
vtkStandardNewMacro ( vtkSlicerFiducialListWidget );
vtkCxxRevisionMacro ( vtkSlicerFiducialListWidget, "$Revision: $");

//---------------------------------------------------------------------------
vtkSlicerFiducialListWidget::vtkSlicerFiducialListWidget ( )
{
  vtkDebugMacro("vtkSlicerFiducialListWidget::Constructor");
  
  this->MainViewer = NULL;
  this->ProcessingMRMLEvent = 0;
  this->RenderPending = 0;
  
  this->DiamondTransformMap.clear();

//  this->Glyph3DList = vtkCollection::New();
  
  // Create the 3d diamond glyphs
  vtkPoints * diamondGlyphPoints = vtkPoints::New();
  diamondGlyphPoints->SetNumberOfPoints(6);
  diamondGlyphPoints->InsertPoint(0, 1, 0, 0);
  diamondGlyphPoints->InsertPoint(1, 0, 1, 0);
  diamondGlyphPoints->InsertPoint(2, 0, 0, 1);
  diamondGlyphPoints->InsertPoint(3, -1, 0, 0);
  diamondGlyphPoints->InsertPoint(4, 0, -1, 0);
  diamondGlyphPoints->InsertPoint(5, 0, 0, -1);

  vtkCellArray * diamondGlyphPolys = vtkCellArray::New();
  diamondGlyphPolys->InsertNextCell( 4 );
  diamondGlyphPolys->InsertCellPoint(0);
  diamondGlyphPolys->InsertCellPoint(1);
  diamondGlyphPolys->InsertCellPoint(3);
  diamondGlyphPolys->InsertCellPoint(4);
  
  diamondGlyphPolys->InsertNextCell(4);
  diamondGlyphPolys->InsertCellPoint(1);
  diamondGlyphPolys->InsertCellPoint(2);
  diamondGlyphPolys->InsertCellPoint(4);
  diamondGlyphPolys->InsertCellPoint(5);

  diamondGlyphPolys->InsertNextCell(4);
  diamondGlyphPolys->InsertCellPoint(2);
  diamondGlyphPolys->InsertCellPoint(0);
  diamondGlyphPolys->InsertCellPoint(5);
  diamondGlyphPolys->InsertCellPoint(3);

  vtkCellArray * diamondGlyphLines = vtkCellArray::New(); 
          
  diamondGlyphLines->InsertNextCell(2);
  diamondGlyphLines->InsertCellPoint(0);
  diamondGlyphLines->InsertCellPoint(3);

  diamondGlyphLines->InsertNextCell(2);
  diamondGlyphLines->InsertCellPoint(1);
  diamondGlyphLines->InsertCellPoint(4);

  diamondGlyphLines->InsertNextCell(2);                                         
  diamondGlyphLines->InsertCellPoint(2);
  diamondGlyphLines->InsertCellPoint(5);

  this->DiamondGlyphPolyData = vtkPolyData::New();
  this->DiamondGlyphPolyData->SetPoints(diamondGlyphPoints);
  //diamondGlyphPoints->Delete();
  this->DiamondGlyphPolyData->SetPolys(diamondGlyphPolys);
  this->DiamondGlyphPolyData->SetLines(diamondGlyphLines);
  //diamondGlyphPolys->Delete();
  //diamondGlyphLines->Delete();

  this->SphereSource = vtkSphereSource::New();
  this->SphereSource->SetRadius(0.3);
  this->SphereSource->SetPhiResolution(10);
  this->SphereSource->SetThetaResolution(10);
/*
  this->DiamondTransform = vtkTransform::New();
  this->DiamondTransform->PostMultiply();
*/
//  this->DebugOn();

}

//---------------------------------------------------------------------------
void vtkSlicerFiducialListWidget::RemoveMRMLObservers()
{
  vtkDebugMacro("vtkSlicerFiducialListWidget::RemoveMRMLObservers\n");
  this->RemoveFiducialObservers();
}

//---------------------------------------------------------------------------
vtkSlicerFiducialListWidget::~vtkSlicerFiducialListWidget ( )
{
  vtkDebugMacro("vtkSlicerFiducialListWidget::Destructor\n");

  // let go of the pointer to the main viewer
  this->MainViewer = NULL;

  this->RemoveMRMLObservers();

  this->DiamondGlyphPolyData->Delete();
  this->SphereSource->Delete();

  vtkDebugMacro("\tDeleting actors...");
  std::map<const char *, vtkActor *>::iterator actorIter;
  for(actorIter = this->DisplayedFiducials.begin();
      actorIter != this->DisplayedFiducials.end();
      actorIter++) 
    {
    if (actorIter->second != NULL)
      {
      actorIter->second->Delete();
      }
    }
  this->DisplayedFiducials.clear();

  std::map<const char *, vtkFollower *>::iterator fIter;
  for(fIter = this->DisplayedTextFiducials.begin();
      fIter != this->DisplayedTextFiducials.end();
      fIter++) 
    {
    if (fIter->second != NULL)
      {
      fIter->second->Delete();
      }
    }
  
  //unsigned int i;

  std::map<const char*, vtkTransform * >::iterator transformIter;
  for (transformIter=this->DiamondTransformMap.begin();
       transformIter != this->DiamondTransformMap.end();
       transformIter++)
    {      
    if (transformIter->second != NULL)
      {
      transformIter->second->Delete();
      }
    }
  this->DiamondTransformMap.clear();

  std::map< const char *, vtkPoints * >::iterator gpIter;
  for (gpIter=  this->GlyphPointsMap.begin();
       gpIter != this->GlyphPointsMap.end();
       gpIter++)
    {
    if (gpIter->second != NULL)
      {
      gpIter->second->Delete();
      }
    }
  this->GlyphPointsMap.clear();
    
  std::map< const char *, vtkFloatArray * >::iterator gsIter;
  for (gsIter = this->GlyphScalarsMap.begin();
       gsIter != this->GlyphScalarsMap.end();
       gsIter++)
    {
    if (gsIter->second != NULL)
      {
      gsIter->second->Delete();
      }
    }
  this->GlyphScalarsMap.clear();

  std::map< const char *, vtkPolyData * >::iterator pdIter;
  for (pdIter = this->GlyphPolyDataMap.begin();
       pdIter != this->GlyphPolyDataMap.end();
       pdIter++)
    {
    if (pdIter->second != NULL)
      {
      pdIter->second->Delete();
      }
    }
  
  for (transformIter=this->TextTransformMap.begin();
       transformIter != this->TextTransformMap.end();
       transformIter++) 
    {
    if (transformIter->second != NULL)
      {
      transformIter->second->Delete();
      }
    }
  this->TextTransformMap.clear();
  
  for (transformIter=this->SymbolTransformMap.begin();
       transformIter != this->SymbolTransformMap.end();
       transformIter++)
    {
    if (transformIter->second != NULL)
      {
      transformIter->second->Delete();
      }
    }
  this->SymbolTransformMap.clear();
  
  std::map< const char *, vtkTransformPolyDataFilter * >::iterator tfIter;
  for (tfIter = this->TransformFilterMap.begin();
       tfIter != this->TransformFilterMap.end();
       tfIter++)
    {
    if (tfIter->second != NULL)
      {
      tfIter->second->Delete();
      }
    }
  this->TransformFilterMap.clear();

  std::map< const char *, vtkMapper * >::iterator gmIter;
  for (gmIter = this->GlyphMapperMap.begin();
       gmIter != this->GlyphMapperMap.end();
       gmIter++)
    {
    if (gmIter->second != NULL)
      {
      gmIter->second->Delete();
      }
    }
  this->GlyphMapperMap.clear();
  
  /*
    if (this->Glyph3DList)
    {
    this->Glyph3DList->RemoveAllItems();
    this->Glyph3DList->Delete();
    }
  */
  std::map< const char *, vtkGlyph3D * >::iterator g3dIter;
  for (g3dIter = this->Glyph3DMap.begin();
       g3dIter != this->Glyph3DMap.end();
       g3dIter++)
    {
    if (g3dIter->second != NULL)
      {
      g3dIter->second->Delete();
      }
    }
  this->Glyph3DMap.clear();
}
//---------------------------------------------------------------------------
void vtkSlicerFiducialListWidget::PrintSelf ( ostream& os, vtkIndent indent )
{
    this->vtkObject::PrintSelf ( os, indent );

    os << indent << "vtkSlicerFiducialListWidget: " << this->GetClassName ( ) << "\n";

    vtkIndent nextIndent;
    nextIndent = indent.GetNextIndent();
    if (this->GetMainViewer() != NULL)
      {
      os << indent << "Main Viewer:\n";
      this->GetMainViewer()->PrintSelf(os, nextIndent);
      }
    
    std::map<const char *, vtkActor *>::iterator iter;
    for(iter=this->DisplayedFiducials.begin(); iter != this->DisplayedFiducials.end(); iter++) 
      {
      os << indent << "Actor " << iter->first << "\n";
      if (iter->second != NULL)
        {
        iter->second->PrintSelf(os, nextIndent);
        }
      else
        {
        os << indent << "vtkActor is null\n";
        }
      }

    std::map<const char *, vtkFollower *>::iterator titer;
    for(titer=this->DisplayedTextFiducials.begin(); titer != this->DisplayedTextFiducials.end(); titer++) 
      {
      os << indent << "Text Actor " << titer->first << "\n";
      if (iter->second != NULL)
        {
        titer->second->PrintSelf(os, nextIndent);
           }
      else
        {
        os << indent << "vtkActor is null\n";
        }
      }


    os << indent << "Maps:\n";
    os << indent << "DiamondTransformMap: size = " << DiamondTransformMap.size() << "\n";
    std::map<const char*, vtkTransform * >::iterator transformIter;
    for (transformIter=this->DiamondTransformMap.begin();
         transformIter != this->DiamondTransformMap.end();
         transformIter++)
      {      
      if (transformIter->second != NULL)
        {
        transformIter->second->PrintSelf(os, nextIndent);
        }
      else
        {
        os << nextIndent << "NULL\n";
        }
      }

    os << indent << "GlyphPointsMap: size = " << this->GlyphPointsMap.size() << "\n";
    std::map< const char *, vtkPoints * >::iterator gpIter;
    for (gpIter=  this->GlyphPointsMap.begin();
         gpIter != this->GlyphPointsMap.end();
         gpIter++)
      {
      if (gpIter->second != NULL)
        {
        gpIter->second->PrintSelf(os, nextIndent);
        }
      else
        {
        os << nextIndent << "NULL\n";
        }
      }

    os << indent << "GlyphScalarsMap: size = " << this->GlyphScalarsMap.size() << "\n";
    std::map< const char *, vtkFloatArray * >::iterator gsIter;
    for (gsIter = this->GlyphScalarsMap.begin();
         gsIter != this->GlyphScalarsMap.end();
         gsIter++)
      {
      if (gsIter->second != NULL)
        {
        gsIter->second->PrintSelf(os, nextIndent);
        }
      else
        {
        os << nextIndent << "NULL\n";
        }
      }
    
    os << indent << "GlyphPolyDataMap: size = " << this->GlyphPolyDataMap.size() << "\n";
    std::map< const char *, vtkPolyData * >::iterator pdIter;
    for (pdIter = this->GlyphPolyDataMap.begin();
         pdIter != this->GlyphPolyDataMap.end();
         pdIter++)
      {
      if (pdIter->second != NULL)
        {
        pdIter->second->PrintSelf(os, nextIndent);
        }
      else
        {
        os << nextIndent << "NULL\n";
        }
      }
    
    os << indent << "TextTransformMap: size = " << this->TextTransformMap.size() << "\n";
    for (transformIter=this->TextTransformMap.begin();
         transformIter != this->TextTransformMap.end();
         transformIter++) 
      {
      if (transformIter->second != NULL)
        {
        transformIter->second->PrintSelf(os, nextIndent);
        }
      else
        {
        os << nextIndent << "NULL\n";
        }
      }
    
    os << indent << "SymbolTransformMap: size = " << this->SymbolTransformMap.size() << "\n";
    for (transformIter=this->SymbolTransformMap.begin();
         transformIter != this->SymbolTransformMap.end();
         transformIter++)
      {
      if (transformIter->second != NULL)
        {
        transformIter->second->PrintSelf(os, nextIndent);
        }
      else
        {
        os << nextIndent << "NULL\n";
        }      
      }
    
    os << indent << "TransformFilterMap: size = " << this->TransformFilterMap.size() << "\n";
    std::map< const char *, vtkTransformPolyDataFilter * >::iterator tfIter;
    for (tfIter = this->TransformFilterMap.begin();
         tfIter != this->TransformFilterMap.end();
         tfIter++)
      {
      if (tfIter->second != NULL)
        {
        tfIter->second->PrintSelf(os, nextIndent);
        }
      else
        {
        os << nextIndent << "NULL\n";
        }
      }
    /*
    os << indent << "Glyph3DMap: size = " << this->Glyph3DMap.size() << "\n";
    for (i = 0; i < this->Glyph3DMap.size(); i++)
      {
      if (this->Glyph3DMap[i] != NULL)
        {
        this->Glyph3DMap[i]->PrintSelf(os, nextIndent);
        }
      else
        {
        os << nextIndent << "NULL\n";
        }
      }
    os << indent << "Glyph3DList: size = " << this->Glyph3DList->GetNumberOfItems() << "\n";

    int i;
    for (i = 0; i < this->Glyph3DList->GetNumberOfItems(); i++)
      {
      ((vtkGlyph3D*)this->Glyph3DList->GetItemAsObject(i))->PrintSelf(os,nextIndent);
      }
    */
    os << indent << "Glyph3DMap: size = " << this->Glyph3DMap.size() << "\n";
    std::map< const char *, vtkGlyph3D * >::iterator g3dIter;
    for (g3dIter = this->Glyph3DMap.begin();
         g3dIter != this->Glyph3DMap.end();
         g3dIter++)
      {
      if (g3dIter->second != NULL)
        {
        g3dIter->second->PrintSelf(os, nextIndent);
        }
      else
        {
        os << nextIndent << "NULL\n";
        }  
      }
    
    os << indent << "GlyphMapperMap: size = " << this->GlyphMapperMap.size() << "\n";
    std::map< const char *, vtkMapper * >::iterator gmIter;
    for (gmIter = this->GlyphMapperMap.begin();
         gmIter != this->GlyphMapperMap.end();
         gmIter++)
      {
      if (gmIter->second != NULL)
        {
        gmIter->second->PrintSelf(os, nextIndent);
        }
      else
        {
        os << nextIndent << "NULL\n";
        }      
      }
    
}

//---------------------------------------------------------------------------
void vtkSlicerFiducialListWidget::ProcessWidgetEvents ( vtkObject *caller,
                                                  unsigned long event, 
                                                  void *callData )
{
  vtkDebugMacro("vtkSlicerFiducialListWidget::ProcessWidgetEvents: event = " << event);
} 

//---------------------------------------------------------------------------
void vtkSlicerFiducialListWidget::ProcessMRMLEvents ( vtkObject *caller,
                                                unsigned long event, 
                                                void *callData )
{
  vtkDebugMacro("vtkSlicerFiducialListWidget::ProcessMRMLEvents: processing = " << this->ProcessingMRMLEvent << ", event = " << event);
  
  if (this->ProcessingMRMLEvent != 0 )
    {
    vtkDebugMacro("Returning because already processing an event, " << this->ProcessingMRMLEvent);
    return;
    }

  this->ProcessingMRMLEvent = event;

  vtkDebugMacro("processing event " << event);
//std::cout <<"ProcessMRML Events: just updating mrml...\n";

  this->UpdateFromMRML();
  this->ProcessingMRMLEvent = 0;
  return;





  if (event == vtkMRMLScene::NodeAddedEvent ||
      event == vtkMRMLScene::NodeRemovedEvent)
    {
    vtkDebugMacro("ProcessMRMLEvents: got a node added or removed event\n");
    }
  else if (event == vtkMRMLFiducialListNode::DisplayModifiedEvent)
    {
    // do a more lightweight update on the fiducial list nodes
    vtkDebugMacro("vtkSlicerFiducialListWidget::ProcessMRMLEvents got a vtkMRMLFiducialListNode::DisplayModifiedEvent, just calling update fids from mrml\n");
    this->UpdateFiducialsFromMRML();
    }
  else if (event == vtkMRMLFiducialListNode::FiducialModifiedEvent)
    {
    vtkDebugMacro("vtkSlicerFiducialListWidget::ProcessMRMLEvents got a FiducialModifiedEvent, removing props and updating from mrml...\n");
    this->RemoveFiducialProps ( );
    this->UpdateFiducialsFromMRML();
    }
  else 
    {
    this->UpdateFromMRML();
    }
  this->ProcessingMRMLEvent = 0;
}

//---------------------------------------------------------------------------
void vtkSlicerFiducialListWidget::CreateWidget ( )
{
  vtkDebugMacro("vtkSlicerFiducialListWidget::CreateWidget...\n");
  
  if (this->IsCreated())
    {
    vtkErrorMacro(<< this->GetClassName() << " already created");
    return;
    }
  
  // Call the superclass to create the whole widget
  this->Superclass::CreateWidget();
}

//---------------------------------------------------------------------------
void vtkSlicerFiducialListWidget::UpdateFromMRML()
{
  vtkDebugMacro("vtkSlicerFiducialListWidget::UpdateFromMRML");
 
  this->RemoveFiducialProps ( );

  this->UpdateFiducialsFromMRML();
}

//---------------------------------------------------------------------------
void vtkSlicerFiducialListWidget::AddList(vtkMRMLFiducialListNode *flist)
{
  if (flist == NULL)
    {
    return;
    }
  
  vtkDebugMacro("\n\n***\nvtkSlicerFiducialListWidget::AddList: starting...\n");

  const char* fid = flist->GetID();
//  std::cout << "\n\n***\nvtkSlicerFiducialListWidget::AddList: fid = " << fid << endl;
  
  float textScale = flist->GetTextScale();
  float symbolScale = flist->GetSymbolScale();
  int glyphType = flist->GetGlyphType();

  // are all the current vectors the same size?

//  std::cout << "vtkSlicerFiducialListWidget::AddList: have glyph type of " << glyphType << endl;
  //this->GlyphSymbolMap.push_back(glyphType);
    
  // where in RAS space will the fiducials be displayed?
  vtkPoints * glyphPoints = vtkPoints::New();
  glyphPoints->Initialize();
  this->GlyphPointsMap[fid] = glyphPoints;
  // let the vector remember it
  //glyphPoints->Delete();

  vtkDebugMacro("...added the new glyph points..., vector size = " << this->GlyphPointsMap.size());
  
  // the scalar array is used to determine the colour of the fiducial,
  // selected or not
  vtkFloatArray * glyphScalars = vtkFloatArray::New();
  this->GlyphScalarsMap[fid] = glyphScalars;
  //glyphScalars->Delete();

  vtkDebugMacro("...added the new scalars...");
  
  // points and scalars are encapsulated in the poly data
  vtkPolyData * glyphPolyData = vtkPolyData::New();
  this->GlyphPolyDataMap[fid] = glyphPolyData;
//  glyphPolyData->Delete();

  vtkDebugMacro("...added the new glyph poly data...");
  
  // short cut for vector access of the new list's structures
  //int n = this->GlyphPolyDataMap.size() - 1;

  vtkDebugMacro("...using n = " << fid << ", poly data size = " << this->GlyphPolyDataMap.size() << ", glyph points vector size = " <<  this->GlyphPointsMap.size() << ", glyph scalars vector size = " << this->GlyphScalarsMap.size() );

  if (this->GlyphPolyDataMap[fid] == NULL) std::cerr << "poly data is null\n"; else std::cerr << "poly data is not null\n";
  if (this->GlyphPointsMap[fid] == NULL) std::cerr << "points is null\n"; else std::cerr << "points is not null\n";

  this->GlyphPolyDataMap[fid]->SetPoints(this->GlyphPointsMap[fid]);
  this->GlyphPolyDataMap[fid]->GetPointData()->SetScalars(this->GlyphScalarsMap[fid]);

  vtkDebugMacro("...set the points, and the scalars...");
  
  // the default size for the text
  vtkTransform *textTransform = vtkTransform::New();
  textTransform->AddObserver(vtkCommand::WarningEvent, this->MRMLCallbackCommand );
  int textPush = 10;
  textTransform->Translate(0, 0, textPush);
  textTransform->GetMatrix()->SetElement(0, 1, 0.333);
  textTransform->Scale(textScale, textScale, 1);
  this->TextTransformMap[fid] = textTransform;
//  textTransform->Delete();
  
  vtkDebugMacro("...added the new text transform...");

  // default size for symbols
  vtkTransform *symbolTransform = vtkTransform::New();
  symbolTransform->AddObserver(vtkCommand::WarningEvent, this->MRMLCallbackCommand );
  symbolTransform->Scale(symbolScale, symbolScale, symbolScale);
  this->SymbolTransformMap[fid] = symbolTransform;
//  symbolTransform->Delete();

  vtkDebugMacro("...added the new symbol transform..");
  
  // set up the shape of the glyph
  vtkTransformPolyDataFilter * transformFilter = vtkTransformPolyDataFilter::New();
  // use the shape built in the constructor
  if (glyphType == vtkMRMLFiducialListNode::Diamond3D)
    {
    transformFilter->SetInput(this->DiamondGlyphPolyData);
    }
  else
    {
//    std::cout << "Using the sphere source\n";
    transformFilter->SetInput(this->SphereSource->GetOutput());
    }
  transformFilter->SetTransform(this->SymbolTransformMap[fid]);
  this->TransformFilterMap[fid] = transformFilter;
//  transformFilter->Delete();

  vtkDebugMacro("...added the new transofrm filter...");
  
  // now set up the glyph
  vtkGlyph3D *glyph3D = vtkGlyph3D::New();
  glyph3D->SetSource(this->TransformFilterMap[fid]->GetOutput());
  glyph3D->SetInput(this->GlyphPolyDataMap[fid]);
  glyph3D->SetScaleFactor(1.0);
  glyph3D->ClampingOn();
  glyph3D->ScalingOff();
  glyph3D->SetRange(0, 1);
  this->Glyph3DMap[fid] = glyph3D;
  //this->Glyph3DList->vtkCollection::AddItem(glyph3D);
  //glyph3D->Delete();

  vtkDebugMacro("...added the new glyph...");
  
  // now set up the mapper
  vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
  mapper->SetInput(this->Glyph3DMap[fid]->GetOutput());
  //mapper->SetInput(((vtkGlyph3D*)this->Glyph3DList->GetItemAsObject(n))->GetOutput());
  vtkLookupTable *lut = vtkLookupTable::SafeDownCast(mapper->GetLookupTable());
  lut->SetNumberOfTableValues(2);
  // set the selected/unselected colours
  lut->SetTableValue(0, 1, 0, 0, 1.0);
  lut->SetTableValue(1, 0, 0, 1, 1.0);
  this->GlyphMapperMap[fid] = mapper;
//  lut->Delete();
//  mapper->Delete();

  vtkDebugMacro("...added the new mapper...");
  
  // set up the list's transform?
  
  /*
  
    this->DiamondGlyphPolyDataMap->SetPoints(glyphPoints);
    vtkFloatArray * glyphScalars = vtkFloatArray::New();
    this->DiamondGlyphPolyDataMap->GetPointData()->SetScalars(glyphScalars);
    GlyphPoints[listNumber]->SetNumberOfPoints(0);
    GlyphScalars[listNumber]->SetNumberOfTuples(0);
  */
  vtkDebugMacro("AddPoints Done...");
}

//---------------------------------------------------------------------------
void vtkSlicerFiducialListWidget::RemoveList(vtkMRMLFiducialListNode * flist)
{
  std::cout << "WARNING: RemoveList not fully implemented yet!\n";
  if (flist != NULL)
    {
    this->RemoveList(flist->GetID());
    }
}

//---------------------------------------------------------------------------
void vtkSlicerFiducialListWidget::RemoveList(const char *id)
{
  if (id == NULL)
    {
    return;
    }
  /*
  std::map<const char *, vtkActor *>::iterator iter;
  iter = this->DisplayedFiducials.find(id);
  if (iter != this->DisplayedFiducials.end())
    {
    this->DisplayedFiducials.erase(iter);
    }
  */
//  std::cout << "Deleting disp fid at " << id << endl;
  this->DisplayedFiducials[id]->Delete();
//  std::cout << "erasing from map...\n";
  this->DisplayedFiducials.erase(id);
  
  this->DisplayedTextFiducials[id]->Delete();
  this->DisplayedTextFiducials.erase(id);

  this->DiamondTransformMap[id]->Delete();
  this->DiamondTransformMap.erase(id);

  this->GlyphPointsMap[id]->Delete();
  this->GlyphPointsMap.erase(id);

  this->GlyphScalarsMap[id]->Delete();
  this->GlyphScalarsMap.erase(id);

  
}

//---------------------------------------------------------------------------
void vtkSlicerFiducialListWidget::RequestRender()
{
    if (this->GetRenderPending())
    {
    return;
    }

  this->SetRenderPending(1);
  this->Script("after idle \"%s Render\"", this->GetTclName());
}

//---------------------------------------------------------------------------
void vtkSlicerFiducialListWidget::Render()
{
  this->MainViewer->Render();
  this->SetRenderPending(0);
}

//---------------------------------------------------------------------------
void vtkSlicerFiducialListWidget::UpdateFiducialsFromMRML()
{
  
  vtkMRMLScene *scene = this->GetMRMLScene();
  if (scene == NULL)
    {
    vtkDebugMacro("... but the scene is null...");
    }
  vtkMRMLNode *node = NULL;

  vtkDebugMacro("UpdateFiducialsFromMRML: Starting to update the viewer's actors, glyphs for the fid lists.");
  
  scene->InitTraversal();
//  std::cout << "\n\nvtkSlicerFiducialListWidget::UpdateFiducialsFromMRML: Start\n";
  //int listNumber = 0;
  while (node=scene->GetNextNodeByClass("vtkMRMLFiducialListNode"))
    {
    vtkMRMLFiducialListNode *flist = vtkMRMLFiducialListNode::SafeDownCast(node);
//    std::cout << "vtkSlicerFiducialListWidget::UpdateFiducialsFromMRML: got a list " << flist->GetID() << ", have " << this->GlyphPointsMap.size() << " glyph points in the vector\n";
    if (flist->HasObserver ( vtkCommand::ModifiedEvent, this->MRMLCallbackCommand ) == 0)
      {
      flist->AddObserver ( vtkCommand::ModifiedEvent, this->MRMLCallbackCommand );
      }
    if (flist->HasObserver ( vtkMRMLTransformableNode::TransformModifiedEvent, this->MRMLCallbackCommand ) == 0)
      {
      flist->AddObserver ( vtkMRMLTransformableNode::TransformModifiedEvent, this->MRMLCallbackCommand );
      }
    // observe display node
    if (flist->HasObserver ( vtkMRMLFiducialListNode::DisplayModifiedEvent, this->MRMLCallbackCommand ) == 0)
      {     
      flist->AddObserver ( vtkMRMLFiducialListNode::DisplayModifiedEvent, this->MRMLCallbackCommand );
      }
    // fiducial point modified?
    if (flist->HasObserver ( vtkMRMLFiducialListNode::FiducialModifiedEvent, this->MRMLCallbackCommand ) == 0)
      {
      flist->AddObserver( vtkMRMLFiducialListNode::FiducialModifiedEvent, this->MRMLCallbackCommand );
      }

    // set up the points at which the glyphs will be shown
    double* selectedColor =  flist->GetSelectedColor();
    double* unselectedColor = flist->GetColor();
    bool use3DSymbols;

    if (0)
      {
      std::cout << "\tList Unselected color: r = " << unselectedColor[0]
                << ", g = " << unselectedColor[1]
                << ", b = " << unselectedColor[2] << endl;
      std::cout << "\tList Selected color:   r = " << selectedColor[0]
                << ", g = " << selectedColor[1]
                << ", b = " << selectedColor[2] << endl;
      std::cout << "\tList Glyph type = " << flist->GetGlyphTypeAsString() << endl;
      }
    
    if (flist->GlyphTypeIs3D())
      {
//      std::cout << "\tusing 3d glyphs\n";
      use3DSymbols = true;
      }
    else
      {
//      std::cout << "\tusing 2d glyphs\n";
      use3DSymbols = false;
      }

    // do we already have the structurse for this list?
    std::map< const char *, vtkPoints * >::iterator gpIter;
    gpIter = this->GlyphPointsMap.find(flist->GetID());
    if (gpIter == this->GlyphPointsMap.end())
      {
      // this id isn't used as a key yet, so add the data structures need for
      // a new list
//      std::cout << "\tAding a new list, list = " << flist->GetID() << endl;
      AddList(flist);      
      vtkDebugMacro("Added the list");
      }
    else
      {
//      std::cout << "\tClearning out the points and scalars for this list " << flist->GetID() << endl;
      // clear out the points list
      this->GlyphPointsMap[flist->GetID()]->SetNumberOfPoints(0);
      this->GlyphScalarsMap[flist->GetID()]->SetNumberOfTuples(0);
      // is this list 2d or 3d?
      if (use3DSymbols)
        {
//        std::cout << "\tDidn't add a new list, but do have 3d symbols, glyph points map now has " <<
//          this->GlyphPointsMap[flist->GetID()]->GetNumberOfPoints() << " points\n";
        }
      }

    if (use3DSymbols)
      {
      
      // make sure that we've got the right glyph for the 3d case
      if (flist->GetGlyphType() == vtkMRMLFiducialListNode::Diamond3D)
        {
        if ( this->TransformFilterMap[flist->GetID()]->GetInput() != this->DiamondGlyphPolyData)
          {
          this->TransformFilterMap[flist->GetID()]->SetInput(this->DiamondGlyphPolyData);
          }
        }
      else
        {
        if (this->TransformFilterMap[flist->GetID()]->GetInput() != this->SphereSource->GetOutput())
          {
          this->TransformFilterMap[flist->GetID()]->SetInput(this->SphereSource->GetOutput());
          }
        }
    
      // set up the selected/unselected colours
//      std::cout << "Setting the 3D glyph mapper's lookup table...table value 0 for unselected, 1 for selected\n";
      if (unselectedColor != NULL)
        {
        vtkLookupTable::SafeDownCast(this->GlyphMapperMap[flist->GetID()]->GetLookupTable())->SetTableValue(0,
                                                                                                            unselectedColor[0],
                                                                                                            unselectedColor[1],
                                                                                                            unselectedColor[2],
                                                                                                            1.0);        
        }
      if (selectedColor != NULL)
        {
        vtkLookupTable::SafeDownCast(this->GlyphMapperMap[flist->GetID()]->GetLookupTable())->SetTableValue(1,
                                                                                                            selectedColor[0],
                                                                                                            selectedColor[1],
                                                                                                            selectedColor[2],
                                                                                                            1.0);
        }
      } // end of 3d symbols

    
    this->GlyphScalarsMap[flist->GetID()]->SetNumberOfTuples(flist->GetNumberOfFiducials());
//    std::cout << "\tAbout to loop through the fid list, for " << flist->GetNumberOfFiducials() << ", set the glyph scalars map to have " << this->GlyphScalarsMap[flist->GetID()]->GetNumberOfTuples() << " tuples\n";
    
    // set up each fiducial on this list 
    for (int f=0; f<flist->GetNumberOfFiducials(); f++)
      {

      // check to see if this fiducial has actors in the DisplayFiducials
      // map
      int actorExists = 0;

      vtkActor * actor = GetFiducialActorByID(flist->GetNthFiducialID(f));
      if (actor != NULL)
        {
        actorExists = 1;
        }
      else
        {
        if (0)
          {
          std::cout << "\tNo actor exists (checked for id " << flist->GetNthFiducialID(f) << "), the displayed fiducials list has " << this->DisplayedFiducials.size() << " and these keys:\n";
          std::map<const char *, vtkActor *>::iterator iter;
          for(iter=this->DisplayedFiducials.begin(); iter != this->DisplayedFiducials.end(); iter++) 
            {
            std::cout << "\t\t" << iter->first << endl;
            }
          }
        }
      //std::cout << "fid = " << f << ": Actor for this point  exists = " << actorExists << " (checked for id " << flist->GetNthFiducialID(f) << ")" << endl;
      
      vtkPolyDataMapper *mapper = NULL;
//      vtkActor *actor = NULL;
      vtkGlyph3D * glyph3d = NULL;
      vtkGlyphSource2D *glyph2d = vtkGlyphSource2D::New();

      // if this list is using 3d symbols
      if (use3DSymbols)
        {

        // add this point to the list of points
        this->GlyphPointsMap[flist->GetID()]->InsertNextPoint(flist->GetNthFiducialXYZ(f));
//        std::cout << "\tadded the next point to the glyph points map, " << f << " = " << flist->GetNthFiducialXYZ(f) << ", glyph points map now has " << this->GlyphPointsMap[flist->GetID()]->GetNumberOfPoints() << " points." << endl;
        
//        glyph3d = this->Glyph3DMap[flist->GetID()];
        glyph3d = (vtkGlyph3D*)this->Glyph3DMap[flist->GetID()]; // ->GetItemAsObject(listNumber);
        
        if (glyph3d != NULL)
          {
          if (actorExists)
            {
            // actor is in the list, get it and the mapper
            //actor = iter->second;
            glyph3d->SetOutput(actor->GetMapper()->GetInput());
//            std::cout << "\tthe glyph actor exists, set the glyphoutput to the actors mapper input\n";
            }
          else
            {
            // no actor, allocate vars and set up the pipeline
            actor = vtkActor::New ( );
//            std::cout << "\tthe glyph actor didn't exist, setting it's mapper to the glyph mapper vector element at " << flist->GetID() << endl;
            actor->SetMapper ( this->GlyphMapperMap[flist->GetID()] );
            this->MainViewer->AddViewProp ( actor );
            }
          
//          std::cout << "\tSetting scalars vector tuples at fid # " << f << " ...number of tuples = " << this->GlyphScalarsMap[flist->GetID()]->GetNumberOfTuples() << "\n";

          
          // reset the glyph actor's colours
          vtkLookupTable::SafeDownCast(actor->GetMapper()->GetLookupTable())->SetTableValue(0, unselectedColor[0],
                                                                                            unselectedColor[1],
                                                                                            unselectedColor[2],
                                                                                            1.0);
          vtkLookupTable::SafeDownCast(actor->GetMapper()->GetLookupTable())->SetTableValue(1,
                                                                                            selectedColor[0],
                                                                                            selectedColor[1],
                                                                                            selectedColor[2],
                                                                                            1.0);
          if (0)
          {
            std::cout << "\tactor's mapper lut values:\n";
            double rgba[4];
            vtkLookupTable::SafeDownCast(actor->GetMapper()->GetLookupTable())->GetTableValue(0, rgba);
            std::cout << "\t0 unselected: r = " << rgba[0] << ", g = " << rgba[1] << ", b = " << rgba[2] << endl;
          vtkLookupTable::SafeDownCast(actor->GetMapper()->GetLookupTable())->GetTableValue(1, rgba);
            std::cout << "\t1 selected:  r = " << rgba[0] << ", g = " << rgba[1] << ", b = " << rgba[2] << endl;
          }
          if (flist->GetNthFiducialSelected(f))
            {
          //  std::cout << "\tfid " << f << " is Selected, setting scalar tuple " << f << " to 1\n";
            this->GlyphScalarsMap[flist->GetID()]->SetTuple1(f, 1.0);
            }
          else
            {
//            std::cout << "\tfid " << f << " is unselected, setting tuple " << f << " to 0\n";
            this->GlyphScalarsMap[flist->GetID()]->SetTuple1(f, 0.0);
            }
//          std::cout << "\tafter setting the tuple it's = " << this->GlyphScalarsMap[flist->GetID()]->GetTuple1(f) << endl;
          }
        }

      //if (this->GetSymbolDimension() == this->Use2DSymbols)
      if (!use3DSymbols)
        {
        if (actorExists)
          {
          // actor is in the list, get it and the mapper
          //actor = iter->second;
          glyph2d->SetOutput(actor->GetMapper()->GetInput());              
          }
        else
          {
          // no actor, allocate vars and set up the pipeline
          mapper = vtkPolyDataMapper::New ();
          mapper->SetInput ( glyph2d->GetOutput() );
          actor = vtkActor::New ( );
          actor->SetMapper ( mapper );
          this->MainViewer->AddViewProp ( actor );
          }
        
        if (glyph2d != NULL)
          {
//          glyph2d->SetGlyphTypeToDiamond();
          glyph2d->SetGlyphType(flist->GetGlyphTypeAsVTKEnum());
          if (flist->GetNthFiducialSelected(f))
            {
            glyph2d->SetColor(flist->GetSelectedColor());
            }
          else
            {
            glyph2d->SetColor(flist->GetColor());
            }
          }         
        } // end of 2d
      
      // handle text
      // check to see if this fiducial follower has actors in the
      // DisplayedTextFiducials map
      int textActorExists = 0;
      std::map<const char *, vtkFollower *>::iterator titer;

      titer = this->DisplayedTextFiducials.find(flist->GetNthFiducialID(f));
      if (titer != this->DisplayedTextFiducials.end())
        {
        textActorExists = 1;
        }
//      std::cout << "\tText actor exists = " << textActorExists << endl;
      
      vtkVectorText *vtext = vtkVectorText::New();
      vtkPolyDataMapper *textMapper;
      vtkFollower *textActor;
      if (textActorExists)
        {
        // get it out of the map
        textActor = titer->second;
        vtext->SetOutput(titer->second->GetMapper()->GetInput());
        }
      else
        {
        textMapper = vtkPolyDataMapper::New ();
        textMapper->SetInput ( vtext->GetOutput() );
        
        textActor = vtkFollower::New();
        textActor->SetCamera(this->MainViewer->GetRenderer()->GetActiveCamera());
        textActor->SetMapper(textMapper);
        
        this->MainViewer->AddViewProp ( textActor );
        }
      vtext->SetText(flist->GetNthFiducialLabelText(f));

      // set the display properties on the actor and the text actor
      this->SetFiducialDisplayProperty(flist, f, actor, textActor);

      // save the actors and clean up, if necessary
      if (!actorExists)
        {
        vtkDebugMacro("After setting disp props, actor vis = " << actor->GetVisibility() << ", adding actor to displayed fid list with first string " << flist->GetNthFiducialID(f) << ", current size = " << this->DisplayedFiducials.size());

        
        // need to use a constant string as the key to the map, use the
        // fiducial point's id
//        std::cout << "\tadding an actor to the displayed fiducials list at id " << flist->GetNthFiducialID(f) << ", the list is now keyed with these strings" << endl;
        this->DisplayedFiducials[flist->GetNthFiducialID(f)] = actor;
        std::map<const char *, vtkActor *>::iterator iter;
        for(iter=this->DisplayedFiducials.begin(); iter != this->DisplayedFiducials.end(); iter++) 
          {
//          std::cout << "\t\t" << iter->first << endl;
          }
//        std::cout << "\tnew size of displayed fids = " << this->DisplayedFiducials.size() << endl;

        //if (this->GetSymbolDimension() == this->Use2DSymbols)
       
        // only call delete if made them new, they didn't exist before
        if (glyph2d != NULL)
          {
          glyph2d->Delete();
          }
        if (!use3DSymbols)
          {
          if (actor != NULL)
            {
            actor->Delete();
            }
          if (mapper != NULL)
            {
            mapper->Delete();
            }
          } // end of 2d
        }

      if (!textActorExists)
        {
        this->DisplayedTextFiducials[flist->GetNthFiducialID(f)] = textActor;
        if (1)
          //if (!use3DSymbols)
          {
          // only delete them if made them new
          if (vtext != NULL)
            {
            vtext->Delete();
            }
          if (textMapper != NULL)
            {
            textMapper->Delete();
            }
          if (0) {
          if (textActor != NULL)
            {
            textActor->Delete();
            }
          }
          } // end of 2d
        }
      } // end for

    // go to the next fiducial list 
    //listNumber++;

    // reset the poly data
    if (0)
      {
      std::cout << "Done looping over fids, for list " << flist->GetID() << ", now resetting points and scalars on the poly data map\n";
      for (int s = 0; s < this->GlyphScalarsMap[flist->GetID()]->GetNumberOfTuples(); s++)
        {
        std::cout << "scalar " << s << " = " << this->GlyphScalarsMap[flist->GetID()]->GetTuple1(s) << endl;
        }
    }
    this->GlyphPolyDataMap[flist->GetID()]->SetPoints(this->GlyphPointsMap[flist->GetID()]);
    this->GlyphPolyDataMap[flist->GetID()]->GetPointData()->SetScalars(this->GlyphScalarsMap[flist->GetID()]);

    } // end while
  
  if (this->MainViewer != NULL)
    {
    // render
    this->RequestRender();
    }
}

//---------------------------------------------------------------------------
void vtkSlicerFiducialListWidget::RemoveFiducialProps()
{
  vtkDebugMacro("vtkSlicerFiducialListWidget::RemoveFiducialProps\n");
  // glyph actors

  int idx = 0;
  std::map<const char *, vtkActor *>::iterator iter;
  for(iter=this->DisplayedFiducials.begin(); iter != this->DisplayedFiducials.end(); iter++) 
    {
    this->MainViewer->RemoveViewProp(iter->second);
    iter->second->Delete();
    }
  this->DisplayedFiducials.clear();

  // text actors
  std::map<const char *, vtkFollower *>::iterator titer;
  for(titer=this->DisplayedTextFiducials.begin(); titer != this->DisplayedTextFiducials.end(); titer++) 
    {
    this->MainViewer->RemoveViewProp(titer->second);
    titer->second->Delete();
    }
  this->DisplayedTextFiducials.clear();

}

//---------------------------------------------------------------------------
void vtkSlicerFiducialListWidget::RemoveFiducialObservers()
{
  vtkDebugMacro("vtkSlicerFiducialListWidget::RemoveFiducialObservers\n");
    if (this->GetMRMLScene() == NULL)
    {
    vtkDebugMacro("vtkSlicerFiducialListWidget::RemoveFiducialObservers: no scene, returning...");
        return;
    }
    // remove the observers on all the fiducial lists
    vtkMRMLFiducialListNode *flist;
    this->GetMRMLScene()->InitTraversal();
    while ((flist = vtkMRMLFiducialListNode::SafeDownCast(this->GetMRMLScene()->GetNextNodeByClass("vtkMRMLFiducialListNode"))) != NULL)
    {
        vtkDebugMacro("Removing observers on fiducial list " << flist->GetID());
        if (flist->HasObserver (vtkCommand::ModifiedEvent, this->MRMLCallbackCommand ) == 1)
          {
          flist->RemoveObservers ( vtkCommand::ModifiedEvent, this->MRMLCallbackCommand );
          }
        if (flist->HasObserver( vtkMRMLFiducialListNode::DisplayModifiedEvent, this->MRMLCallbackCommand ) == 1)
          {
          flist->RemoveObservers ( vtkMRMLFiducialListNode::DisplayModifiedEvent, this->MRMLCallbackCommand );
          }
        if (flist->HasObserver( vtkMRMLTransformableNode::TransformModifiedEvent, this->MRMLCallbackCommand ) == 1)
          {
          flist->RemoveObservers ( vtkMRMLTransformableNode::TransformModifiedEvent, this->MRMLCallbackCommand );
          }
        if (flist->HasObserver(vtkMRMLFiducialListNode::FiducialModifiedEvent, this->MRMLCallbackCommand ) == 1)
          {
          flist->RemoveObservers ( vtkMRMLFiducialListNode::FiducialModifiedEvent, this->MRMLCallbackCommand );
          }
    }
}

//---------------------------------------------------------------------------
void vtkSlicerFiducialListWidget::SetFiducialDisplayProperty(vtkMRMLFiducialListNode *flist, 
                                                       int n,
                                                       vtkActor *actor, vtkFollower *textActor)
{
  vtkDebugMacro("vtkSlicerFiducialListWidget::SetFiducialDisplayProperty: n = " << n);
  float *xyz = flist->GetNthFiducialXYZ(n);
  int selected = flist->GetNthFiducialSelected(n);

  actor->SetPosition(xyz[0], xyz[1], xyz[2]);
  actor->SetScale(flist->GetSymbolScale());

  textActor->SetPosition(xyz[0], xyz[1], xyz[2]);
  textActor->SetScale(flist->GetTextScale());

  vtkMRMLTransformNode* tnode = flist->GetParentTransformNode();
  if (tnode != NULL && tnode->IsLinear())
    {
    vtkMRMLLinearTransformNode *lnode = vtkMRMLLinearTransformNode::SafeDownCast(tnode);
    vtkMatrix4x4* transformToWorld = vtkMatrix4x4::New();
    transformToWorld->Identity();
    lnode->GetMatrixTransformToWorld(transformToWorld);
    actor->SetUserMatrix(transformToWorld);
    textActor->SetUserMatrix(transformToWorld);
    transformToWorld->Delete();
    }

  actor->SetVisibility(flist->GetVisibility());
  textActor->SetVisibility(flist->GetVisibility());
  if (selected)
    {
    actor->GetProperty()->SetColor(flist->GetSelectedColor());
    textActor->GetProperty()->SetColor(flist->GetSelectedColor());
    }
  else
    {
    actor->GetProperty()->SetColor(flist->GetColor());
    textActor->GetProperty()->SetColor(flist->GetColor());
    }
  actor->GetProperty()->SetOpacity(flist->GetOpacity());
  textActor->GetProperty()->SetOpacity(flist->GetOpacity());
  actor->GetProperty()->SetAmbient(flist->GetAmbient());
  textActor->GetProperty()->SetAmbient(flist->GetAmbient());
  actor->GetProperty()->SetDiffuse(flist->GetDiffuse());
  textActor->GetProperty()->SetDiffuse(flist->GetDiffuse());
  actor->GetProperty()->SetSpecular(flist->GetSpecular());
  textActor->GetProperty()->SetSpecular(flist->GetSpecular());
  actor->GetProperty()->SetSpecularPower(flist->GetPower());
  textActor->GetProperty()->SetSpecularPower(flist->GetPower());
  actor->SetTexture(NULL);

  vtkDebugMacro("vtkSlicerFiducialListWidget::SetFiducialDisplayProperty: done setting for " << n << "\n");

}

//---------------------------------------------------------------------------
std::string
vtkSlicerFiducialListWidget::GetFiducialNodeID (const char *actorid, int &index)
{
  // take the index off the actor id to get the fiducial node's id
  std::string actorString = actorid;
  std::stringstream ss;
  std::string sid;
  ss << actorid;
  ss >> sid;
  ss >> index;
  return sid;
}

//---------------------------------------------------------------------------
vtkActor *
vtkSlicerFiducialListWidget::GetFiducialActorByID (const char *id)
{
  vtkDebugMacro("vtkSlicerFiducialListWidget::GetFiducialActorByID: id = " << id);
  if ( !id )
    {
    return (NULL);
    }
  std::string sid = id;

  std::map<const char *, vtkActor *>::iterator iter;
  // search for matching string (can't use find, since it would look for 
  // matching pointer not matching content)
  for(iter=this->DisplayedFiducials.begin(); iter != this->DisplayedFiducials.end(); iter++) 
    {
//std::cout << "GetFiducialActorBy ID " << id << " (as string = " << sid.c_str() <<  "), now checking " << iter->first << endl;
    if ( iter->first && !strcmp( iter->first, sid.c_str() ) )
      {
      return (iter->second);
      }
    }
  return (NULL);
}
