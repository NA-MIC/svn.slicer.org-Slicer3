/*=========================================================================

Program:   MIMX Meshing Toolkit
Module:    $RCSfile: vtkMimxMeshActor.h,v $
Language:  C++
Date:      $Date: 2008/05/05 19:29:58 $
Version:   $Revision: 1.10 $

 Musculoskeletal Imaging, Modelling and Experimentation (MIMX)
 Center for Computer Aided Design
 The University of Iowa
 Iowa City, IA 52242
 http://www.ccad.uiowa.edu/mimx/
 
Copyright (c) The University of Iowa. All rights reserved.
See MIMXCopyright.txt or http://www.ccad.uiowa.edu/mimx/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

// .NAME vtkMimxMeshActor - 
// .SECTION Description
// vtkMimxMeshActor is the class to handle display and manipulation
// of meshes. Previously this was handled by vtkMimxUnstructuredGridActor
// but this is shared by the Bounding Boxes and additional features
// required prompted the move to a specific class for this type of
// object.

#ifndef __vtkMimxMeshActor_h
#define __vtkMimxMeshActor_h

#include <string>
#include <list>

#include "vtkMimxActorBase.h"


class vtkActor;
class vtkDataSetMapper;
class vtkIdList;
class vtkUnstructuredGrid;
class vtkShrinkFilter;
class vtkPolyDataMapper;
class vtkTubeFilter;
class vtkRenderer;
class vtkGeometryFilter;
class vtkExtractGeometry;
class vtkScalarBarActor;
class vtkLookupTable;
class vtkPlaneWidget;
class vtkPlane;
class vtkRenderWindowInteractor;
class vtkExtractCells;
class vtkPointSet;
class vtkDoubleArray;

class MeshDisplayProperty
{
public:
  std::string name;
  bool IsVisible;
  int DisplayType;
  vtkExtractCells *ExtractCellsFilter;
  vtkActor *SurfaceActor;
  vtkActor *OutlineActor;
  vtkActor *InteriorActor;
  vtkDataSetMapper *SurfaceMapper;
  vtkPolyDataMapper *OutlineMapper;
  vtkDataSetMapper *InteriorMapper;
  vtkShrinkFilter *ShrinkFilter;
  vtkGeometryFilter *GeometryFilter;
  vtkTubeFilter *TubeFilter;
  vtkShrinkFilter *InteriorShrinkFilter;
  std::string activeAttribute;
};

class vtkMimxMeshActor : public vtkMimxActorBase
{
public:
  static vtkMimxMeshActor *New();
  vtkTypeRevisionMacro(vtkMimxMeshActor,vtkMimxActorBase);
  void PrintSelf(ostream& os, vtkIndent indent);
  
  vtkUnstructuredGrid* GetDataSet();
  void SetDataSet( vtkUnstructuredGrid *ugrid);

  enum { 
    DisplaySurface                = 1,
    DisplayOutline                = 2,
    DisplaySurfaceAndOutline      = 3
  };
  
  enum { 
    DisplayMesh                   = 1,
    DisplayElementSets            = 2
  };
  
  vtkSetMacro(ElementSetName, char*);
  vtkGetMacro(ElementSetName, char*);
  
  vtkGetMacro(OutlineActor, vtkActor*);
  vtkGetMacro(InteriorActor, vtkActor*);
  
  vtkGetMacro(NumberOfElementSets, int);
 
  vtkGetMacro(IsAverageEdgeLengthCalculated, int);
  vtkGetMacro(AverageEdgeLength, double);

  void DeleteNodeSet(const char *Name);
  void DeleteElementSet(const char *Name);
  
  void SetRenderer(vtkRenderer *renderer);
  vtkRenderer* GetRenderer();
  
  void SetInteractor(vtkRenderWindowInteractor *interactor);
  vtkRenderWindowInteractor *GetInteractor( );
  
  int  GetDisplayMode( );
  void SetDisplayMode( int mode );
  
  void ShowMesh();
  void HideMesh();
  
  void SetLegendTextColor( double color[3] );
  double *GetLegendTextColor( );
  
  void SetMeshDisplayType( int mode );
  int  GetMeshDisplayType();
  void SetMeshOutlineColor(double red, double green, double blue);
  void SetMeshOutlineColor(double rgb[3]);
  void GetMeshOutlineColor(double &red, double &green, double &blue);
  void GetMeshOutlineColor(double rgb[3]);
  void SetMeshOutlineRadius(double radius);
  double GetMeshOutlineRadius(); 
  void SetMeshShrinkFactor(double shrinkFactor);
  double GetMeshShrinkFactor();
  void SetMeshColor(double red, double green, double blue);
  void SetMeshColor(double rgb[3]);
  void GetMeshColor(double &red, double &green, double &blue);
  void GetMeshColor(double rgb[3]);
  void SetMeshOpacity(double opacity);
  double GetMeshOpacity();
  void SetMeshScalarVisibility(bool visibility);
  bool GetMeshScalarVisibility( );
  void SetMeshLegendVisibility(bool visible);
  bool GetMeshLegendVisibility( );
  void SetLegendRange(double min, double max);
  void SetMeshScalarName(std::string scalarName);
  std::string GetMeshScalarName( );
  void EnableMeshCuttingPlane();
  void DisableMeshCuttingPlane();
  
  void SetInvertCuttingPlane( bool invert );
  bool GetInvertCuttingPlane( );

    
  void SetElementSetDisplayType( std::string name, int mode );
  int  GetElementSetDisplayType( std::string name );
  void SetElementSetOutlineColor(std::string name, double red, double green, double blue);
  void SetElementSetOutlineColor(std::string name, double rgb[3]);
  void GetElementSetOutlineColor(std::string name, double &red, double &green, double &blue);
  void GetElementSetOutlineColor(std::string name, double rgb[3]);
   
  void SetElementSetShrinkFactor(std::string name, double shrinkFactor);
  double GetElementSetShrinkFactor( std::string name );
  void SetElementSetColor(std::string name, double red, double green, double blue);
  void SetElementSetColor(std::string name, double rgb[3]);
  void GetElementSetColor(std::string name, double &red, double &green, double &blue);
  void GetElementSetColor(std::string name, double rgb[3]);
  void SetElementSetOpacity(std::string name, double opacity);
  double GetElementSetOpacity( std::string name );
  void ShowElementSet(std::string name);
  void HideElementSet(std::string name);
  bool GetElementSetVisibility(std::string name);
  void SetElementSetOutlineRadius(std::string name, double radius);
  double GetElementSetOutlineRadius(std::string name); 
  
  void SetElementSetScalarName(std::string setName, std::string scalarName);
  std::string GetElementSetScalarName( std::string setName );
  void EnableElementSetCuttingPlane( std::string setName );
  void DisableElementSetCuttingPlane( std::string setName );
  void SetElementSetScalarVisibility(std::string setName, bool visibility);
  bool GetElementSetScalarVisibility( std::string setName );
  
  void DeleteElementSetListItem( std::string setName );
  void AddElementSetListItem( std::string setName );
  void CreateElementSetList( );
  void DeleteBoundaryConditionStep(int StepNum);
  void ConcatenateStrings(const char* Step, const char* Num, 
          const char* NodeSetName, const char* Type, const char* Direction, char *Name);
  void ChangeElementSetNumbers(const char* ElSetName, int StartEleNum);
  void ChangeNodeSetNumbers(const char* NodeSetName, int StartNodeNum);
  void CalculateAverageEdgeLength();
  vtkPointSet* GetPointSetOfNodeSet(const char* NodeSetName);
  void StoreImageBasedMaterialProperty(const char* ElSetName);
  void StoreImageBasedMaterialPropertyReBin(const char* ElSetName);
  void StoreConstantMaterialProperty(const char* ElSetName, double YoungMod);
  void StoreConstantPoissonsRatio(const char* ElSetName, double PoissonRatio);
  double* ComputeElementSetScalarRange(const char* ElSetName, const char* ArrayName);
protected:
  vtkMimxMeshActor();
  ~vtkMimxMeshActor();
  vtkUnstructuredGrid *UnstructuredGrid;
  vtkDataSetMapper *UnstructuredGridMapper;
  void UpdateMeshDisplay();
  void UpdateElementSetDisplay();
   
private:
  vtkActor *OutlineActor;
  vtkActor *InteriorActor;
  vtkPolyDataMapper *OutlineMapper;
  vtkDataSetMapper *InteriorMapper;
  vtkShrinkFilter *ShrinkFilter;
  vtkShrinkFilter *InteriorShrinkFilter;
  vtkGeometryFilter *OutlineGeometryFilter;
  vtkTubeFilter *TubeFilter;
  vtkRenderer *Renderer;
  vtkRenderWindowInteractor *Interactor;
  vtkPlaneWidget *CuttingPlaneWidget;
  vtkExtractGeometry *ClipPlaneGeometryFilter;
  vtkPlane *CuttingPlane;
  vtkScalarBarActor *LegendActor;
  vtkLookupTable *lutFilter;
  double ElementShrinkFactor;
  bool IsVisible;
        char *ElementSetName;
        int DisplayMode;
        int DisplayType;
        int NumberOfElementSets;
        bool CuttingPlaneEnabled;
        double TextColor[3];
        std::string activeAttribute;
        std::list<MeshDisplayProperty*> ElementSetDisplayList;
        int IsAverageEdgeLengthCalculated;
        double AverageEdgeLength;
        vtkPointSet *PointSetOfNodeSet;
  vtkMimxMeshActor(const vtkMimxMeshActor&);  // Not implemented.
  void operator=(const vtkMimxMeshActor&);  // Not implemented.
};

#endif

