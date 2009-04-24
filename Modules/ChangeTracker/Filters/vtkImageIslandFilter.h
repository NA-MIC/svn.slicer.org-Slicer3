/*=========================================================================

  Program:   Slicer3

  Copyright (c) Brigham and Women's Hospital (BWH) All Rights Reserved.

  See License.txt or http://www.slicer.org/copyright/copyright.txt for details.

==========================================================================*/

// .NAME vtkImageIslandFilter

#ifndef __vtkImageIslandFilter_h
#define __vtkImageIslandFilter_h

#define IMAGEISLANDFILTER_NEIGHBORHOOD_3D 6
#define IMAGEISLANDFILTER_NEIGHBORHOOD_2D 4

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include "assert.h"
 
#include "vtkImageData.h"
#include "vtkImageToImageFilter.h"
#include "vtkChangeTracker.h"
//BTX
template<class T> class IslandMemory {
public:
  // If you provide NewID, which means it is unequal to -1, then make sure that it is unique 
  int AddIsland(int NewStartVoxel, int NewSize, T NewLabel, int NewID, int MaxSize = -1);
  int AddIsland(IslandMemory* NewIslandMem) {return this->AddIsland(NewIslandMem->StartVoxel, NewIslandMem->Size, NewIslandMem->Label, NewIslandMem->ID);}
  
  IslandMemory<T>* DeleteIsland(int DelID);
  IslandMemory* GetIsland(int getID) {
     IslandMemory* Ptr = this; 
     while ((Ptr) && (getID != Ptr->ID)) Ptr = Ptr->Next;  
     return Ptr;
  }
  IslandMemory* GetNext() {return this->Next;}

  int GetSize() {return this->Size;}
  void SetSize(int newSize, IslandMemory* SetIsland, int MaxSize = -1);

  void SetLabel(T newLabel) {this->Label = newLabel;} 
  T GetLabel() {return this->Label ;}
  
  int GetID() {return this->ID;}
  int GetStartVoxel() {return this->StartVoxel;}
  void Print();

  // Prints out paramteres in a line and returns the number of islands in the stack
  int PrintLine ();
  void PrintLine2 ();

  int NumberOfIslands() { return (this->Next ? this->Next->NumberOfIslands() +1 : (this->ID > -1)); }

  ~IslandMemory(){if (this->Next) delete this->Next; this->Next =NULL;}
  IslandMemory(){this->CreateVariables(); }


protected:
  void CreateVariables ();
  int   ID;
  int   StartVoxel;
  int   Size;
  T  Label;
  IslandMemory* Next;
};

// This is so that Islands are grouped by size -> Makes searching a lot faster
template<class T> class IslandMemoryGroup {
public:
  // If you provide NewID, which means it is unequal to -1, then make sure that it is unique 
  int AddIsland(int NewStartVoxel, int NewSize, T NewLabel, int NewID);
  int AddIsland(IslandMemory<T>* NewIslandMem) {return this->AddIsland(NewIslandMem->StartVoxel, NewIslandMem->Size, NewIslandMem->Label, NewIslandMem->ID);}
  
  IslandMemory<T>* DeleteIsland(int DelID, int DelSize);
  IslandMemory<T>* GetIsland(int GetID, int GetSize); 

  IslandMemory<T>* GetNext() {
    if (!this->List) return NULL;
    IslandMemory<T>* result = this->List->GetNext();
    if (result || !this->Next) return result;
    return this->Next->GetNext(); 
  }
  IslandMemoryGroup<T>* GetNextGroup() {return this->Next;} 
  IslandMemory<T>* GetList() { return this->List;} 


  int GetMaxSize() {return this->MaxSize;}
  void SetMaxSize(int initMaxSize) {this->MaxSize = initMaxSize;}

  // Have to do it this way bc of MaxSize
  int GetSize();
  void SetSize(int newSize, IslandMemory<T>* SetIsland);

  void SetLabel(T newLabel) {assert(this->List); this->List->SetLabel(newLabel);} 
  T GetLabel() {assert(this->List);return this->List->GetLabel();}
  
  int GetID() {assert(this->List); return this->List->GetID();}
  int GetStartVoxel() {assert(this->List); return this->List->GetStartVoxel();}

  // Prints out paramteres in a line and returns the number of islands in the stack
  int PrintLine (); 

  int NumberOfIslands() {  
    if (this->Size == -1) return 0; 
    IslandMemoryGroup<T>* ptr = this;
    int result = 0;
    while (ptr) {
      result += ptr->List->NumberOfIslands();
      ptr = ptr->Next;
    }
    return result;
  }

  IslandMemoryGroup* GetGroup(int getSize) {
     IslandMemoryGroup<T>* ptr = this;
     getSize = (this->MaxSize < getSize ? this->MaxSize : getSize);
     while(ptr && ptr->Size != getSize) ptr = ptr->Next;
     return ptr;
  }

  ~IslandMemoryGroup(){
    if (this->List) delete List; List = NULL; 
    if (this->Next) delete this->Next; this->Next =NULL;
  }
  IslandMemoryGroup(){this->CreateVariables(); }

protected:
  IslandMemoryGroup<T>* AddGroup(int NewSize);

  void CreateVariables () {this->Size = -1; this->List = NULL; this->Next = NULL; this->MaxSize = -1;}
  // Minimum IslandSize current stored in stack 
  int   Size;
  int  MaxSize;
  IslandMemory<T>* List;
  IslandMemoryGroup* Next;
};


template<class T> class EMStack {
  public:
  bool Pop(T& Result) { 
     Result = this->Entry; 
     if (this->Next) {
       EMStack *Del = this->Next;  
       this->Entry  = this->Next->Entry; 
       this->Next   = this->Next->Next;
       Del->Next = NULL; 
       delete Del; 
       return true;
     } 
     if (this->Valid) {
        this->Valid = false;
        return true;
     }
     return false;
  }

  void Push(T NewEntry) {
    if (this->Valid) {
      EMStack<T> *Ptr = this->Next;
      this->Next = new EMStack<T>;
      this->Next->Valid = true;
      this->Next->Entry = this->Entry;
      this->Next->Next = Ptr;
      this->Entry = NewEntry;
    } else {
      this->Entry = NewEntry; 
      this->Valid = true;   
    } 
  }
  
  void Print();

  ~EMStack () {if (this->Next) delete this->Next; this->Next = NULL;}
  EMStack() {this->Valid = false; this->Next = NULL; }

  private : 
    T Entry;
    bool Valid;
    EMStack* Next;
};
//ETX

class VTK_CHANGETRACKER_EXPORT vtkImageIslandFilter : public vtkImageToImageFilter
{
  public:
  static vtkImageIslandFilter *New();
  vtkTypeMacro(vtkImageIslandFilter,vtkImageToImageFilter);
  void PrintSelf(ostream& os, vtkIndent indent);

  vtkSetMacro(IslandMinSize,int);
  vtkGetMacro(IslandMinSize,int);
  
  //Description:
  // If set only islands with that label will be deleted - runs a lot faster
  vtkSetMacro(IslandInputLabelMin,int);
  vtkGetMacro(IslandInputLabelMin,int);

  //Description:
  // If set only islands with that label will be deleted - runs a lot faster
  vtkSetMacro(IslandInputLabelMax,int);
  vtkGetMacro(IslandInputLabelMax,int);

  //Description:
  // If set only islands with that label will be deleted - runs a lot faster
  void SetIslandInputLabel(int init) {this->SetIslandInputLabelMin(init);this->SetIslandInputLabelMax(init);}
 
  // Description:
  // If IslandInputLabelMin and ...Max is set then islands below minimum size will be set to this label
  // If this variable is not set then island will be assigned to largest island around 
  vtkSetMacro(IslandOutputLabel,int);
  vtkGetMacro(IslandOutputLabel,int);

  // Description:
  // Set IslandROI has to be of type short and same dimension as input 
  // if IslandROI is set then only those islands will be removed that atleast have one voxel in the IslandROI
  // If the variable is not set then the entire image is considered 
  // ROI is defined by values unequal to 0
  void SetIslandROI(vtkImageData * init) {IslandROI = init;}

  // Description:
  // Print information throughout processing the pipeline
  // 0 = No information is printed out 
  // 1 = Important information is printed out
  // 2 = Comprehensive information
  vtkSetMacro(PrintInformation,int); 
  vtkGetMacro(PrintInformation,int); 


  // Description:
  // If this flag is true then only islands are removed that are holes, i.e. neighbors all have the same label and it is no placed along the edge.
  // Currently only works for dynamic assignment 
  vtkSetMacro(RemoveHoleOnlyFlag,int); 
  vtkGetMacro(RemoveHoleOnlyFlag,int); 

  // Description:
  // Returns the size of the largest island in the label range of  [IslandInputLabelMin,IslandInputLabelMax]
  int GetMaxIslandSize(vtkImageData *InputData);

  // Description:
  // Do you want to do island removal slice by slice or in 3D 
  void SetNeighborhoodDim3D()  {this->NeighborhoodDim = IMAGEISLANDFILTER_NEIGHBORHOOD_3D;} 
  void SetNeighborhoodDim2D()  {this->NeighborhoodDim = IMAGEISLANDFILTER_NEIGHBORHOOD_2D;} 
  vtkGetMacro(NeighborhoodDim,int);


protected:

  vtkImageIslandFilter();
  vtkImageIslandFilter(const vtkImageIslandFilter&) {};
  ~vtkImageIslandFilter();

  void operator=(const vtkImageIslandFilter&) {};
  
  // When it works on parallel machines use : 
  //  void ThreadedExecute(vtkImageData *inData, vtkImageData *outData,int outExt[6], int id);
  // If you do not want to have it multi threaded 
  void ExecuteData(vtkDataObject *);
  void ExecuteInformation(){this->vtkImageToImageFilter::ExecuteInformation();};
  void ExecuteInformation(vtkImageData *inData,vtkImageData *outData);
  void ComputeInputUpdateExtent(int inExt[6], int outExt[6]);

  int IslandMinSize;        // Smalles size of islands allowed, otherwise will be erased
  int IslandInputLabelMin;
  int IslandInputLabelMax;
  int IslandOutputLabel;

  int PrintInformation;

  int NeighborhoodDim;
 
  int RemoveHoleOnlyFlag; 
  // this is not so clean but works 
  vtkImageData *IslandROI;

};
#endif



 







