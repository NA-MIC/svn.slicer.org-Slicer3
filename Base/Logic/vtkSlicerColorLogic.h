/*=auto=========================================================================

  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See Doc/copyright/copyright.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: vtkSlicerColorLogic.h,v $
  Date:      $Date$
  Version:   $Revision$

=========================================================================auto=*/

///  vtkSlicerColorLogic - slicer logic class for color manipulation
/// 
/// This class manages the logic associated with reading, saving,
/// and changing propertied of the colors


#ifndef __vtkSlicerColorLogic_h
#define __vtkSlicerColorLogic_h

#include <stdlib.h>

#include "vtkSlicerBaseLogic.h"
#include "vtkSlicerLogic.h"

#include "vtkMRML.h"

class vtkMRMLColorNode;
class VTK_SLICER_BASE_LOGIC_EXPORT vtkSlicerColorLogic : public vtkSlicerLogic 
{
  public:
  
  /// The Usual vtk class functions
  static vtkSlicerColorLogic *New();
  vtkTypeRevisionMacro(vtkSlicerColorLogic,vtkSlicerLogic);
  void PrintSelf(ostream& os, vtkIndent indent);

  /// 
  /// Update logic state when MRML scene chenges
  void ProcessMRMLEvents(vtkObject * caller,
                         unsigned long event,
                         void * callData);

  /// 
  /// Add a series of color nodes, setting the types to the defaults, so that
  /// they're accessible to the rest of Slicer
  void AddDefaultColorNodes();

  /// 
  /// Remove the colour nodes that were added
  void RemoveDefaultColorNodes();

  /// 
  /// Return the default color table node id for a given type
  const char * GetDefaultColorTableNodeID(int type);

  /// 
  /// Return the default freesurfer color node id for a given type
  const char * GetDefaultFreeSurferColorNodeID(int type);

  /// 
  /// Return the default dGEMRIC color node id for a given type
  const char * GetDefaultdGEMRICColorNodeID(int type);

  /// 
  /// Return the default PET color node id for a given type
  const char * GetDefaultPETColorNodeID(int type);

  /// 
  /// Return a default color node id for a freesurfer label map volume
  const char * GetDefaultFreeSurferLabelMapColorNodeID();
  
  /// 
  /// Return a default color node id for a volume
  const char * GetDefaultVolumeColorNodeID();

  /// 
  /// Return a default color node id for a label map
  const char * GetDefaultLabelMapColorNodeID();

  /// 
  /// Return a default color node id for a model
  const char * GetDefaultModelColorNodeID();

  /// 
  /// return a default color node id for a procedural color node
  char * GetDefaultProceduralColorNodeID(const char *name);

  /// 
  /// return a default color node id for a file based node, based on the file name
  char * GetDefaultFileColorNodeID(const char *fileName);

  /// 
  /// look for color files in the Base/Logic/Resources/ColorFiles directory and
  /// put their names in the ColorFiles list
  void FindColorFiles();

  /// 
  /// Add a file to the ColorFiles list, checking first for null, duplicates
  void AddColorFile(const char *fileName);
  
  /// 
  /// load in a color file and add the name to the names list, creating a
  /// storage node
  vtkMRMLColorNode * LoadColorFile(const char *fileName, const char *nodeName = NULL);

  /// 
  /// Get/Set the user defined paths where to look for extra colour files
  vtkGetStringMacro(UserColorFilePaths);
  vtkSetStringMacro(UserColorFilePaths);

protected:
  vtkSlicerColorLogic();
  ~vtkSlicerColorLogic();
  vtkSlicerColorLogic(const vtkSlicerColorLogic&);
  void operator=(const vtkSlicerColorLogic&);

  /// 
  /// a vector holding discovered colour files, white space separated with:
  /// int name r g b a
  /// with rgba in the range 0-255
  //BTX
  std::vector<std::string> ColorFiles;
  //ETX
  /// 
  /// a string holding delimiter separated (; on win32, : else) paths where to
  /// look for extra colour files, set from the return value of
  /// vtkSlicerApplication::GetColorFilePaths
  char *UserColorFilePaths;
};

#endif

