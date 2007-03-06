/*=auto=========================================================================

  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See Doc/copyright/copyright.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: vtkSlicerColorLogic.h,v $
  Date:      $Date: 2006/01/08 04:48:05 $
  Version:   $Revision: 1.45 $

=========================================================================auto=*/

// .NAME vtkSlicerColorLogic - slicer logic class for color manipulation
// .SECTION Description
// This class manages the logic associated with reading, saving,
// and changing propertied of the colors


#ifndef __vtkSlicerColorLogic_h
#define __vtkSlicerColorLogic_h

#include <stdlib.h>

#include "vtkSlicerBaseLogic.h"
#include "vtkSlicerLogic.h"

#include "vtkMRML.h"

class VTK_SLICER_BASE_LOGIC_EXPORT vtkSlicerColorLogic : public vtkSlicerLogic 
{
  public:
  
  // The Usual vtk class functions
  static vtkSlicerColorLogic *New();
  vtkTypeRevisionMacro(vtkSlicerColorLogic,vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Update logic state when MRML scene chenges
  void ProcessMRMLEvents(vtkObject * caller,
                         unsigned long event,
                         void * callData);

  // Description:
  // Add a series of color nodes, setting the types to the defaults, so that
  // they're accessible to the rest of Slicer
  void AddDefaultColorNodes();

  // Description:
  // Remove the colour nodes that were added
  void RemoveDefaultColorNodes();

  // Description:
  // Remove the colour nodes that were added
  void RemoveDefaultColorNodesFromScene();

  // Description:
  // Return the default color table node id for a given type
  const char * GetDefaultColorTableNodeID(int type);

  // Description:
  // Return the default freesurfer color node id for a given type
  const char * GetDefaultFreeSurferColorNodeID(int type);

  // Description:
  // Return a default color node id for a freesurfer label map volume
  const char *GetDefaultFreeSurferLabelMapColorNodeID();
  
  // Description:
  // Return a default color node id for a volume
  const char * GetDefaultVolumeColorNodeID();

  // Description:
  // Return a default color node id for a label map
  const char * GetDefaultLabelMapColorNodeID();

  // Description:
  // Return a default color node id for a model
  const char * GetDefaultModelColorNodeID();
  
protected:
  vtkSlicerColorLogic();
  ~vtkSlicerColorLogic();
  vtkSlicerColorLogic(const vtkSlicerColorLogic&);
  void operator=(const vtkSlicerColorLogic&);

};

#endif

