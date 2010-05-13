/*=auto=========================================================================

  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) 
  All Rights Reserved.

  See Doc/copyright/copyright.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer

=========================================================================auto=*/

#include "vtkMRMLLabelMapVolumeDisplayNode.h"

#include <stdlib.h>
#include <iostream>

#include "TestingMacros.h"

int vtkMRMLLabelMapVolumeDisplayNodeTest1(int , char * [] )
{
  vtkSmartPointer< vtkMRMLLabelMapVolumeDisplayNode > node1 = vtkSmartPointer< vtkMRMLLabelMapVolumeDisplayNode >::New();

  EXERCISE_BASIC_OBJECT_METHODS( node1 );

  EXERCISE_BASIC_DISPLAY_MRML_METHODS(vtkMRMLLabelMapVolumeDisplayNode, node1);
  
  return EXIT_SUCCESS;
}
