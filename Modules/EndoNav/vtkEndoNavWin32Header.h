/*==========================================================================

Portions (c) Copyright 2008 Brigham and Women's Hospital (BWH) All Rights Reserved.

See Doc/copyright/copyright.txt
or http://www.slicer.org/copyright/copyright.txt for details.

Program:   3D Slicer
Module:    $HeadURL: http://svn.slicer.org/Slicer3/branches/EndoTracking/Modules/EndoNav/vtkEndoNavWin32Header.h $
Date:      $Date: 2009-01-05 13:28:20 -0500 (Mon, 05 Jan 2009) $
Version:   $Revision: 8267 $

==========================================================================*/

#ifndef __vtkEndoNavWin32Header_h
#define __vtkEndoNavWin32Header_h

#include <vtkEndoNavConfigure.h>

#if defined(WIN32) && !defined(VTKSLICER_STATIC)
#if defined(EndoNav_EXPORTS)
#define VTK_OPENIGTLINKIF_EXPORT __declspec( dllexport ) 
#else
#define VTK_OPENIGTLINKIF_EXPORT __declspec( dllimport ) 
#endif
#else
#define VTK_OPENIGTLINKIF_EXPORT 
#endif
#endif
