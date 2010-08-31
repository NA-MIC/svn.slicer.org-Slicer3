/*==============================================================================

  Program: 3D Slicer

  Copyright (c) 2010 Kitware Inc.

  See Doc/copyright/copyright.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
  and was partially funded by NIH grant 3P41RR013218-12S1

==============================================================================*/

#ifndef __qSlicerEMSegmentRunSegmentationPanel_h
#define __qSlicerEMSegmentRunSegmentationPanel_h 

// CTK includes
#include <ctkPimpl.h>

// EMSegment includes
#include "qSlicerEMSegmentWidget.h"

// EMSegment/MRML includes
class vtkMRMLROINode;

#include "qSlicerEMSegmentModuleExport.h"

class qSlicerEMSegmentRunSegmentationPanelPrivate;

class Q_SLICER_QTMODULES_EMSEGMENT_EXPORT qSlicerEMSegmentRunSegmentationPanel :
    public qSlicerEMSegmentWidget
{ 
  Q_OBJECT

public:

  typedef qSlicerEMSegmentWidget Superclass;
  qSlicerEMSegmentRunSegmentationPanel(QWidget *newParent=0);

  void updateWidgetFromMRML();

  void updateMRMLFromWidget();

private:

  void setMRMLROINode(vtkMRMLROINode* node);

private slots:

  void display2DVOI(bool show);

private:
  CTK_DECLARE_PRIVATE(qSlicerEMSegmentRunSegmentationPanel);
  friend class qSlicerEMSegmentRunSegmentationStep; // for access to setMRMLROINode
};

#endif
