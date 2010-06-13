#ifndef __qSlicerScalarVolumeDisplayWidget_h
#define __qSlicerScalarVolumeDisplayWidget_h

// Qt includes
#include <QWidget>

// CTK includes
#include <ctkPimpl.h>
#include <ctkVTKObject.h>

// SlicerQt includes
#include <qSlicerWidget.h>

#include "qSlicerVolumesModuleExport.h"

class vtkMRMLNode;
class vtkMRMLScalarVolumeNode;
class qSlicerScalarVolumeDisplayWidgetPrivate;

class Q_SLICER_QTMODULES_VOLUMES_EXPORT qSlicerScalarVolumeDisplayWidget : public qSlicerWidget
{
  Q_OBJECT

public:
  /// Constructors
  typedef qSlicerWidget Superclass;
  explicit qSlicerScalarVolumeDisplayWidget(QWidget* parent);
  virtual ~qSlicerScalarVolumeDisplayWidget(){}


public slots:

  /// 
  /// Set the MRML node of interest
  void setMRMLVolumeNode(vtkMRMLScalarVolumeNode* volumeNode);
  void setMRMLVolumeNode(vtkMRMLNode* node);


protected slots:

protected:

private:
  CTK_DECLARE_PRIVATE(qSlicerScalarVolumeDisplayWidget);
};

#endif
