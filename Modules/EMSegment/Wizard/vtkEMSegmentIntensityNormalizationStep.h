/*=auto=======================================================================

  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights
  Reserved.

  See Doc/copyright/copyright.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: vtkEMSegmentIntensityNormalizationStep.h,v$
  Date:      $Date: 2006/01/06 17:56:51 $
  Version:   $Revision: 1.6 $
  Author:    $Nicolas Rannou (BWH), Sylvain Jaume (MIT)$

=======================================================================auto=*/

#ifndef __vtkEMSegmentIntensityNormalizationStep_h
#define __vtkEMSegmentIntensityNormalizationStep_h

#include "vtkEMSegmentStep.h"

class vtkKWMenuButton;
class vtkKWMenuButtonWithLabel;
class vtkKWScaleWithEntry;
class vtkKWFrameWithLabel;
class vtkKWEntryWithLabel;
class vtkKWCheckButtonWithLabel;
class vtkKWHistogram;
class vtkKWPiecewiseFunctionEditor;

class VTK_EMSEGMENT_EXPORT vtkEMSegmentIntensityNormalizationStep :
  public vtkEMSegmentStep
{
public:
  static vtkEMSegmentIntensityNormalizationStep *New();
  vtkTypeRevisionMacro(vtkEMSegmentIntensityNormalizationStep,
      vtkEMSegmentStep);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Reimplement the superclass's method (see vtkKWWizardStep).
  virtual void ShowUserInterface();

  // Description:
  // Callbacks.
  virtual void NormalizationTargetSelectionChangedCallback(vtkIdType VolId);
  virtual void NormalizationEnableCallback(vtkIdType VolId, int state);
  virtual void NormalizationNormTypeCallback(vtkIdType VolId, int enumType);
  virtual void NormalizationPrintInfoCallback(vtkIdType VolId, int checked);
  virtual void NormalizationValueCallback(vtkIdType VolId, double dValue);
  virtual void NormalizationSmoothingWidthCallback(vtkIdType VolId, int
      iValue);
  virtual void NormalizationMaxSmoothingWidthCallback(vtkIdType VolId, int
      iValue);
  virtual void NormalizationRelativeMaxVoxelNumCallback(vtkIdType VolId,
      double dValue);

  virtual void NormalizationHistogramChangedCallback(vtkIdType VolId);

  //Get histogram Value
  virtual void GetHistogramValue();

  // Description:
  // Observers
  virtual void AddCursorMovingGUIEvents();
  virtual void RemoveCursorMovingGUIEvents();
  virtual void ProcessCursorMovingGUIEvents(vtkObject *caller, unsigned long
      event, void *callData);

  //BTX
  enum
    {
    NormalizationDefaultT1SPGR = 0,
    NormalizationDefaultT2
    };
  //ETX

  // Description:
  // Callbacks.
  virtual void MaskTargetSelectionChangedCallback(vtkIdType VolId);
  //virtual void MaskEnableCallback(vtkIdType VolId, int state);
  //virtual void MaskNormTypeCallback(vtkIdType VolId, int enumType);
  //virtual void MaskPrintInfoCallback(vtkIdType VolId, int checked);
  //virtual void MaskValueCallback(vtkIdType VolId, double dValue);

  virtual void MaskHistogramChangedCallback(vtkIdType VolId);

  // Get histogram value
  //virtual void GetMaskHistogramValue();

protected:
  vtkEMSegmentIntensityNormalizationStep();
  ~vtkEMSegmentIntensityNormalizationStep();

  virtual void PopulateNormalizationTargetVolumeSelector();
  virtual void PopulateNormalizationHistogramSelector();

  virtual void PopulateMaskTargetVolumeSelector();
  virtual void PopulateMaskHistogramSelector();

  virtual void ResetDefaultParameters(vtkIdType targetVolId);
  virtual void HideUserInterface();

  vtkKWMenuButtonWithLabel     *NormalizationTargetVolumeMenuButton;
  vtkKWFrameWithLabel          *NormalizationParametersFrame;
  vtkKWCheckButtonWithLabel    *NormalizationEnableCheckButton;
  vtkKWMenuButton              *NormalizationDefaultsMenuButton;
  vtkKWCheckButtonWithLabel    *NormalizationPrintCheckButton;
  vtkKWEntryWithLabel          *NormalizationNormValueEntry;
  vtkKWEntryWithLabel          *NormalizationSmoothingWidthEntry;
  vtkKWEntryWithLabel          *NormalizationMaxSmoothingWidthEntry;
  vtkKWScaleWithEntry          *NormalizationRelativeMaxVoxelScale;

  vtkKWMenuButtonWithLabel     *MaskHistogramMenuButton;
  vtkKWHistogram               *MaskHistogram;
  vtkKWFrameWithLabel          *MaskHistogramFrame;
  vtkKWPiecewiseFunctionEditor *MaskPiecewiseFunctionEditor;

  vtkKWMenuButtonWithLabel     *NormalizationHistogramMenuButton;
  vtkKWHistogram               *NormalizationHistogram;
  vtkKWFrameWithLabel          *NormalizationHistogramFrame;
  vtkKWPiecewiseFunctionEditor *NormalizationPiecewiseFunctionEditor;

  vtkKWEntryWithLabel          *NormalizationValueRecommendedEntry;
  vtkKWFrameWithLabel          *RecommendationFrame;

  vtkIdType IdEvent;

private:
  vtkEMSegmentIntensityNormalizationStep(const
      vtkEMSegmentIntensityNormalizationStep&);
  void operator=(const vtkEMSegmentIntensityNormalizationStep&);
};

#endif

