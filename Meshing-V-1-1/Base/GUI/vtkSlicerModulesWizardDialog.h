#ifndef __vtkSlicerModulesWizardDialog_h
#define __vtkSlicerModulesWizardDialog_h

#include "vtkKWWizardDialog.h"

class vtkSlicerRepositoryStep;
class vtkSlicerModulesStep;
class vtkSlicerProgressStep;

class vtkSlicerModulesWizardDialog : public vtkKWWizardDialog
{
public:
  static vtkSlicerModulesWizardDialog* New();
  vtkTypeRevisionMacro(vtkSlicerModulesWizardDialog,vtkKWWizardDialog);

  // Description:
  // Access to the steps.
  vtkGetObjectMacro(RepositoryStep, vtkSlicerRepositoryStep);
  vtkGetObjectMacro(ModulesStep, vtkSlicerModulesStep);
  vtkGetObjectMacro(ProgressStep, vtkSlicerProgressStep);

protected:
  vtkSlicerModulesWizardDialog();
  ~vtkSlicerModulesWizardDialog() {};

  // Description:
  // Create the widget.
  virtual void CreateWidget();

  // Description:
  // Steps
  vtkSlicerRepositoryStep *RepositoryStep;
  vtkSlicerModulesStep *ModulesStep;
  vtkSlicerProgressStep *ProgressStep;

private:
  vtkSlicerModulesWizardDialog(const vtkSlicerModulesWizardDialog&);   // Not implemented.
  void operator=(const vtkSlicerModulesWizardDialog&);  // Not implemented.
};

#endif
