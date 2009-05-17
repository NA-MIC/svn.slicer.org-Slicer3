#ifndef __vtkIGTPlanningLoadingPreoperativeDataStep_h
#define __vtkIGTPlanningLoadingPreoperativeDataStep_h

#include "vtkIGTPlanningStep.h"

class vtkKWFrameWithLabel;
class vtkKWMenuButtonWithLabel;
class vtkKWMenuButtonWithLabel;
class vtkKWPushButton;

class VTK_IGT_EXPORT vtkIGTPlanningLoadingPreoperativeDataStep : public vtkIGTPlanningStep
{
public:
  static vtkIGTPlanningLoadingPreoperativeDataStep *New();
  vtkTypeRevisionMacro(vtkIGTPlanningLoadingPreoperativeDataStep,vtkIGTPlanningStep);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Reimplement the superclass's method (see vtkKWWizardStep).
  virtual void ShowUserInterface();

protected:
  vtkIGTPlanningLoadingPreoperativeDataStep();
  ~vtkIGTPlanningLoadingPreoperativeDataStep();

  virtual void PopulatePreoperativeImageDataSelector();
  virtual void PopulateToolModelSelector();
  void ProcessGUIEvents(vtkObject *caller, unsigned long event, void *callData);

  vtkKWMenuButtonWithLabel   *PreoperativeImageDataMenuButton; 
  vtkKWMenuButtonWithLabel   *ToolModelMenuButton; 

  vtkKWPushButton            *AddVolumeButton; 

private:
  vtkIGTPlanningLoadingPreoperativeDataStep(const vtkIGTPlanningLoadingPreoperativeDataStep&);
  void operator=(const vtkIGTPlanningLoadingPreoperativeDataStep&);
};

#endif
