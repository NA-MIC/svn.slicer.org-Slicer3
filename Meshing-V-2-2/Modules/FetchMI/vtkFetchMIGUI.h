#ifndef __vtkFetchMIGUI_h
#define __vtkFetchMIGUI_h

#include "vtkSlicerBaseGUIWin32Header.h"
#include "vtkSlicerModuleGUI.h"

#include "vtkMRMLScene.h"
#include "vtkFetchMILogic.h"
#include "vtkMRMLFetchMINode.h"
#include "vtkIntArray.h"


class vtkSlicerNodeSelectorWidget;
class vtkFetchMIQueryTermWidget;
class vtkFetchMIFlatResourceWidget;
class vtkFetchMIResourceUploadWidget;
class vtkFetchMITagViewWidget;
class vtkXNDHandler;

class vtkKWFrame;
class vtkKWNotebook;
class vtkKWLabel;
class vtkKWFrameWithLabel;
class vtkKWPushButton;
class vtkKWEntry;
class vtkKWEntryWithLabel;
class vtkKWTopLevel;
class vtkKWMenuButton;
class vtkFetchMIIcons;


class VTK_FETCHMI_EXPORT vtkFetchMIGUI : public vtkSlicerModuleGUI
{
  public:
  static vtkFetchMIGUI *New();
  vtkTypeMacro(vtkFetchMIGUI,vtkSlicerModuleGUI);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the logic pointer from parent class pointer.
  // Overloads implementation in vtkSlicerModulesGUI
  // to allow loadable modules.
  virtual void SetModuleLogic ( vtkSlicerLogic *logic )
  {
  this->SetLogic(reinterpret_cast<vtkFetchMILogic*> (logic)); 
  }

   // Description: Get/Set MRML node
  vtkGetObjectMacro (Logic, vtkFetchMILogic);
  vtkSetObjectMacro (Logic, vtkFetchMILogic);
  
  // Description: Get/Set MRML node
  vtkGetObjectMacro (FetchMINode, vtkMRMLFetchMINode);
  vtkSetObjectMacro (FetchMINode, vtkMRMLFetchMINode);

  vtkGetObjectMacro ( QueryList, vtkFetchMIQueryTermWidget );
  vtkGetObjectMacro ( ResourceList, vtkFetchMIFlatResourceWidget );
  vtkGetObjectMacro ( TaggedDataList, vtkFetchMIResourceUploadWidget );
  vtkGetObjectMacro ( AddServerButton, vtkKWPushButton );
  vtkGetObjectMacro ( ServerMenuButton, vtkKWMenuButton );
  vtkGetObjectMacro ( NewServerLabel, vtkKWLabel );
  vtkGetObjectMacro ( RefreshButton, vtkKWPushButton );
  vtkGetObjectMacro ( AddServerEntry, vtkKWEntry );
  vtkGetObjectMacro ( FetchMIIcons, vtkFetchMIIcons );
  vtkGetObjectMacro ( QueryTagsButton, vtkKWPushButton );
  vtkGetObjectMacro ( TagViewer, vtkFetchMITagViewWidget);
  vtkGetObjectMacro ( ServerTypeMenuButton, vtkKWMenuButton );
  vtkGetObjectMacro ( CloseNewServerButton, vtkKWPushButton );
  vtkGetObjectMacro ( NewServerWindow, vtkKWTopLevel );
  vtkGetObjectMacro ( Notebook, vtkKWNotebook );

  vtkGetMacro (GUIWidth, int );
  vtkSetMacro (GUIWidth, int );

  // Description:
  // Create widgets

  virtual void BuildGUI ( );
  virtual void TearDownGUI ( );
  void LoadTclPackage ( );

  vtkIntArray *NewObservableEvents();

  // Description:
  // Add obsereves to GUI widgets
  virtual void AddGUIObservers ( );

  // Description:
  // Remove obsereves to GUI widgets
  virtual void RemoveGUIObservers ( );
  virtual void RemoveMRMLNodeObservers ( );
  virtual void RemoveLogicObservers ( );
  
  // Description:
  // Pprocess events generated by Logic
  virtual void ProcessLogicEvents ( vtkObject *caller, unsigned long event,
                                  void *callData ){};

  // Description:
  // Pprocess events generated by GUI widgets
  virtual void ProcessGUIEvents ( vtkObject *caller, unsigned long event,
                                  void *callData );

  // Description:
  // Pprocess events generated by MRML
  virtual void ProcessMRMLEvents ( vtkObject *caller, unsigned long event, 
                                  void *callData);
  // Description:
  // Pops up a window into which the name of a new server
  // can be added. There is currently no check for whether
  // the server is valid. 
  void RaiseNewServerWindow();
  void WithdrawNewServerWindow();
  void DestroyNewServerWindow();

  // Description:
  // Describe behavior at module startup and exit.
  virtual void Enter ( );
  virtual void Exit ( );
  virtual void Init ( );

  // Description:
  // Get the categorization of the module.  The category is used for
  // grouping modules together into menus.
  const char *GetCategory() const {return "Informatics";}

protected:
  vtkFetchMIGUI();
  ~vtkFetchMIGUI();
  vtkFetchMIGUI(const vtkFetchMIGUI&);
  void operator=(const vtkFetchMIGUI&);

  // Description:
  // Updates GUI widgets based on parameters values in MRML node
  void UpdateGUI();
  void UpdateTagTableFromGUI();
  void PopulateQueryListFromServer();
  
  // Description:
  // If the selected server is re-selected, then the
  // appropriate tagtable's RestoreSelectionState variable is
  // set, and Logic repopulates its list of available Tags and
  // values for query. This method restores
  // the values that a user had previously selected
  // for tags in the GUI as a user convenience.
  void RestoreSelectedValuesForTagsFromMRML();
  
  // Description:
  //--- Method restores user's selection state for all tags.
  void UpdateTagTableFromMRML();
  void UpdateResourceTableFromMRML();
  void InitializeSceneTable();
  void UpdateSceneTableFromMRML();
  void AddMRMLSceneRow();
  void UpdateVolumeNodes();
  void UpdateModelNodes();
  void UpdateFiducialListNodes();
  void UpdateUnstructuredGridNodes();
  void UpdateColorTableNodes();

  // Description:
  // Apply/remove current  tag to selected data.
  void TagSelectedData();
  void RemoveTagFromSelectedData();

  // Description:
  //--- This method loops through all resources selected in the ResourceList.
  //--- For each selected resource, the DeleteResourceFromServer method on the Logic
  //--- is called, which in turn requests the DeleteResource method on the appropriate handler.
  //--- Then handler returns 0 if an error occurs, and the logic passes this back.
  //--- For each resource successfully deleted, the ResourceList is updated; and if any
  //--- resource deletion causes an error, a message dialog is posted warning the user
  //--- that not all resources may have been deleted.
  void DeleteSelectedResourcesFromServer();

  void ShowAllTagView();
  void ShowSelectionTagView();

  // Description:
  // updates status text in the slicer window.
  void SetStatusText( const char *txt);


  vtkFetchMILogic *Logic;
  vtkMRMLFetchMINode* FetchMINode;
  vtkFetchMIQueryTermWidget *QueryList;
  vtkFetchMIFlatResourceWidget *ResourceList;
  vtkFetchMIResourceUploadWidget *TaggedDataList;
  vtkKWPushButton *AddServerButton;
  vtkKWMenuButton *ServerTypeMenuButton;
  vtkKWPushButton *CloseNewServerButton;
  vtkKWTopLevel *NewServerWindow;
  vtkKWLabel *NewServerLabel;
  vtkKWMenuButton *ServerMenuButton;
  vtkKWPushButton *RefreshButton;
  vtkKWEntry *AddServerEntry;
  vtkKWPushButton *QueryTagsButton;
  vtkFetchMITagViewWidget *TagViewer;
  vtkFetchMIIcons *FetchMIIcons;
  vtkKWNotebook *Notebook;

  int UpdatingGUI;
  int UpdatingMRML;
  int GUIWidth;
  // Description:
  // This flag gets set when the module is entered,
  // and unset when the module is exited. It is used
  // to put observers on and off the MRML Scene.
  bool Raised;

};

#endif

