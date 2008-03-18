/*=auto=========================================================================

  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See Doc/copyright/copyright.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: vtkSlicerApplicationGUI.cxx,v $
  Date:      $Date: 2006/01/08 04:48:05 $
  Version:   $Revision: 1.45 $

=========================================================================auto=*/

#include <sstream>
#include <string>
#include <vtksys/SystemTools.hxx> 
#include <itksys/SystemTools.hxx> 

#include "vtkCommand.h"
#include "vtkCornerAnnotation.h"
#include "vtkObjectFactory.h"
#include "vtkToolkits.h"

// things for temporary MainViewer display.
#include "vtkCubeSource.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkCamera.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderWindow.h"

#include "vtkKWApplication.h"
#include "vtkKWTclInteractor.h"
#include "vtkKWWidget.h"
#include "vtkKWFrame.h"
#include "vtkKWMenu.h"
#include "vtkKWMenuButton.h"
#include "vtkKWNotebook.h"
#include "vtkKWRenderWidget.h"
#include "vtkKWUserInterfacePanel.h"
#include "vtkKWResourceUtilities.h"
#include "vtkKWSplitFrame.h"
#include "vtkKWUserInterfaceManagerNotebook.h"
#include "vtkKWMessageDialog.h"
#include "vtkKWToolbarSet.h"
#include "vtkKWMessageDialog.h"
#include "vtkKWProgressDialog.h"
#include "vtkKWEntry.h"
#include "vtkKWLabel.h"

#include "vtkSlicerWindow.h"
#include "vtkSlicerApplication.h"
#include "vtkSlicerApplicationGUI.h"
#include "vtkSlicerApplicationLogic.h"
#include "vtkSlicerModuleGUI.h"
#include "vtkSlicerModuleNavigator.h"
#include "vtkSlicerGUILayout.h"
#include "vtkSlicerTheme.h"
#include "vtkSlicerFont.h"
#include "vtkSlicerColor.h"
#include "vtkSlicerMRMLSaveDataWidget.h"
#include "vtkSlicerApplicationSettingsInterface.h"
#include "vtkSlicerSliceControllerWidget.h"
#include "vtkSlicerViewerInteractorStyle.h"

#include "vtkSlicerSlicesGUI.h"

#include "vtkSlicerFiducialListWidget.h"

#ifdef USE_PYTHON
#ifdef _DEBUG
#undef _DEBUG
#include <Python.h>
#define _DEBUG
#else
#include <Python.h>
#endif
#endif

//---------------------------------------------------------------------------
vtkStandardNewMacro (vtkSlicerApplicationGUI);
vtkCxxRevisionMacro(vtkSlicerApplicationGUI, "$Revision: 1.0 $");
vtkCxxSetObjectMacro(vtkSlicerApplicationGUI, SlicesGUI, vtkSlicerSlicesGUI);

// temporary crud for vtkDebugLeak hunting. Will remove
// these and other related #ifndefs-#endifs throughout.

//#define LOGODISPLAY_DEBUG
//#define TOOLBAR_DEBUG
//#define VIEWCONTROL_DEBUG
//#define SLICEVIEWER_DEBUG

//#define MENU_DEBUG
//#define SLICESCONTROL_DEBUG
//#define MODULECHOOSE_DEBUG

//---------------------------------------------------------------------------
vtkSlicerApplicationGUI::vtkSlicerApplicationGUI (  )
{

  this->MRMLScene = NULL;

    //---  
    // widgets used in the Slice module
    //---
  
    //--- slicer main window
    this->MainSlicerWindow = vtkSlicerWindow::New ( );

    // Frames that comprise the Main Slicer GUI

    this->TopFrame = vtkKWFrame::New();
    this->LogoFrame = vtkKWFrame::New();
    this->SlicesControlFrame = vtkSlicerModuleCollapsibleFrame::New();
    this->ViewControlFrame = vtkSlicerModuleCollapsibleFrame::New();
    this->DropShadowFrame = vtkKWFrame::New();
    this->GridFrame1 = vtkKWFrame::New ( );
    this->GridFrame2 = vtkKWFrame::New ( );

    // initialize in case any are not defined.
    this->ApplicationToolbar = NULL;
    this->ViewControlGUI = NULL;
    this->SlicesControlGUI = NULL;
    this->LogoDisplayGUI = NULL;
    
    //--- GUIs containing components packed inside the Frames
#ifndef TOOLBAR_DEBUG
    this->ApplicationToolbar = vtkSlicerToolbarGUI::New ( );
#endif
#ifndef VIEWCONTROL_DEBUG
    this->ViewControlGUI = vtkSlicerViewControlGUI::New ( );
#endif
#ifndef SLICESCONTROL_DEBUG
    this->SlicesControlGUI = vtkSlicerSlicesControlGUI::New ( );
#endif
#ifndef LOGODISPLAY_DEBUG    
    this->LogoDisplayGUI = vtkSlicerLogoDisplayGUI::New ( );
#endif
    
    //--- Main viewer
    this->ViewerWidget = NULL;
    this->FiducialListWidget = NULL;

    // use STL::MAP to hold all main slice viewers
    this->SlicesGUI = NULL;

    //--- Save and load scene dialogs, widgets
    this->LoadSceneDialog = vtkKWLoadSaveDialog::New();

    this->SaveDataWidget = vtkSlicerMRMLSaveDataWidget::New();

    //--- unique tag used to mark all view notebook pages
    //--- so that they can be identified and deleted when 
    //--- viewer is reformatted.
    this->ViewerPageTag = 1999;
}



//---------------------------------------------------------------------------
vtkSlicerApplicationGUI::~vtkSlicerApplicationGUI ( )
{
    if (this->SaveDataWidget)
      {
      this->SaveDataWidget->SetParent(NULL);
      this->SaveDataWidget->Delete();
      this->SaveDataWidget=NULL;
      }

    // TODO: fix this method to remove all actors, etc.
    // from Main Viewer, slice viewers before deleting.
    this->DestroyMain3DViewer ( );
    this->DestroyMainSliceViewers ( );

    // Delete frames
    if ( this->TopFrame )
      {
      this->TopFrame->SetParent ( NULL );
      this->TopFrame->Delete ( );
      this->TopFrame = NULL;
      }
    if ( this->LogoFrame )
      {
      this->LogoFrame->SetParent ( NULL );
      this->LogoFrame->Delete ();
      this->LogoFrame = NULL;
      }
    if ( this->DropShadowFrame )
      {
      this->DropShadowFrame->SetParent ( NULL );
      this->DropShadowFrame->Delete ( );
      this->DropShadowFrame = NULL;
      }
    if ( this->SlicesControlFrame )
      {
      this->SlicesControlFrame->SetParent ( NULL );
      this->SlicesControlFrame->Delete ( );
      this->SlicesControlFrame = NULL;
      }
    if ( this->ViewControlFrame )
      {
      this->ViewControlFrame->SetParent ( NULL );
      this->ViewControlFrame->Delete ( );
      this->ViewControlFrame = NULL;
      }

    if ( this->LoadSceneDialog )
      {
      this->LoadSceneDialog->SetParent ( NULL );
      this->LoadSceneDialog->Delete();
      this->LoadSceneDialog = NULL;
      }
    if ( this->MainSlicerWindow )
      {
      if ( this->GetApplication() )
        {
        this->GetApplication()->RemoveWindow ( this->MainSlicerWindow );
        vtkSlicerApplication *app = vtkSlicerApplication::SafeDownCast( this->GetApplication() );
        if ( app->GetMainLayout() )
          {
          app->GetMainLayout()->SetMainSlicerWindow(NULL);
          }
        }
      this->MainSlicerWindow->SetParent ( NULL );
      this->MainSlicerWindow->Delete ( );
      this->MainSlicerWindow = NULL;
      }

    if (this->SaveDataWidget)
      {
      this->SaveDataWidget->SetParent ( NULL );
      this->SaveDataWidget->Delete();
      }

    if (this->SlicesGUI)
      {
        this->SlicesGUI->Delete();
        this->SlicesGUI = NULL;
      }
    
    this->SetApplication(NULL);
    this->SetApplicationLogic ( NULL );
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI:: DeleteComponentGUIs()
{
#ifndef VIEWCONTROL_DEBUG
    if ( this->ViewControlGUI )
      {
//      this->ViewControlGUI->TearDownGUI ( );
      this->ViewControlGUI->RemoveSliceEventObservers();
      this->ViewControlGUI->SetAndObserveMRMLScene ( NULL );
      this->ViewControlGUI->SetApplicationGUI ( NULL);
      this->ViewControlGUI->SetApplication ( NULL );
      this->ViewControlGUI->Delete ( );
      this->ViewControlGUI = NULL;
      }
#endif
#ifndef LOGODISPLAY_DEBUG
    if ( this->LogoDisplayGUI )
      {
      this->LogoDisplayGUI->Delete ( );
      this->LogoDisplayGUI = NULL;
      }
#endif
#ifndef SLICESCONTROL_DEBUG
    if ( this->SlicesControlGUI )
      {
      this->SlicesControlGUI->TearDownGUI ( );
      this->SlicesControlGUI->Delete ( );
      this->SlicesControlGUI = NULL;
      }
#endif
#ifndef TOOLBAR_DEBUG
    if ( this->ApplicationToolbar )
      {
      this->ApplicationToolbar->Delete ( );
      this->ApplicationToolbar = NULL;
      }
#endif
}



//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::PrintSelf ( ostream& os, vtkIndent indent )
{
    this->vtkObject::PrintSelf ( os, indent );

    os << indent << "SlicerApplicationGUI: " << this->GetClassName ( ) << "\n";
    os << indent << "MainSlicerWindow: " << this->GetMainSlicerWindow ( ) << "\n";
    // print widgets?
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::ProcessLoadSceneCommand()
{
    this->LoadSceneDialog->RetrieveLastPathFromRegistry("OpenPath");

    this->LoadSceneDialog->Invoke();
    // If a file has been selected for loading...
    const char *fileName = this->LoadSceneDialog->GetFileName();
    if ( fileName ) 
      {
        std::string fl(fileName);
        if (this->GetMRMLScene() && fl.find(".mrml") != std::string::npos ) 
          {
          vtkKWProgressDialog *progressDialog = vtkKWProgressDialog::New();
          progressDialog->SetParent( this->MainSlicerWindow );
          progressDialog->SetMasterWindow( this->MainSlicerWindow );
          progressDialog->Create();
          std::string message("Loading Scene...\n");
          message += std::string(fileName);
          progressDialog->SetMessageText( message.c_str() );
          // don't observe the scene, to avoid getting render updates
          // during load.  TODO: make a vtk-based progress bar that doesn't
          // call the tcl update method
          //progressDialog->SetObservedObject( this->GetMRMLScene() );
          progressDialog->Display();
          this->GetMRMLScene()->SetURL(fileName);
          this->GetMRMLScene()->Connect();
          progressDialog->SetParent(NULL);
          progressDialog->Delete();
          this->LoadSceneDialog->SaveLastPathToRegistry("OpenPath");
          }
        else if (this->GetMRMLScene() && fl.find(".xml") != std::string::npos ) 
          {
          this->Script ( "ImportSlicer2Scene \"%s\"", fileName);
          this->LoadSceneDialog->SaveLastPathToRegistry("OpenPath");
          }
        else if ( this->GetMRMLScene() && fl.find(".xcat") != std::string::npos )
          {
          this->Script ( "XcedeCatalogImport %s", fileName);
          this->LoadSceneDialog->SaveLastPathToRegistry("OpenPath");
          }

        if (  this->GetMRMLScene()->GetErrorCode() != 0 ) 
          {
          vtkKWMessageDialog *dialog = vtkKWMessageDialog::New();
          dialog->SetParent (  this->MainSlicerWindow );
          dialog->SetStyleToMessage();
          std::string msg = this->GetMRMLScene()->GetErrorMessage();
          dialog->SetText(msg.c_str());
          dialog->Create ( );
          dialog->Invoke();
          dialog->Delete();
          }
      }
    return;
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::ProcessImportSceneCommand()
{
    this->LoadSceneDialog->RetrieveLastPathFromRegistry(
      "OpenPath");

    this->LoadSceneDialog->Invoke();
    // If a file has been selected for loading...
    const char *fileName = this->LoadSceneDialog->GetFileName();
    if ( fileName ) 
      {
        std::string fl(fileName);
        if (this->GetMRMLScene() && fl.find(".mrml") != std::string::npos ) 
          {
          this->GetMRMLScene()->SetURL(fileName);
          this->GetMRMLScene()->Import();
          this->LoadSceneDialog->SaveLastPathToRegistry("OpenPath");
          }
        else if (this->GetMRMLScene() && fl.find(".xml") != std::string::npos ) 
          {
          this->Script ( "ImportSlicer2Scene %s", fileName);
          this->LoadSceneDialog->SaveLastPathToRegistry("OpenPath");
          }
        else if ( this->GetMRMLScene() && fl.find(".xcat") != std::string::npos )
          {
          this->Script ( "XcedeCatalogImport %s", fileName);
          this->LoadSceneDialog->SaveLastPathToRegistry("OpenPath");
          }


        if (  this->GetMRMLScene()->GetErrorCode() != 0 ) 
          {
          vtkKWMessageDialog *dialog = vtkKWMessageDialog::New();
          dialog->SetParent (  this->MainSlicerWindow );
          dialog->SetStyleToMessage();
          std::string msg = this->GetMRMLScene()->GetErrorMessage();
          dialog->SetText(msg.c_str());
          dialog->Create ( );
          dialog->Invoke();
          dialog->Delete();
          }
      }
    return;
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::ProcessAddDataCommand()
{
  this->GetApplication()->Script("::Loader::ShowDialog");
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::ProcessAddVolumeCommand()
{
  this->GetApplication()->Script("::LoadVolume::ShowDialog");
}


//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::ProcessCloseSceneCommand()
{
  vtkKWMessageDialog *dialog = vtkKWMessageDialog::New();
  dialog->SetParent ( this->MainSlicerWindow );
  dialog->SetStyleToOkCancel();
  dialog->SetText("Are you sure you want to close scene?");
  dialog->Create ( );
  if (dialog->Invoke())
  {
    if (this->GetMRMLScene()) 
      {
      //this->MRMLScene->Clear();
      this->MRMLScene->Clear(false);
      }
  }
  dialog->Delete();
}  

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::ProcessSaveSceneAsCommand()
{
    this->SaveDataWidget->SetAndObserveMRMLScene(this->GetMRMLScene());
    this->SaveDataWidget->AddObserver ( vtkSlicerMRMLSaveDataWidget::DataSavedEvent,  (vtkCommand *)this->GUICallbackCommand );
    this->SaveDataWidget->Invoke();  

    this->SaveDataWidget->RemoveObservers ( vtkSlicerMRMLSaveDataWidget::DataSavedEvent,  (vtkCommand *)this->GUICallbackCommand );
    return;
}    

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::AddGUIObservers ( )
{
        vtkSlicerApplication::SafeDownCast ( this->GetApplication() )->AddObserver ( vtkCommand::ModifiedEvent, (vtkCommand *)this->GUICallbackCommand );
        this->GetMainSlicerWindow()->GetFileMenu()->AddObserver (vtkKWMenu::MenuItemInvokedEvent, (vtkCommand *)this->GUICallbackCommand );
        // keep track of changes to Home Module set from Application Settings interface;
        // try trapping the View menu events, and just updating the home module from registry...
        this->LoadSceneDialog->AddObserver ( vtkCommand::ModifiedEvent, (vtkCommand *)this->GUICallbackCommand );
#ifndef TOOLBAR_DEBUG
        this->GetApplicationToolbar()->AddGUIObservers ( );
#endif
#ifndef VIEWCONTROL_DEBUG
        this->GetViewControlGUI()->AddGUIObservers ( );
#endif
#ifndef SLICESCONTROL_DEBUG
        this->GetSlicesControlGUI ( )->AddGUIObservers ( );
#endif
#ifndef LOGODISPLAY_DEBUG
        this->GetLogoDisplayGUI ( )->AddGUIObservers ( );
#endif

        if (this->SlicesGUI)
          {
            vtkSlicerSliceGUI *g = NULL;
            char *layoutname = NULL;
            int nSliceGUI = this->SlicesGUI->GetNumberOfSliceGUI();
            for (int i = 0; i < nSliceGUI; i++)
              {
                if (i == 0)
                  {
                    g = this->SlicesGUI->GetFirstSliceGUI();
                    layoutname = this->SlicesGUI->GetFirstSliceGUILayoutName();
                  }
                else
                  {
                    g = this->SlicesGUI->GetNextSliceGUI(layoutname);
                    layoutname = this->SlicesGUI->GetNextSliceGUILayoutName(layoutname);
                  }
                
                g->GetSliceController()->AddObserver(vtkSlicerSliceControllerWidget::ExpandEvent, 
                                                     (vtkCommand *)this->GUICallbackCommand);
                g->GetSliceController()->AddObserver(vtkSlicerSliceControllerWidget::ShrinkEvent, 
                                                     (vtkCommand *)this->GUICallbackCommand);
              }
          }

        if (this->SaveDataWidget)
          {
            this->SaveDataWidget->AddObserver ( vtkSlicerMRMLSaveDataWidget::DataSavedEvent,  (vtkCommand *)this->GUICallbackCommand );
          }
}


//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::RemoveGUIObservers ( )
{

  vtkSlicerApplication::SafeDownCast ( this->GetApplication() )->RemoveObservers ( vtkCommand::ModifiedEvent, (vtkCommand *)this->GUICallbackCommand );
    this->GetMainSlicerWindow()->GetFileMenu()->RemoveObservers ( vtkKWMenu::MenuItemInvokedEvent, (vtkCommand *)this->GUICallbackCommand );
    this->LoadSceneDialog->RemoveObservers ( vtkCommand::ModifiedEvent, (vtkCommand *) this->GUICallbackCommand );

#ifndef TOOLBAR_DEBUG
    this->GetApplicationToolbar()->RemoveGUIObservers ( );
#endif
#ifndef VIEWCONTROL_DEBUG
    this->GetViewControlGUI ( )->RemoveGUIObservers ( );
#endif
#ifndef SLICESCONTROL_DEBUG
    this->GetSlicesControlGUI ( )->RemoveGUIObservers ( );
#endif
#ifndef LOGODISPLAY_DEBUG
    this->GetLogoDisplayGUI ( )->RemoveGUIObservers ( );
#endif    

    if (this->SlicesGUI)
      {
        vtkSlicerSliceGUI *g = NULL;
        char *layoutname = NULL;
        int nSliceGUI = this->SlicesGUI->GetNumberOfSliceGUI();
        for (int i = 0; i < nSliceGUI; i++)
          {
            if (i == 0)
              {
                g = this->SlicesGUI->GetFirstSliceGUI();
                layoutname = this->SlicesGUI->GetFirstSliceGUILayoutName();
              }
            else
              {
                g = this->SlicesGUI->GetNextSliceGUI(layoutname);
                layoutname = this->SlicesGUI->GetNextSliceGUILayoutName(layoutname);
              }
            
            g->GetSliceController()->RemoveObservers(vtkSlicerSliceControllerWidget::ExpandEvent, 
                                                     (vtkCommand *)this->GUICallbackCommand);
            g->GetSliceController()->RemoveObservers(vtkSlicerSliceControllerWidget::ShrinkEvent, 
                                                     (vtkCommand *)this->GUICallbackCommand);
          }
      }
    
    this->RemoveMainSliceViewerObservers ( );
    
    if (this->SaveDataWidget)
      {
        this->SaveDataWidget->RemoveObservers ( vtkSlicerMRMLSaveDataWidget::DataSavedEvent,  (vtkCommand *)this->GUICallbackCommand );
      }
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::ProcessGUIEvents ( vtkObject *caller,
                                                 unsigned long event, void *callData )
{
        // This code is just a placeholder until the logic is set up to use properly:
        // For now, the GUI controls the GUI instead of going thru the logic...
        // TODO:
        // Actually, these events want to set "activeModule" in the logic;
        // using this->Logic->SetActiveModule ( ) which is currently commented out.
        // Observers on that logic should raise and lower the appropriate page.
        // So for now, the GUI is controlling the GUI instead of going thru the logic.
        //---

        vtkKWLoadSaveDialog *filebrowse = vtkKWLoadSaveDialog::SafeDownCast(caller);
        vtkSlicerApplication *app = vtkSlicerApplication::SafeDownCast( this->GetApplication() );
        vtkSlicerGUILayout *layout = app->GetMainLayout ( );
        vtkSlicerMRMLSaveDataWidget *saveDataWidget = vtkSlicerMRMLSaveDataWidget::SafeDownCast(caller);

        // catch changes to ApplicationSettings, and update the ApplicationSettingsInterface
        if ( event == vtkCommand::ModifiedEvent && vtkSlicerApplication::SafeDownCast(caller) == app )
        {
                this->MainSlicerWindow->GetApplicationSettingsInterface()->Update();
        }

        if (event == vtkSlicerModuleGUI::ModuleSelectedEvent) 
        {
                this->SelectModule((const char*)callData);
                return;
        }

        if (saveDataWidget == this->SaveDataWidget && event == vtkSlicerMRMLSaveDataWidget::DataSavedEvent)
        {
        }

        vtkSlicerSliceControllerWidget *controller_caller = 
                vtkSlicerSliceControllerWidget::SafeDownCast(caller);
        if (caller)
        {
                if (event == vtkSlicerSliceControllerWidget::ExpandEvent) 
                {
                        if (this->SlicesGUI)
                        {
                                vtkSlicerSliceGUI *g = NULL;
                                char *layoutname = NULL;
                                int nSliceGUI = this->SlicesGUI->GetNumberOfSliceGUI();
                                for (int i = 0; i < nSliceGUI; i++)
                                {
                                        if (i == 0)
                                        {
                                                g = this->SlicesGUI->GetFirstSliceGUI();
                                                layoutname = this->SlicesGUI->GetFirstSliceGUILayoutName();
                                        }
                                        else
                                        {
                                                g = this->SlicesGUI->GetNextSliceGUI(layoutname);
                                                layoutname = this->SlicesGUI->GetNextSliceGUILayoutName(layoutname);
                                        }

                                        g->GetSliceController()->Expand();
                                }
                        }
                }
                else if (event == vtkSlicerSliceControllerWidget::ShrinkEvent) 
                {
                        if (this->SlicesGUI)
                        {
                                vtkSlicerSliceGUI *g = NULL;
                                char *layoutname = NULL;
                                int nSliceGUI = this->SlicesGUI->GetNumberOfSliceGUI();
                                for (int i = 0; i < nSliceGUI; i++)
                                {
                                        if (i == 0)
                                        {
                                                g = this->SlicesGUI->GetFirstSliceGUI();
                                                layoutname = this->SlicesGUI->GetFirstSliceGUILayoutName();
                                        }
                                        else
                                        {
                                                g = this->SlicesGUI->GetNextSliceGUI(layoutname);
                                                layoutname = this->SlicesGUI->GetNextSliceGUILayoutName(layoutname);
                                        }

                                        g->GetSliceController()->Shrink();
                                }
                        }
                }
        }

        // check to see if any caller belongs to SliceControllerGUI
        // and if so, call this->SliceControllerGUI->ProcessGUIEvents ( )
        // with same params. hand it off.

        // check to see if any caller belongs to ViewControlGUI
        // and if so, call this->ViewControlGUI->ProcessGUIEvents ( )
        // with same params. hand it off.


}


//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::ProcessLogicEvents ( vtkObject *caller,
                                                   unsigned long event, void *callData )
{
    // Fill in
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::ProcessMRMLEvents ( vtkObject *caller,
                                                  unsigned long event, void *callData )
{
    // Fill in
}


//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::Enter ( )
{
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::Exit ( )
{
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::SelectModule ( const char *moduleName )
{
#ifndef TOOLBAR_DEBUG
  this->GetApplicationToolbar()->GetModuleChooseGUI()->SelectModule(moduleName);
#endif
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::BuildGUI ( )
{
    int i;

    this->SaveDataWidget->SetParent ( this->MainSlicerWindow);
    this->SaveDataWidget->SetAndObserveMRMLScene(this->GetMRMLScene());
    
    // Set up the conventional window: 3Dviewer, slice widgets, UI panel for now.
    if ( this->GetApplication() != NULL ) {

        vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();
        vtkSlicerGUILayout *layout = app->GetMainLayout ( );
        
        // Set a pointer to the MainSlicerWindow in vtkSlicerGUILayout, and
        // Set default sizes for all main frames (UIpanel and viewers) in GUI
        layout->SetMainSlicerWindow ( this->MainSlicerWindow );
        layout->InitializeLayoutDimensions ( app->GetApplicationWindowWidth(),
                                             app->GetApplicationWindowHeight(),
                                             app->GetApplicationSlicesFrameHeight());

        if ( this->MainSlicerWindow != NULL ) {

            // set up Slicer's main window
            this->MainSlicerWindow->SecondaryPanelVisibilityOn ( );
            this->MainSlicerWindow->MainPanelVisibilityOn ( );
            app->AddWindow ( this->MainSlicerWindow );

            // Create the console before the window
            // - this will make the console independent of the main window
            //   so it can be raised/lowered independently
            this->MainSlicerWindow->GetTclInteractor()->SetApplication(app);
            this->MainSlicerWindow->GetTclInteractor()->Create();

            // TODO: it would be nice to make this a menu option on the tkcon itself,
            // but for now just up the font size

            this->MainSlicerWindow->Create ( );        
            
//            app->GetTclInteractor()->SetFont("Courier 12");

            // configure initial GUI layout
            layout->InitializeMainSlicerWindowSize ( );
            layout->ConfigureMainSlicerWindowPanels ( );

            // Build main GUI frames and components that fill them
            this->BuildGUIFrames ( );

            // Build Logo GUI panel
#ifndef LOGODISPLAY_DEBUG
            vtkSlicerLogoDisplayGUI *logos = this->GetLogoDisplayGUI ( );
            logos->SetApplicationGUI ( this );
            logos->SetApplication ( app );
            logos->BuildGUI ( this->LogoFrame );
#endif            
            // Build toolbar
#ifndef TOOLBAR_DEBUG
            vtkSlicerToolbarGUI *appTB = this->GetApplicationToolbar ( );
            appTB->SetApplicationGUI ( this );
            appTB->SetApplication ( app );
            appTB->SetApplicationLogic ( this->GetApplicationLogic());
            appTB->BuildGUI ( );
#endif

            // Build SlicesControl panel
#ifndef SLICESCONTROL_DEBUG            
            vtkSlicerSlicesControlGUI *scGUI = this->GetSlicesControlGUI ( );
            scGUI->SetApplicationGUI ( this );
            scGUI->SetApplication ( app );
            scGUI->SetAndObserveMRMLScene ( this->MRMLScene );
            scGUI->BuildGUI ( this->SlicesControlFrame->GetFrame() );
#endif

            // Build 3DView Control panel
#ifndef VIEWCONTROL_DEBUG
            vtkSlicerViewControlGUI *vcGUI = this->GetViewControlGUI ( );
            vcGUI->SetApplicationGUI ( this );
            vcGUI->SetApplication ( app );
            vcGUI->SetAndObserveMRMLScene ( this->MRMLScene );
            vcGUI->BuildGUI ( this->ViewControlFrame->GetFrame() );
#endif

            this->MainSlicerWindow->GetMainNotebook()->SetUseFrameWithScrollbars ( 1 );
            this->MainSlicerWindow->GetMainNotebook()->SetEnablePageTabContextMenu ( 0 );
            
            // Build 3DViewer and Slice Viewers

#ifndef SLICEVIEWER_DEBUG
            // restore view layout from application registry...
            //this->BuildMainViewer ( app->GetApplicationLayoutType());
            this->BuildMainViewer (vtkSlicerGUILayout::SlicerLayoutDefaultView );
#endif

            // after SliceGUIs are created, the ViewControlGUI
            // needs to observe them to feed its magnifier
            // Zoom Widget.

#ifndef MENU_DEBUG
            // Construct menu bar and set up global key bindings
            // 
            // File Menu
            //
            i = this->GetMainSlicerWindow()->GetFileMenu()->InsertCommand (
                      this->GetMainSlicerWindow()->GetFileMenuInsertPosition(),
                                      "Load Scene...", this, "ProcessLoadSceneCommand");
            this->MainSlicerWindow->GetFileMenu()->SetItemAccelerator ( i, "Ctrl-O");
            this->MainSlicerWindow->GetFileMenu()->SetBindingForItemAccelerator ( i, this->MainSlicerWindow);

            this->GetMainSlicerWindow()->GetFileMenu()->InsertCommand (
                      this->GetMainSlicerWindow()->GetFileMenuInsertPosition(),
                                      "Import Scene...", this, "ProcessImportSceneCommand");

            i = this->GetMainSlicerWindow()->GetFileMenu()->InsertCommand (
                      this->GetMainSlicerWindow()->GetFileMenuInsertPosition(),
                                      "Add Data...", this, "ProcessAddDataCommand");
            this->MainSlicerWindow->GetFileMenu()->SetItemAccelerator ( i, "Ctrl-A");
            this->MainSlicerWindow->GetFileMenu()->SetBindingForItemAccelerator ( i, this->MainSlicerWindow);

            i = this->GetMainSlicerWindow()->GetFileMenu()->InsertCommand (
                      this->GetMainSlicerWindow()->GetFileMenuInsertPosition(),
                                      "Add Volume...", this, "ProcessAddVolumeCommand");

            i = this->GetMainSlicerWindow()->GetFileMenu()->InsertCommand (this->GetMainSlicerWindow()->GetFileMenuInsertPosition(),
                                               "Save", this, "ProcessSaveSceneAsCommand");
            this->MainSlicerWindow->GetFileMenu()->SetItemAccelerator ( i, "Ctrl-S");
            this->MainSlicerWindow->GetFileMenu()->SetBindingForItemAccelerator ( i, this->MainSlicerWindow);

            i = this->GetMainSlicerWindow()->GetFileMenu()->InsertCommand (this->GetMainSlicerWindow()->GetFileMenuInsertPosition(),
                                               "Close Scene", this, "ProcessCloseSceneCommand");
            this->MainSlicerWindow->GetFileMenu()->SetItemAccelerator ( i, "Ctrl-W");
            this->MainSlicerWindow->GetFileMenu()->SetBindingForItemAccelerator ( i, this->MainSlicerWindow);

            this->GetMainSlicerWindow()->GetFileMenu()->InsertSeparator (
                this->GetMainSlicerWindow()->GetFileMenuInsertPosition());

            // don't need the 'close command'
            this->GetMainSlicerWindow()->GetFileMenu()->DeleteItem ( 
                this->GetMainSlicerWindow()->GetFileMenu()->GetIndexOfItem(
                    this->GetMainSlicerWindow()->GetFileCloseMenuLabel()));
            //
            // Edit Menu
            //
            i = this->MainSlicerWindow->GetEditMenu()->AddCommand ("Set Home", NULL, "$::slicer3::ApplicationGUI SetCurrentModuleToHome");
            this->MainSlicerWindow->GetEditMenu()->SetItemAccelerator ( i, "Ctrl+H");
            this->MainSlicerWindow->GetEditMenu()->SetBindingForItemAccelerator ( i, this->MainSlicerWindow);
            i = this->MainSlicerWindow->GetEditMenu()->AddCommand ( "Undo", NULL, "$::slicer3::MRMLScene Undo" );
            this->MainSlicerWindow->GetEditMenu()->SetItemAccelerator ( i, "Ctrl+Z");
            this->MainSlicerWindow->GetEditMenu()->SetBindingForItemAccelerator ( i, this->MainSlicerWindow);
            i = this->MainSlicerWindow->GetEditMenu()->AddCommand ( "Redo", NULL, "$::slicer3::MRMLScene Redo" );
            this->MainSlicerWindow->GetEditMenu()->SetItemAccelerator ( i, "Ctrl+Y");
            this->MainSlicerWindow->GetEditMenu()->SetBindingForItemAccelerator ( i, this->MainSlicerWindow);
            
            this->GetMainSlicerWindow()->GetEditMenu()->InsertSeparator (this->GetMainSlicerWindow()->GetEditMenu()->GetNumberOfItems());
            i = this->MainSlicerWindow->GetEditMenu()->AddCommand ( "Edit Box", NULL, "::EditBox::ShowDialog" );
            this->MainSlicerWindow->GetEditMenu()->SetItemAccelerator ( i, "space");
            this->MainSlicerWindow->GetEditMenu()->SetBindingForItemAccelerator ( i, this->MainSlicerWindow);

#ifdef USE_PYTHON
            i = this->MainSlicerWindow->GetWindowMenu()->AddCommand ( "Python console", NULL, "$::slicer3::ApplicationGUI PythonConsole" );
            this->MainSlicerWindow->GetWindowMenu()->SetItemAccelerator ( i, "Ctrl+P");
            this->MainSlicerWindow->GetWindowMenu()->SetBindingForItemAccelerator ( i, this->MainSlicerWindow);
#endif

#ifndef FIDUCIALS_DEBUG
            this->GetMainSlicerWindow()->GetEditMenu()->InsertSeparator (this->GetMainSlicerWindow()->GetEditMenu()->GetNumberOfItems());
            // make the new fiducial list, but delete the returned node as
            // it's held onto by the scene
            i = this->MainSlicerWindow->GetEditMenu()->AddCommand ( "New Fiducial List", NULL, "[$::slicer3::FiducialsGUI GetLogic] AddFiducialListSelected" );
            this->MainSlicerWindow->GetEditMenu()->SetItemAccelerator ( i, "Ctrl+L");
            this->MainSlicerWindow->GetEditMenu()->SetBindingForItemAccelerator ( i, this->MainSlicerWindow);
#endif

            //
            // View Menu
            //

            this->MainSlicerWindow->GetViewMenu()->InsertCascade( 2, "Font size",  this->MainSlicerWindow->GetFontSizeMenu());
            int zz;
            const char *size;
            vtkSlicerFont *f = app->GetSlicerTheme()->GetSlicerFonts();
            int rindex;
            for ( zz = 0; zz < f->GetNumberOfFontSizes(); zz++ )
              {
              size = f->GetFontSize( zz );
              this->MainSlicerWindow->GetFontSizeMenu()->AddRadioButton ( size, NULL, "$::slicer3::ApplicationGUI SetApplicationFontSize");
              rindex = this->MainSlicerWindow->GetFontSizeMenu()->GetIndexOfItem ( size );
              this->MainSlicerWindow->GetFontSizeMenu()->SetItemSelectedValue ( rindex, size );
              //-- initialize interface to match application settings
              if ( ! (strcmp (app->GetApplicationFontSize(), size )))
                {
                this->MainSlicerWindow->GetFontSizeMenu()->SetItemSelectedState ( rindex, 1 );
                }
              else
                {
                this->MainSlicerWindow->GetFontSizeMenu()->SetItemSelectedState ( rindex, 0 );
                }
              }
            
            this->MainSlicerWindow->GetViewMenu()->InsertCascade( 2, "Font family",  this->MainSlicerWindow->GetFontFamilyMenu());
            const char *font;
            for ( zz = 0; zz < f->GetNumberOfFontFamilies(); zz++ )
              {
              font = f->GetFontFamily ( zz );
              this->MainSlicerWindow->GetFontFamilyMenu()->AddRadioButton ( font, NULL, "$::slicer3::ApplicationGUI SetApplicationFontFamily");
              rindex = this->MainSlicerWindow->GetFontFamilyMenu()->GetIndexOfItem ( font );
              this->MainSlicerWindow->GetFontFamilyMenu()->SetItemSelectedValue ( rindex, font );
              //-- initialize interface to match application settings
              if ( ! (strcmp (app->GetApplicationFontFamily( ), font )))
                {
                this->MainSlicerWindow->GetFontFamilyMenu()->SetItemSelectedState ( rindex, 1 );
                }
              else
                {
                this->MainSlicerWindow->GetFontFamilyMenu()->SetItemSelectedState ( rindex, 0 );
                }
              }

            this->GetMainSlicerWindow()->GetViewMenu()->InsertCommand (
                                                                       this->GetMainSlicerWindow()->GetViewMenuInsertPosition(),
                                                                       "Cache & Remote I/O Manager", NULL, "$::slicer3::RemoteIOGUI DisplayManagerWindow");

            //
            // Help Menu
            //
            this->GetMainSlicerWindow()->GetHelpMenu()->InsertCommand (
                                                                       this->GetMainSlicerWindow()->GetHelpMenuInsertPosition(),
                                                                       "Browse tutorials (not yet available)", NULL, "$::slicer3::ApplicationGUI OpenTutorialsLink");
            //
            // Feedback Menu
            //
            this->GetMainSlicerWindow()->GetFeedbackMenu()->AddCommand ("Feedback: report a bug (www)", NULL, "$::slicer3::ApplicationGUI OpenBugLink");
            this->GetMainSlicerWindow()->GetFeedbackMenu()->AddCommand ("Feedback: report usability issue (www)", NULL, "$::slicer3::ApplicationGUI OpenUsabilityLink");
            this->GetMainSlicerWindow()->GetFeedbackMenu()->AddCommand ("Feedback: make a feature request (www)", NULL, "$::slicer3::ApplicationGUI OpenFeatureLink");
            this->GetMainSlicerWindow()->GetFeedbackMenu()->AddCommand ("Community: Slicer Visual Blog (www)", NULL, "$::slicer3::ApplicationGUI PostToVisualBlog");            

            
            this->LoadSceneDialog->SetParent ( this->MainSlicerWindow );
            this->LoadSceneDialog->Create ( );
            this->LoadSceneDialog->SetFileTypes("{ {Scenes} {.mrml .xml} } { {MRML Scene} {.mrml} } { {Slicer2 Scene} {.xml} } { {Xcede Catalog} {.xcat} } { {All} {.*} }");
            this->LoadSceneDialog->RetrieveLastPathFromRegistry("OpenPath");

#endif
        }

    }
}


//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::SetApplicationFontFamily ( )
{
  vtkSlicerApplication *app = vtkSlicerApplication::SafeDownCast(this->GetApplication());
  if ( app )
    {
    vtkSlicerTheme *theme = app->GetSlicerTheme();
    if ( theme )
      {
      vtkSlicerFont *font = theme->GetSlicerFonts();
      if ( font )
        {
        vtkKWMenu *m = this->GetMainSlicerWindow()->GetFontFamilyMenu();

        const char *val;
        // what's the selected family?
        for ( int i=0; i < m->GetNumberOfItems(); i++ )
          {
          val = m->GetItemSelectedValue ( i );
          if ( m->GetItemSelectedState ( i ))
            {
            // i guess we found it check to see if m has a valid value:        
            if ( font->IsValidFontFamily ( val ) )
              {
              // save the change, and update the font
              theme->SetFontFamily ( val );
              app->SetApplicationFontFamily (val );
              app->Script ( "font configure %s -family %s", theme->GetApplicationFont2(), val);
              app->Script ( "font configure %s -family %s", theme->GetApplicationFont1(), val);
              app->Script ( "font configure %s -family %s", theme->GetApplicationFont0(), val);
              //--- and handle GUI widgets who aren't captured by the theme
              this->GetApplicationToolbar()->ReconfigureGUIFonts();
              this->GetViewControlGUI()->ReconfigureGUIFonts();
              this->GetMainSlicerWindow()->GetApplicationSettingsInterface()->Update();
              }
            break;
            }
          }
        
        }
      }
    }
}
    



//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::SetApplicationFontSize ( )
{
  vtkSlicerApplication *app = vtkSlicerApplication::SafeDownCast(this->GetApplication());
  if ( app )
    {
    vtkSlicerTheme *theme = app->GetSlicerTheme();
    if (theme)
      {
      vtkSlicerFont *font = theme->GetSlicerFonts();
      if ( font)
        {
        vtkKWMenu *m = this->GetMainSlicerWindow()->GetFontSizeMenu();
        
        const char *val;
        int f2, f1, f0;
        // what's the selected family?
        for ( int i=0; i < m->GetNumberOfItems(); i++ )
          {
          val = m->GetItemSelectedValue( i );
          if ( m->GetItemSelectedState ( i ))
            {
            // check to see if m has a valid value:
            if ( font->IsValidFontSize ( val) )
              {
              f2 = font->GetFontSize2( val );
              f1 = font->GetFontSize1( val );
              f0 = font->GetFontSize0( val );
              app->SetApplicationFontSize ( val );
              app->Script ( "font configure %s -size %d", theme->GetApplicationFont2(), f2);
              app->Script ( "font configure %s -size %d", theme->GetApplicationFont1(), f1);
              app->Script ( "font configure %s -size %d", theme->GetApplicationFont0(), f0);
              //--- and handle GUI widgets who aren't captured by the theme
              this->GetApplicationToolbar()->ReconfigureGUIFonts();
              this->GetViewControlGUI()->ReconfigureGUIFonts();
              this->GetMainSlicerWindow()->GetApplicationSettingsInterface()->Update();
              }
            }
          }
        }
      }
    }
}



//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::UpdateFontSizeMenu()
{
  vtkSlicerApplication *app = vtkSlicerApplication::SafeDownCast(this->GetApplication());
  if ( app )
    {
    // check to see if menu matches the Application Settings interface.
    // if so, do nothing.
    // if not, update.
    vtkKWMenu *m = this->GetMainSlicerWindow()->GetFontSizeMenu();
    if ( m )
      {
      const char *size = app->GetApplicationFontSize ( );
      if ( !( m->GetItemSelectedState ( size ) ))
        {
        m->SetItemSelectedState ( size, 1 );
        }
      }
    }
}


//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::UpdateFontFamilyMenu()
{
  vtkSlicerApplication *app = vtkSlicerApplication::SafeDownCast(this->GetApplication());
  if ( app )
    {
    // check to see if menu matches the Application Settings interface.
    // if so, do nothing.
    // if not, update.
    vtkKWMenu *m = this->GetMainSlicerWindow()->GetFontFamilyMenu();
    if ( m )
      {
      const char *font = app->GetApplicationFontFamily();
      if (!( m->GetItemSelectedState ( font ) ))
        {
        m->SetItemSelectedState ( font, 1 );
        }
      }    
    }
}



//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::InitializeSlicesControlGUI (  )
{
#ifndef SLICESCONTROL_DEBUG
  vtkSlicerSlicesControlGUI *scGUI = this->GetSlicesControlGUI ( );
  scGUI->UpdateFromMRML();
  scGUI->UpdateSliceGUIInteractorStyles ( );
#endif
}


//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::InitializeViewControlGUI (  )
{
#ifndef VIEWCONTROL_DEBUG
  vtkSlicerViewControlGUI *vcGUI = this->GetViewControlGUI ( );
  vcGUI->UpdateFromMRML();
  vcGUI->UpdateSliceGUIInteractorStyles();
  vcGUI->UpdateMainViewerInteractorStyles( );
  vcGUI->InitializeNavigationWidgetCamera( );
  vcGUI->UpdateNavigationWidgetViewActors ( );
  vcGUI->ConfigureNavigationWidgetRender ( );
#endif
}


//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::PythonConsole (  )
{
  
#ifdef USE_PYTHON
  PyObject* d = vtkSlicerApplication::GetInstance()->GetPythonDictionary();
  if ( d == NULL )
    {
    vtkSlicerApplication::GetInstance()->RequestDisplayMessage ( "Error", "Failed to startup python interpreter: dictionary null" );
    return;
    }
    
  PyObject* v = PyRun_StringFlags ( "import sys;\n"
                                    "try:\n"
                                    "  import Slicer;\n"
                                    "  reload ( Slicer );\n"
                                    "  Slicer.StartConsole();\n"
                                    "except Exception, e:\n"
                                    "  print 'Failed to import Slicer', e\n"
                                    "sys.stdout.flush();\n"
                                    "sys.stderr.flush();\n",
                                    Py_file_input,
                                    d,
                                    d,
                                    NULL);

  if (v == NULL)
    {
    PyErr_Print();
    vtkSlicerApplication::GetInstance()->RequestDisplayMessage ( "Error", "Failed to startup python interpreter" );
    return;
    }
  Py_DECREF ( v );
  if (Py_FlushLine())
    {
    PyErr_Clear();
    }
#endif
}


//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::SetCurrentModuleToHome (  )
{
#ifndef TOOLBAR_DEBUG
  if ( this->GetApplication() != NULL )
      {
      if ( this->GetApplicationToolbar()->GetModuleChooseGUI() )
        {
        if ( this->GetApplicationToolbar()->GetModuleChooseGUI()->GetModuleNavigator() )
          {
          vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();
          const char *name = this->GetApplicationToolbar()->GetModuleChooseGUI()->GetModuleNavigator()->GetCurrentModuleName ( );
          //--- save to registry.
          if (name != NULL)
            {
            app->SetHomeModule ( name );
            }
          this->GetMainSlicerWindow()->GetApplicationSettingsInterface()->Update();
          }
        }
      }
#endif
}


//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::OpenTutorialsLink ()
{
  if ( this->GetApplication() != NULL )
    {
    vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();
    //app->OpenLink ("" );
    }
}
//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::OpenBugLink ()
{
  if ( this->GetApplication() != NULL )
    {
    vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();
    app->OpenLink ("http://www.na-mic.org/Bug/index.php" );
    }
}
//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::OpenUsabilityLink ()
{
  if ( this->GetApplication() != NULL )
    {
    vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();
    app->OpenLink ("http://www.na-mic.org/Bug/index.php" );
    }

}
//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::OpenFeatureLink ()
{
  if ( this->GetApplication() != NULL )
    {
    vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();
    app->OpenLink ("http://www.na-mic.org/Bug/index.php" );
    }
}
//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::PostToVisualBlog ()
{
  if ( this->GetApplication() != NULL )
    {
    vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();
    app->OpenLink ("http://www.na-mic.org/Wiki/index.php/Slicer3:VisualBlog" );
    }
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::DestroyMainSliceViewers ( )
{

  if ( this->GetApplication() != NULL )
    {
      vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();
      vtkSlicerGUILayout *layout = app->GetMainLayout ( );

      //
      // Destroy 3 main slice viewers
      //
          if (this->SlicesGUI)
          {
                  vtkSlicerSliceGUI *g = NULL;
                  char *layoutname = NULL;
                  int nSliceGUI = this->SlicesGUI->GetNumberOfSliceGUI();
                  for (int i = 0; i < nSliceGUI; i++)
                  {
                          if (i == 0)
                          {
                                  g = this->SlicesGUI->GetFirstSliceGUI();
                                  layoutname = this->SlicesGUI->GetFirstSliceGUILayoutName();
                          }
                          else
                          {
                                  g = this->SlicesGUI->GetNextSliceGUI(layoutname);
                                  layoutname = this->SlicesGUI->GetNextSliceGUILayoutName(layoutname);
                          }
                          g->SetAndObserveMRMLScene(NULL);
                          g->SetAndObserveModuleLogic(NULL);
                          g->RemoveGUIObservers();
                          g->SetApplicationLogic(NULL);

                          if (layout->GetCurrentViewArrangement() == vtkSlicerGUILayout::SlicerLayoutFourUpView)
                                  g->UngridGUI();
                          else
                                  g->UnpackGUI();
                          g->Delete();
                          g = NULL;
                  }
          }

      if ( this->GridFrame1 )
        {
          app->Script ("pack forget %s ", this->GridFrame1->GetWidgetName ( ) );
          this->GridFrame1->Delete ( );
          this->GridFrame1 = NULL;
        }
      if ( this->GridFrame2 )
        {
          app->Script ("pack forget %s ", this->GridFrame2->GetWidgetName ( ) );
          this->GridFrame2->Delete ( );
          this->GridFrame2 = NULL;
        }
    }
}


//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::DestroyMain3DViewer ( )
{
  //

  if ( this->GetApplication() != NULL )
    {
      vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();
      vtkSlicerGUILayout *layout = app->GetMainLayout ( );

      // Destroy fiducial list
      if ( this->FiducialListWidget )
      {
          this->FiducialListWidget->RemoveMRMLObservers ();
          this->FiducialListWidget->SetParent(NULL);
          this->FiducialListWidget->Delete();
          this->FiducialListWidget = NULL;
      }
      
      // Destroy main 3D viewer
      //
      if ( this->ViewerWidget )
        {
        this->ViewerWidget->RemoveMRMLObservers ( );
        if ( layout->GetCurrentViewArrangement() == vtkSlicerGUILayout::SlicerLayoutFourUpView )
          {
            this->ViewerWidget->UngridWidget ( );
          }
        else
          {
            this->ViewerWidget->UnpackWidget ( );
          }
        this->ViewerWidget->SetApplicationLogic ( NULL );
        this->ViewerWidget->SetParent ( NULL );
        this->ViewerWidget->Delete ( );
        this->ViewerWidget = NULL;
        }
    }
}



//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::DisplayMainSlicerWindow ( )
{
  if ( this->GetApplication() != NULL )
    {
      vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();
      
      this->MainSlicerWindow->Display ( );
      int w = this->MainSlicerWindow->GetWidth ( );
      int h = this->MainSlicerWindow->GetHeight ( );
      int vh = app->GetMainLayout()->GetDefault3DViewerHeight();
      int sh = app->GetMainLayout()->GetDefaultSliceGUIFrameHeight();
      int sfh = this->MainSlicerWindow->GetSecondarySplitFrame()->GetFrame1Size();
      int sf2h = this->MainSlicerWindow->GetSecondarySplitFrame()->GetFrame2Size();
      
      // pop up a warning dialog here if the computer's
      // display resolution in x is less than 1000 pixels.
      /*
      const char *wstr = app->Script ("winfo screenwidth .");
      const char *hstr = app->Script ("winfo screenheight .");
      int screenwidth = atoi (wstr);
      int screenheight = atoi ( hstr );
      if (screenwidth < 1000 )
        {
        vtkKWMessageDialog *message = vtkKWMessageDialog::New();
        message->SetParent ( this->MainSlicerWindow );
        message->SetStyleToMessage();
        message->SetDialogName("WarningScreenResolution");
        message->SetText ("Slicer requires a horizontal screen resolution of at least 1024 pixels to display it's user interface. Some GUI elements may not be visible.");
        message->Create();
        message->Invoke();
        message->Delete();
        }
      */
    }
}

    

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::BuildMainViewer ( int arrangementType)
{

  if ( this->GetApplication() != NULL )
    {
      vtkSlicerApplication *app = vtkSlicerApplication::SafeDownCast(this->GetApplication());
      vtkSlicerColor *color = app->GetSlicerTheme()->GetSlicerColors ( );
      vtkSlicerGUILayout *layout = app->GetMainLayout ( );
      vtkSlicerWindow *win = this->MainSlicerWindow;
        
      this->GridFrame1->SetParent ( this->MainSlicerWindow->GetViewFrame ( ) );
      this->GridFrame1->Create ( );            
      this->GridFrame2->SetParent ( this->MainSlicerWindow->GetSecondaryPanelFrame ( ) );
      this->GridFrame2->Create ( );            
      this->CreateMainSliceViewers ( arrangementType );
      this->CreateMain3DViewer (arrangementType );
      this->PackMainViewer ( arrangementType , NULL );
    }
}

vtkSlicerSliceGUI* vtkSlicerApplicationGUI::GetMainSliceGUI(char *layoutName)
{
  if (this->SlicesGUI)
    return this->SlicesGUI->GetSliceGUI(layoutName);
  else
    return NULL;
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::CreateMainSliceViewers ( int arrangementType )
{
        if ( this->GetApplication() != NULL )
        {
          vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();
          vtkSlicerColor *color = app->GetSlicerTheme()->GetSlicerColors ( );
          
          //
          // 3 Slice Viewers
          //
          vtkSlicerSliceGUI *g = vtkSlicerSliceGUI::New ( );
          this->SlicesGUI->AddSliceGUI("Red", g);
          g = vtkSlicerSliceGUI::New ( );
          this->SlicesGUI->AddSliceGUI("Yellow", g);
          g = vtkSlicerSliceGUI::New ( );
          this->SlicesGUI->AddSliceGUI("Green", g);
          if (this->SlicesGUI)
            {
              char *layoutname = NULL;
              int nSliceGUI = this->SlicesGUI->GetNumberOfSliceGUI();
              for (int i = 0; i < nSliceGUI; i++)
                {
                  if (i == 0)
                    {
                      g = this->SlicesGUI->GetFirstSliceGUI();
                      layoutname = this->SlicesGUI->GetFirstSliceGUILayoutName();
                    }
                  else
                    {
                      g = this->SlicesGUI->GetNextSliceGUI(layoutname);
                      layoutname = this->SlicesGUI->GetNextSliceGUILayoutName(layoutname);
                    }
                  g->SetApplication(app);
                  g->SetApplicationLogic(this->ApplicationLogic);
                  
                  switch (layoutname[0])
                    {
                    case 'R':
                      //if ( this->MainSliceLogic0 )
                      //      g->SetAndObserveModuleLogic( this->MainSliceLogic0 );
                      if ( this->GetApplicationLogic()->GetSlicerSliceLogic("Red") )
                        g->SetAndObserveModuleLogic( this->GetApplicationLogic()->GetSlicerSliceLogic("Red") );
                      
                      g->BuildGUI( this->MainSlicerWindow->GetMainSplitFrame( ), color->SliceGUIRed );
                      break;
                    case 'Y':
                      if ( this->GetApplicationLogic()->GetSlicerSliceLogic("Yellow") )
                        g->SetAndObserveModuleLogic( this->GetApplicationLogic()->GetSlicerSliceLogic("Yellow") );
                      
                      g->BuildGUI( this->MainSlicerWindow->GetMainSplitFrame( ), color->SliceGUIYellow );
                      break;
                    case 'G':
                      if ( this->GetApplicationLogic()->GetSlicerSliceLogic("Green") )
                        g->SetAndObserveModuleLogic( this->GetApplicationLogic()->GetSlicerSliceLogic("Green") );
                      
                      g->BuildGUI( this->MainSlicerWindow->GetMainSplitFrame( ), color->SliceGUIGreen );
                      break;
                    }
                  g->AddGUIObservers();
                  g->SetAndObserveMRMLScene( this->MRMLScene );
                }
            }
        }
}


//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::CreateMain3DViewer ( int arrangementType )
{
  if ( this->GetApplication() != NULL )
    {
      vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();
      vtkSlicerColor *color = app->GetSlicerTheme()->GetSlicerColors ( );
      vtkSlicerGUILayout *layout = app->GetMainLayout ( );
      //
      // Make 3D Viewer
      //
      this->ViewerWidget = vtkSlicerViewerWidget::New ( );
      this->ViewerWidget->SetApplication( app );
      this->ViewerWidget->SetParent(this->MainSlicerWindow->GetViewPanelFrame());
      
      // add events
      vtkIntArray *events = vtkIntArray::New();
      events->InsertNextValue(vtkMRMLScene::SceneCloseEvent);
      events->InsertNextValue(vtkMRMLScene::NewSceneEvent);
      events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
      events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
      events->InsertNextValue(vtkCommand::ModifiedEvent);
      this->ViewerWidget->SetAndObserveMRMLSceneEvents (this->MRMLScene, events );
      
      // use the events for the fiducial list widget as well
      //events->Delete();
      this->ViewerWidget->Create();
      this->ViewerWidget->GetMainViewer()->SetRendererBackgroundColor (app->GetSlicerTheme()->GetSlicerColors()->ViewerBlue );
      this->ViewerWidget->UpdateFromMRML();
      this->ViewerWidget->SetApplicationLogic ( this->GetApplicationLogic () );
      // add the fiducial list widget
      this->FiducialListWidget = vtkSlicerFiducialListWidget::New();
      this->FiducialListWidget->SetApplication( app );
      this->FiducialListWidget->SetMainViewer(this->ViewerWidget->GetMainViewer());
      this->FiducialListWidget->SetViewerWidget(this->ViewerWidget);
      this->FiducialListWidget->SetInteractorStyle(vtkSlicerViewerInteractorStyle::SafeDownCast(this->ViewerWidget->GetMainViewer()->GetRenderWindowInteractor()->GetInteractorStyle()));
      this->FiducialListWidget->Create();
      this->FiducialListWidget->SetAndObserveMRMLSceneEvents (this->MRMLScene, events );
      events->Delete();
      this->FiducialListWidget->UpdateFromMRML();
    }
}



//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::PackMainViewer ( int arrangmentType, const char *whichSlice)
{
  if ( this->GetApplication() != NULL )
    {
      vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();
      // parent the sliceGUI  based on selected view arrangement & build
      
      switch ( arrangmentType)
        {
        case vtkSlicerGUILayout::SlicerLayoutInitialView:
          this->PackConventionalView ( );
          break;
        case vtkSlicerGUILayout::SlicerLayoutDefaultView:
          this->PackConventionalView ( );
          break;
        case vtkSlicerGUILayout::SlicerLayoutFourUpView:
          this->PackFourUpView ( );
          break;
        case vtkSlicerGUILayout::SlicerLayoutOneUp3DView:
          this->PackOneUp3DView ( );
          break;
        case vtkSlicerGUILayout::SlicerLayoutOneUpRedSliceView:
          this->PackOneUpSliceView ("Red");
          break;
        case vtkSlicerGUILayout::SlicerLayoutOneUpYellowSliceView:
          this->PackOneUpSliceView ("Yellow");
          break;
        case vtkSlicerGUILayout::SlicerLayoutOneUpGreenSliceView:
          this->PackOneUpSliceView ("Green");
          break;
        case vtkSlicerGUILayout::SlicerLayoutOneUpSliceView:
          this->PackOneUpSliceView ( whichSlice );
          break;
        case vtkSlicerGUILayout::SlicerLayoutTabbed3DView:
          this->PackTabbed3DView ( );
          break;
        case vtkSlicerGUILayout::SlicerLayoutTabbedSliceView:
          this->PackTabbedSliceView ( );
          break;
        case vtkSlicerGUILayout::SlicerLayoutCompareView:
          this->PackCompareView();
          //this->CompareViewSelect();
          break;
        case vtkSlicerGUILayout::SlicerLayoutLightboxView:
          this->PackLightboxView ( );
          break;
        default:
          this->PackConventionalView ( );
          break;
        }
    }
}


//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::RepackMainViewer ( int arrangementType, const char *whichSlice )
{
  this->UnpackMainSliceViewers ( );
  this->UnpackMain3DViewer ( );
  // Since I can't find a way to re-title this main page titled "View",
  // we make sure it's visible, and then 'hide' it only when we want to
  // show tabs that say things other than "View".
  this->MainSlicerWindow->GetViewNotebook()->ShowPage ( "View");
  this->PackMainViewer ( arrangementType, whichSlice );
}




//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::UnpackMainSliceViewers (  )
{
  if ( this->GetApplication() != NULL )
    {
      vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();
      vtkSlicerGUILayout *layout = app->GetMainLayout ( );
      //
      // Unpack or Ungrid main slice viewers
      //
      if (this->SlicesGUI)
        {
          vtkSlicerSliceGUI *g = NULL;
          char *layoutname = NULL;
          int nSliceGUI = this->SlicesGUI->GetNumberOfSliceGUI();
          for (int i = 0; i < nSliceGUI; i++)
            {
              if (i == 0)
                {
                  g = this->SlicesGUI->GetFirstSliceGUI();
                  layoutname = this->SlicesGUI->GetFirstSliceGUILayoutName();
                }
              else
                {
                  g = this->SlicesGUI->GetNextSliceGUI(layoutname);
                  layoutname = this->SlicesGUI->GetNextSliceGUILayoutName(layoutname);
                }
              
              if ( layout->GetCurrentViewArrangement() == vtkSlicerGUILayout::SlicerLayoutFourUpView ||
                   layout->GetCurrentViewArrangement() == vtkSlicerGUILayout::SlicerLayoutDefaultView)
                {
                  g->UngridGUI();
                }
              else
                {
                  g->UnpackGUI();
                }
            }
        }
      
      if ( this->GridFrame1 )
        {
          app->Script ("pack forget %s ", this->GridFrame1->GetWidgetName ( ) );
        }
    }
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::UnpackMain3DViewer (  )
{
  if ( this->GetApplication() != NULL )
    {
    vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();
    vtkSlicerColor *color = app->GetSlicerTheme()->GetSlicerColors ( );
    vtkSlicerGUILayout *layout = app->GetMainLayout ( );
    //
    // 3D Viewer
    //
    this->MainSlicerWindow->GetViewNotebook()->RemovePagesMatchingTag ( this->ViewerPageTag );
    if ( layout->GetCurrentViewArrangement() == vtkSlicerGUILayout::SlicerLayoutFourUpView )
      {
      this->ViewerWidget->UngridWidget ( );
      }
    else
      {
      this->ViewerWidget->UnpackWidget ( );
      }
    }
}


//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::PackConventionalView ( )
{
  if ( this->GetApplication() != NULL )
    {
      vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();
      vtkSlicerColor *color = app->GetSlicerTheme()->GetSlicerColors ( );
      vtkSlicerGUILayout *layout = app->GetMainLayout ( );

      // Expose the main panel frame and secondary panel frame.

//      this->MainSlicerWindow->SetMainPanelVisibility ( 1 );
      this->MainSlicerWindow->SetSecondaryPanelVisibility ( 1 );

        // try gridding the sliceGUIs inside a "GridFrame2"
        // which itself is packed in the SecondaryPanelFrame.
        // use column configure to make each GUI resize evenly.
      this->Script ( "pack %s -side top -fill both -expand y -padx 0 -pady 0 ", this->GridFrame2->GetWidgetName ( ) );
      this->Script ("grid columnconfigure %s 0 -weight 1", this->GridFrame2->GetWidgetName() );
      this->Script ("grid columnconfigure %s 1 -weight 1", this->GridFrame2->GetWidgetName() );
      this->Script ("grid columnconfigure %s 2 -weight 1", this->GridFrame2->GetWidgetName() );
      
      this->ViewerWidget->PackWidget(this->MainSlicerWindow->GetViewFrame() );
      vtkSlicerSliceGUI *g = this->SlicesGUI->GetSliceGUI("Red");
      
      int nSlicesGUIInFrame2 = this->SlicesGUI->GetNumberOfSliceGUI() - 1;
      
      if (layout->GetCurrentViewArrangement() == vtkSlicerGUILayout::SlicerLayoutCompareView)
        //|| layout->GetCurrentViewArrangement() == vtkSlicerGUILayout::SlicerLayoutFourUpView )
        {
        g->UngridGUI();
        this->Script ("grid rowconfigure %s 0 -weight 1", this->GridFrame2->GetWidgetName() );
        
        g->GridSpanGUI ( this->GetGridFrame2(), 0, 0, nSlicesGUIInFrame2, 1 );
        g = this->SlicesGUI->GetSliceGUI("Yellow");
        g->GridSpanGUI( this->GetGridFrame2(), 0, 1, nSlicesGUIInFrame2, 1 );
        g = this->SlicesGUI->GetSliceGUI("Green");
        g->GridSpanGUI ( this->GetGridFrame2(), 0, 2, nSlicesGUIInFrame2, 1 );
        }
      else
        {
        this->Script ("grid rowconfigure %s 0 -weight 1", this->GridFrame2->GetWidgetName() );
        g->GridGUI ( this->GetGridFrame2(), 0, 0 );
        g = this->SlicesGUI->GetSliceGUI("Yellow");
        g->GridGUI ( this->GetGridFrame2(), 0, 1 );
        g = this->SlicesGUI->GetSliceGUI("Green");
        g->GridGUI ( this->GetGridFrame2(), 0, 2 );
        }
      
      this->MainSlicerWindow->GetViewNotebook()->SetAlwaysShowTabs ( 0 );
      layout->SetCurrentViewArrangement ( vtkSlicerGUILayout::SlicerLayoutDefaultView );
      if ( layout->GetDefaultSliceGUIFrameHeight() > 0 )
        {
        this->MainSlicerWindow->GetSecondarySplitFrame()->SetFrame1Size ( layout->GetDefaultSliceGUIFrameHeight () );
        }
      else
        {
        this->MainSlicerWindow->GetSecondarySplitFrame()->SetFrame1Size ( layout->GetStandardSliceGUIFrameHeight () );
        }
    }
}


//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::PackOneUp3DView ( )
{
  if ( this->GetApplication() != NULL )
    {
      vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();
      vtkSlicerColor *color = app->GetSlicerTheme()->GetSlicerColors ( );
      vtkSlicerGUILayout *layout = app->GetMainLayout ( );

      // Expose the main panel frame only
//      this->MainSlicerWindow->SetMainPanelVisibility ( 1 );
      this->MainSlicerWindow->SetSecondaryPanelVisibility ( 0 );
      this->MainSlicerWindow->GetViewNotebook()->SetAlwaysShowTabs ( 0 );      
      this->ViewerWidget->PackWidget(this->MainSlicerWindow->GetViewFrame() );
      
      vtkSlicerSliceGUI *g = this->SlicesGUI->GetSliceGUI("Red");
      g->PackGUI(this->MainSlicerWindow->GetSecondaryPanelFrame());
      g = this->SlicesGUI->GetSliceGUI("Yellow");
      g->PackGUI(this->MainSlicerWindow->GetSecondaryPanelFrame());
      g = this->SlicesGUI->GetSliceGUI("Green");
      g->PackGUI(this->MainSlicerWindow->GetSecondaryPanelFrame());
      
      this->MainSlicerWindow->GetViewNotebook()->SetAlwaysShowTabs ( 0 );
      layout->SetCurrentViewArrangement ( vtkSlicerGUILayout::SlicerLayoutOneUp3DView );
    }
}


//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::PackOneUpSliceView ( const char * whichSlice )
{
  if ( this->GetApplication() != NULL )
    {
      vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();
      vtkSlicerColor *color = app->GetSlicerTheme()->GetSlicerColors ( );
      vtkSlicerGUILayout *layout = app->GetMainLayout ( );

      // Expose the main panel frame only
//      this->MainSlicerWindow->SetMainPanelVisibility ( 1 );
      this->MainSlicerWindow->SetSecondaryPanelVisibility ( 0 );

      vtkSlicerSliceGUI *g;

      if ( !strcmp (whichSlice, "Red" ) )
        {
          g = this->SlicesGUI->GetSliceGUI("Red");
          g->PackGUI( this->MainSlicerWindow->GetViewFrame());
          g = this->SlicesGUI->GetSliceGUI("Yellow");
          g->PackGUI(NULL);
          g = this->SlicesGUI->GetSliceGUI("Green");
          g->PackGUI(NULL);
          layout->SetCurrentViewArrangement ( vtkSlicerGUILayout::SlicerLayoutOneUpRedSliceView );
        }
      else if ( !strcmp ( whichSlice, "Yellow" ) )
        {
          g = this->SlicesGUI->GetSliceGUI("Red");
          g->PackGUI(NULL);
          g = this->SlicesGUI->GetSliceGUI("Yellow");
          g->PackGUI( this->MainSlicerWindow->GetViewFrame());
          g = this->SlicesGUI->GetSliceGUI("Green");
          g->PackGUI(NULL);
          layout->SetCurrentViewArrangement ( vtkSlicerGUILayout::SlicerLayoutOneUpYellowSliceView );
        }
      else if ( !strcmp ( whichSlice, "Green" ) )
        {
          g = this->SlicesGUI->GetSliceGUI("Red");
          g->PackGUI(NULL);
          g = this->SlicesGUI->GetSliceGUI("Yellow");
          g->PackGUI(NULL);
          g = this->SlicesGUI->GetSliceGUI("Green");
          g->PackGUI( this->MainSlicerWindow->GetViewFrame());
          layout->SetCurrentViewArrangement ( vtkSlicerGUILayout::SlicerLayoutOneUpGreenSliceView );
        }
      
      //      this->ViewerWidget->PackWidget(this->MainSlicerWindow->GetViewFrame() );
      this->MainSlicerWindow->GetViewNotebook()->SetAlwaysShowTabs ( 0 );
    }
}


//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::PackFourUpView ( )
{
  if ( this->GetApplication() != NULL )
    {
      vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();
      vtkSlicerColor *color = app->GetSlicerTheme()->GetSlicerColors ( );
      vtkSlicerGUILayout *layout = app->GetMainLayout ( );
      
      // Expose both the main panel frame and secondary panel frame
      //      this->MainSlicerWindow->SetMainPanelVisibility ( 1 );
      this->MainSlicerWindow->SetSecondaryPanelVisibility ( 0 );
      this->MainSlicerWindow->GetSecondarySplitFrame()->SetFrame1Size ( 0 );
      
      // Use this frame in MainSlicerWindow's ViewFrame to grid in the various viewers.
      this->Script ( "pack %s -side top -fill both -expand y -padx 0 -pady 0 ", this->GridFrame1->GetWidgetName ( ) );
      this->Script ("grid rowconfigure %s 0 -weight 1", this->GridFrame1->GetWidgetName() );
      this->Script ("grid rowconfigure %s 1 -weight 1", this->GridFrame1->GetWidgetName() );
      this->Script ("grid columnconfigure %s 0 -weight 1", this->GridFrame1->GetWidgetName() );
      this->Script ("grid columnconfigure %s 1 -weight 1", this->GridFrame1->GetWidgetName() );

      this->ViewerWidget->GridWidget ( this->GridFrame1, 0, 1 );
      vtkSlicerSliceGUI *g = this->SlicesGUI->GetSliceGUI("Red");
      g->GridGUI( this->GetGridFrame1(), 0, 0 );
      g = this->SlicesGUI->GetSliceGUI("Yellow");
      g->GridGUI( this->GetGridFrame1(), 1, 0 );
      g = this->SlicesGUI->GetSliceGUI("Green");
      g->GridGUI( this->GetGridFrame1(), 1, 1 );
      
      this->MainSlicerWindow->GetViewNotebook()->SetAlwaysShowTabs ( 0 );
      layout->SetCurrentViewArrangement ( vtkSlicerGUILayout::SlicerLayoutFourUpView );
    }
}




//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::PackTabbed3DView ( )
{

  // TODO: implement multi-tabbed ViewerWidgets
  if ( this->GetApplication() != NULL )
    {
      vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();
      vtkSlicerColor *color = app->GetSlicerTheme()->GetSlicerColors ( );
      vtkSlicerGUILayout *layout = app->GetMainLayout ( );
      if ( layout->GetDefaultSliceGUIFrameHeight() > 0 )
        {
        this->MainSlicerWindow->GetSecondarySplitFrame()->SetFrame1Size ( layout->GetDefaultSliceGUIFrameHeight() );
        }
      else
        {
        this->MainSlicerWindow->GetSecondarySplitFrame()->SetFrame1Size ( layout->GetStandardSliceGUIFrameHeight() );
        }
//      this->MainSlicerWindow->SetMainPanelVisibility ( 1 );
      this->MainSlicerWindow->SetSecondaryPanelVisibility ( 0 );

      vtkSlicerSliceGUI *g = this->SlicesGUI->GetSliceGUI("Red");
      g->PackGUI( this->MainSlicerWindow->GetSecondaryPanelFrame( ));
          g = this->SlicesGUI->GetSliceGUI("Yellow");
      g->PackGUI( this->MainSlicerWindow->GetSecondaryPanelFrame( ));
          g = this->SlicesGUI->GetSliceGUI("Green");
      g->PackGUI( this->MainSlicerWindow->GetSecondaryPanelFrame( ));

      // Add a page for the current view, and each saved view.
      this->MainSlicerWindow->GetViewNotebook()->AddPage("Current view", NULL, NULL, this->ViewerPageTag );
      this->ViewerWidget->PackWidget(this->MainSlicerWindow->GetViewNotebook()->GetFrame ("Current view" ));

      // don't know how to change the title of this one,
      // so just hide it in this configuration, and expose
      // it again when the view configuration changes.
      this->MainSlicerWindow->GetViewNotebook()->HidePage ( "View");
      
      // Tab the 3D view
      this->MainSlicerWindow->GetViewNotebook()->SetAlwaysShowTabs ( 1 );
      layout->SetCurrentViewArrangement ( vtkSlicerGUILayout::SlicerLayoutTabbed3DView );
    }

}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::PackTabbedSliceView ( )
{
  // TODO: implement this and add an icon on the toolbar for it
  if ( this->GetApplication() != NULL )
    {
      vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();
      vtkSlicerColor *color = app->GetSlicerTheme()->GetSlicerColors ( );
      vtkSlicerGUILayout *layout = app->GetMainLayout ( );
      //      this->MainSlicerWindow->SetMainPanelVisibility ( 1 );
      this->MainSlicerWindow->SetSecondaryPanelVisibility ( 0 );
      
      this->MainSlicerWindow->GetViewNotebook()->AddPage("Red slice", NULL, NULL, this->ViewerPageTag );
      
      vtkSlicerSliceGUI *g = this->SlicesGUI->GetSliceGUI("Red");
      
      g->PackGUI ( this->MainSlicerWindow->GetViewNotebook( )->GetFrame ("Red slice") );
      this->MainSlicerWindow->GetViewNotebook()->AddPage("Yellow slice", NULL, NULL, this->ViewerPageTag );
      
      g = this->SlicesGUI->GetSliceGUI("Yellow");
      g->PackGUI ( this->MainSlicerWindow->GetViewNotebook( )->GetFrame ("Yellow slice") );
      this->MainSlicerWindow->GetViewNotebook()->AddPage("Green slice", NULL, NULL, this->ViewerPageTag );
      
      g = this->SlicesGUI->GetSliceGUI("Green");
      g->PackGUI ( this->MainSlicerWindow->GetViewNotebook( )->GetFrame ("Green slice") );
      
      // don't know how to change the title of this one,
      // so just hide it in this configuration, and expose
      // it again when the view configuration changes.
      this->MainSlicerWindow->GetViewNotebook()->HidePage ( "View");
      //      this->ViewerWidget->PackWidget(this->MainSlicerWindow->GetViewFrame() );
      
      // Tab the slices
      this->MainSlicerWindow->GetViewNotebook()->SetAlwaysShowTabs ( 1 );
      layout->SetCurrentViewArrangement ( vtkSlicerGUILayout::SlicerLayoutTabbedSliceView );
    }
}

void vtkSlicerApplicationGUI::PackCompareView()
{
  if ( this->GetApplication() != NULL )
    {
    vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();
    vtkSlicerColor *color = app->GetSlicerTheme()->GetSlicerColors ( );
    vtkSlicerGUILayout *layout = app->GetMainLayout ( );     
    
    this->MainSlicerWindow->SetSecondaryPanelVisibility ( 1 );
    
    // setup the layout for Frame1
    this->Script ( "pack %s -side top -fill both -expand y -padx 0 -pady 0 ", this->GridFrame1->GetWidgetName ( ) );
    this->Script ("grid rowconfigure %s 0 -weight 1", this->GridFrame1->GetWidgetName() );
    this->Script ("grid columnconfigure %s 0 -weight 1", this->GridFrame1->GetWidgetName() );
    this->Script ("grid columnconfigure %s 1 -weight 1", this->GridFrame1->GetWidgetName() );
    
    vtkSlicerSliceGUI *g = this->SlicesGUI->GetSliceGUI("Red");
    
    if (layout->GetCurrentViewArrangement() == vtkSlicerGUILayout::SlicerLayoutFourUpView || 
        layout->GetCurrentViewArrangement() == vtkSlicerGUILayout::SlicerLayoutDefaultView)
      {
      this->ViewerWidget->GridSpanWidget ( this->GridFrame1, 0, 1, 2, 1 );
      g->GridSpanGUI ( this->GetGridFrame1( ), 0, 0, 2, 1 );
      }
    else
      {
      this->ViewerWidget->GridWidget ( this->GridFrame1, 0, 1 );
      g->GridGUI ( this->GetGridFrame1( ), 0, 0 );
      }
    g->GetSliceViewer()->SetWidth(layout->GetDefaultSliceGUIFrameWidth());
    
      
      
    // insert a number of new main slice viewers according to user's input
    this->AddMainSliceGUI("Compare0");
    this->AddMainSliceGUI("Compare1");
//    this->AddMainSliceGUI("Compare2");
//    this->AddMainSliceGUI("Compare3");
    
    // configure the new layout
    this->Script ( "pack %s -side top -fill both -expand y -padx 0 -pady 0 ", this->GridFrame2->GetWidgetName ( ) );
    this->Script ("grid columnconfigure %s 0 -weight 1", this->GridFrame2->GetWidgetName() );
    
    char *layoutname = NULL;
    int nSliceGUI = this->SlicesGUI->GetNumberOfSliceGUI();
    int ncount = 0;
    for (int i = 0; i < nSliceGUI; i++)
      {
      if (i == 0)
        {
        g = this->SlicesGUI->GetFirstSliceGUI();
        layoutname = this->SlicesGUI->GetFirstSliceGUILayoutName();
        }
      else
        {
        g = this->SlicesGUI->GetNextSliceGUI(layoutname);
        layoutname = this->SlicesGUI->GetNextSliceGUILayoutName(layoutname);
        }
          
      if ( strcmp(layoutname, "Red") == 0 ||
           strcmp(layoutname, "Yellow") == 0 ||
           strcmp(layoutname, "Green") == 0)
        continue;
      else
        {
        g->GridSpanGUI( this->GetGridFrame2( ), ncount, 0, 1, 3 );
        g->GetSliceViewer()->SetWidth(layout->GetDefaultSliceGUIFrameWidth());
        this->Script ("grid rowconfigure %s %d -weight 1", this->GridFrame2->GetWidgetName(), ncount );
        ncount++;
        }
      }
      
    this->MainSlicerWindow->GetViewNotebook()->SetAlwaysShowTabs ( 0 );
    layout->SetCurrentViewArrangement( vtkSlicerGUILayout::SlicerLayoutCompareView );
    if ( layout->GetDefaultSliceGUIFrameHeight() > 0 )
      {
      this->MainSlicerWindow->GetSecondarySplitFrame()->SetFrame1Size ( layout->GetDefaultSliceGUIFrameHeight () );
      }
    else
      {
      this->MainSlicerWindow->GetSecondarySplitFrame()->SetFrame1Size ( layout->GetStandardSliceGUIFrameHeight () );
      }
    }
}

void vtkSlicerApplicationGUI::AddMainSliceGUI(char *layoutName)
{
  if ( this->GetApplication() != NULL )
    {
      vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();
      vtkSlicerColor *color = app->GetSlicerTheme()->GetSlicerColors ( );
      // if get "compare0" is NULL, add it
      vtkSlicerSliceLogic *sliceLogic = this->GetApplicationLogic()->GetSlicerSliceLogic(layoutName);
      if (sliceLogic == NULL)
        {
          vtkIntArray *events = vtkIntArray::New();
          events->InsertNextValue(vtkMRMLScene::NewSceneEvent);
          events->InsertNextValue(vtkMRMLScene::SceneCloseEvent);
          events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
          events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
          
          sliceLogic = vtkSlicerSliceLogic::New ( );
          this->GetApplicationLogic()->AddSliceLogic(layoutName, sliceLogic);
          sliceLogic->SetName(layoutName);
          
          sliceLogic->SetMRMLScene( this->MRMLScene );
          sliceLogic->ProcessLogicEvents();
          sliceLogic->ProcessMRMLEvents(this->MRMLScene, vtkCommand::ModifiedEvent, NULL);
          sliceLogic->SetAndObserveMRMLSceneEvents( this->GetMRMLScene(), events );
          
          events->Delete();
        }
      
      vtkSlicerSliceGUI *g = this->SlicesGUI->GetSliceGUI(layoutName);
      if (g == NULL)
        {
          g = vtkSlicerSliceGUI::New( );
          this->SlicesGUI->AddSliceGUI(layoutName, g);
          g->SetApplication(app);
          g->SetApplicationLogic(this->ApplicationLogic);
          g->BuildGUI ( this->MainSlicerWindow->GetMainSplitFrame(), color->SliceGUIYellow );
          //g->GetSliceController()->GetSliceNode()->SetOrientationToSagittal();
          g->AddGUIObservers();
          g->SetAndObserveMRMLScene(this->MRMLScene);
          g->SetAndObserveModuleLogic(sliceLogic);

          // need to check the delete section to prevent memory leak
          g->GetSliceController()->AddObserver(vtkSlicerSliceControllerWidget::ExpandEvent, 
                                      (vtkCommand *)this->GUICallbackCommand);
          g->GetSliceController()->AddObserver(vtkSlicerSliceControllerWidget::ShrinkEvent, 
                                      (vtkCommand *)this->GUICallbackCommand);
        }
    }
}

void vtkSlicerApplicationGUI::CompareViewSelect()
{
    //--- create temporary pop-up frame to display when custom configuration is selected
    vtkKWFrame *popUpFrame1 = vtkKWFrame::New ( );
    popUpFrame1->SetParent(this->TopFrame);
    popUpFrame1->Create ( );
    //popUpFrame1->SetBinding ( "<Leave>", this, "HideLightboxCustomLayoutFrame" );
    this->Script ( "pack %s -side left -anchor w -padx 2 -pady 2 -fill x -fill y -expand n", popUpFrame1->GetWidgetName ( ) );   
    this->CompareViewBoxRowEntry = vtkKWEntry::New ( );
    this->CompareViewBoxRowEntry->SetParent ( popUpFrame1 );
    this->CompareViewBoxRowEntry->Create ( );
    this->CompareViewBoxRowEntry->SetValueAsInt (1);
    this->CompareViewBoxRowEntry->SetWidth ( 3 );
    this->CompareViewBoxColumnEntry = vtkKWEntry::New ( );
    this->CompareViewBoxColumnEntry->SetParent ( popUpFrame1 );
    this->CompareViewBoxColumnEntry->Create ( );
    this->CompareViewBoxColumnEntry->SetWidth ( 3 );
    this->CompareViewBoxColumnEntry->SetValueAsInt (1);
    vtkKWLabel *rowsLabel = vtkKWLabel::New();
    rowsLabel->SetParent ( popUpFrame1 );
    rowsLabel->Create ( );
    rowsLabel->SetText ( "Number of rows:" );
    vtkKWLabel *columnsLabel = vtkKWLabel::New();
    columnsLabel->SetParent ( popUpFrame1 );
    columnsLabel->Create ( );
    columnsLabel->SetText ( "Number of columns:" );
    this->CompareViewBoxApplyButton = vtkKWPushButton::New ( );
    this->CompareViewBoxApplyButton->SetParent ( popUpFrame1 );
    this->CompareViewBoxApplyButton->Create ( );
    this->CompareViewBoxApplyButton->SetText ("Apply");    
    this->Script ( "grid %s -row 0 -column 0 -padx 2 -pady 8", rowsLabel->GetWidgetName());
    this->Script ( "grid %s -row 0 -column 1 -padx 6 -pady 8", this->CompareViewBoxRowEntry->GetWidgetName() );
    this->Script ( "grid %s -row 1 -column 0 -padx 2 -pady 8", columnsLabel->GetWidgetName());
    this->Script ( "grid %s -row 1 -column 1 -padx 6 -pady 8", this->CompareViewBoxColumnEntry->GetWidgetName() );
    this->Script ( "grid %s -row 2 -column 0 -columnspan 2 -pady 8", this->CompareViewBoxApplyButton->GetWidgetName() );
    // delete temporary stuff
    rowsLabel->Delete();
    columnsLabel->Delete();
    popUpFrame1->Delete();
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::PackLightboxView ( )
{
  /*
  // TO DO implement this.
  if ( this->GetApplication() != NULL )
    {
      vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();
      vtkSlicerColor *color = app->GetSlicerTheme()->GetSlicerColors ( );
      vtkSlicerGUILayout *layout = app->GetMainLayout ( );

      this->MainSlicerWindow->GetViewNotebook()->SetAlwaysShowTabs ( 0 );
      layout->SetCurrentViewArrangement ( vtkSlicerGUILayout::SlicerLayoutLightboxView );
    }
  */
}




//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::Save3DViewConfig ( )
{
  if ( this->ViewerWidget )
    {
      // TODO: Save the ViewerWidget's Camera Node
      this->ViewerWidget->GetMainViewer()->GetRenderer()->ComputeVisiblePropBounds ( this->MainRendererBBox );
    }
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::Restore3DViewConfig ( )
{
  if ( this->ViewerWidget )
    {
      // TODO: Restore the ViewerWidget's Camera Node
      this->ViewerWidget->GetMainViewer()->GetRenderer()->ResetCamera ( );
    }
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::ConfigureMainSliceViewers ( )
{
  if (this->SlicesGUI)
    {
      vtkSlicerSliceGUI *g = NULL;
      char *layoutname = NULL;
      int nSliceGUI = this->SlicesGUI->GetNumberOfSliceGUI();
      for (int i = 0; i < nSliceGUI; i++)
        {
          if (i == 0)
            {
              g = this->SlicesGUI->GetFirstSliceGUI();
              layoutname = this->SlicesGUI->GetFirstSliceGUILayoutName();
            }
          else
            {
              g = this->SlicesGUI->GetNextSliceGUI(layoutname);
              layoutname = this->SlicesGUI->GetNextSliceGUILayoutName(layoutname);
            }
          
          switch (layoutname[0])
            {
            case 'R':
              g->GetSliceController()->GetSliceNode()->SetOrientationToAxial();
              break;
            case 'Y':
              g->GetSliceController()->GetSliceNode()->SetOrientationToSagittal();
              break;
            case 'G':
              g->GetSliceController()->GetSliceNode()->SetOrientationToCoronal();
              break;
            default:
              g->GetSliceController()->GetSliceNode()->SetOrientationToSagittal();
            }
        }
    }
}



//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::AddMainSliceViewerObservers ( )
{
  if (this->SlicesGUI)
    {
      vtkSlicerSliceGUI *g = NULL;
      char *layoutname = NULL;
      int nSliceGUI = this->SlicesGUI->GetNumberOfSliceGUI();
      for (int i = 0; i < nSliceGUI; i++)
        {
          if (i == 0)
            {
              g = this->SlicesGUI->GetFirstSliceGUI();
              layoutname = this->SlicesGUI->GetFirstSliceGUILayoutName();
            }
          else
            {
              g = this->SlicesGUI->GetNextSliceGUI(layoutname);
              layoutname = this->SlicesGUI->GetNextSliceGUILayoutName(layoutname);
            }
          
          g->AddGUIObservers();
        }
    }
}



//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::RemoveMainSliceViewerObservers ( )
{
  if (this->SlicesGUI)
    {
      vtkSlicerSliceGUI *g = NULL;
      char *layoutname = NULL;
      int nSliceGUI = this->SlicesGUI->GetNumberOfSliceGUI();
      for (int i = 0; i < nSliceGUI; i++)
        {
          if (i == 0)
            {
              g = this->SlicesGUI->GetFirstSliceGUI();
              layoutname = this->SlicesGUI->GetFirstSliceGUILayoutName();
            }
          else
            {
              g = this->SlicesGUI->GetNextSliceGUI(layoutname);
              layoutname = this->SlicesGUI->GetNextSliceGUILayoutName(layoutname);
            }
          
          g->RemoveGUIObservers();
        }
    }
}


// to do: the three vtkSlicerSliceLogic pointers should be changed with a 
// pointer to a vtkSlicerSliceLogic map.
// current implementation only takes 3 for red, yellow, and green
//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::SetAndObserveMainSliceLogic ( vtkSlicerSliceLogic *l0,
                                                            vtkSlicerSliceLogic *l1,
                                                            vtkSlicerSliceLogic *l2 )
{
        //if (this->SlicesGUI)
        //{
        //      vtkSlicerSliceGUI *g = NULL;
        //      char *layoutname = NULL;
        //      int nSliceGUI = this->SlicesGUI->GetNumberOfSliceGUI();
        //      for (int i = 0; i < nSliceGUI; i++)
        //      {
        //              if (i == 0)
        //              {
        //                      g = this->SlicesGUI->GetFirstSliceGUI();
        //                      layoutname = this->SlicesGUI->GetFirstSliceGUILayoutName();
        //              }
        //              else
        //              {
        //                      g = this->SlicesGUI->GetNextSliceGUI(layoutname);
        //                      layoutname = this->SlicesGUI->GetNextSliceGUILayoutName(layoutname);
        //              }

        //              g->RemoveGUIObservers();
        //      }
        //}
  
  vtkSlicerSliceGUI *g = this->SlicesGUI->GetSliceGUI("Red");
  g->SetAndObserveModuleLogic(l0);
  g = this->SlicesGUI->GetSliceGUI("Yellow");
  g->SetAndObserveModuleLogic(l1);
  g = this->SlicesGUI->GetSliceGUI("Green");
  g->SetAndObserveModuleLogic(l2);
}



//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::PopulateModuleChooseList ( )
{
#ifndef TOOLBAR_DEBUG
  this->GetApplicationToolbar()->GetModuleChooseGUI()->Populate();
#endif
}



//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::PackFirstSliceViewerFrame ( )
{
  /*
  this->Script ("pack %s -side left  -expand y -fill both -padx 0 -pady 0", 
    this->DefaultSlice0Frame->GetWidgetName( ) );
  */
}



//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::BuildGUIFrames ( )
{

  if ( this->GetApplication() != NULL ) {
    // pointers for convenience
    vtkSlicerApplication *app = vtkSlicerApplication::SafeDownCast(this->GetApplication() );
    vtkSlicerGUILayout *layout = app->GetMainLayout ( );
    
    if ( this->MainSlicerWindow != NULL ) {
      
      this->MainSlicerWindow->GetMainPanelFrame()->SetWidth ( layout->GetDefaultGUIPanelWidth() );
      this->MainSlicerWindow->GetMainPanelFrame()->SetHeight ( layout->GetDefaultGUIPanelHeight() );
      this->MainSlicerWindow->GetMainPanelFrame()->SetReliefToSunken();
      
      this->TopFrame->SetParent ( this->MainSlicerWindow->GetMainPanelFrame ( ) );
      this->TopFrame->Create ( );
      this->TopFrame->SetHeight ( layout->GetDefaultTopFrameHeight ( ) );
      
      this->LogoFrame->SetParent ( this->TopFrame );
      this->LogoFrame->Create( );
      this->LogoFrame->SetHeight ( layout->GetDefaultTopFrameHeight ( ) );            
      
      this->DropShadowFrame->SetParent ( this->MainSlicerWindow->GetMainPanelFrame() );
      this->DropShadowFrame->Create ( );
      // why is the theme not setting this???
      this->DropShadowFrame->SetBackgroundColor ( 0.9, 0.9, 1.0);
      
      this->SlicesControlFrame->SetParent ( this->DropShadowFrame );
      this->SlicesControlFrame->Create( );
      this->SlicesControlFrame->ExpandFrame ( );
      this->SlicesControlFrame->SetLabelText ( "Manipulate Slice Views");
      this->SlicesControlFrame->GetFrame()->SetHeight ( layout->GetDefaultSlicesControlFrameHeight ( ) );
      
      this->ViewControlFrame->SetParent ( this->DropShadowFrame );
      this->ViewControlFrame->Create( );
      this->ViewControlFrame->ExpandFrame ( );
      this->ViewControlFrame->SetLabelText ( "Manipulate 3D View" );
      this->ViewControlFrame->GetFrame()->SetHeight (layout->GetDefaultViewControlFrameHeight ( ) );
      
      
      app->Script ( "pack %s -side top -fill x -padx 1 -pady 1", this->TopFrame->GetWidgetName() );
      app->Script ( "pack %s -side left -expand 1 -fill x -padx 1 -pady 1", this->LogoFrame->GetWidgetName() );
      app->Script ( "pack %s -side bottom -expand n -fill x -padx 1 -ipady 1 -pady 0", this->DropShadowFrame->GetWidgetName() );
      app->Script ( "pack %s -side bottom -expand n -fill x -padx 0 -ipady 5 -pady 2", this->ViewControlFrame->GetWidgetName() );
      app->Script ( "pack %s -side bottom -expand n -fill x -padx 0 -ipady 5 -pady 1", this->SlicesControlFrame->GetWidgetName() );
      
    }
  }
}


//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::ConfigureRemoteIOSettings()
{
  vtkMRMLScene *scene = this->GetMRMLScene();
  if ( this->GetApplication() != NULL )
    {
    vtkSlicerApplication *app = vtkSlicerApplication::SafeDownCast(this->GetApplication() );
    if ( scene != NULL )
      {
      //--- update CacheManager
      vtkCacheManager *cm = scene->GetCacheManager();
      if ( cm != NULL )
        {
        if ( strcmp (cm->GetRemoteCacheDirectory(), app->GetRemoteCacheDirectory() ))
          {
          cm->SetRemoteCacheDirectory (app->GetRemoteCacheDirectory() );
          }
        if ( cm->GetEnableForceRedownload() != app->GetEnableForceRedownload() )
          {
          cm->SetEnableForceRedownload (app->GetEnableForceRedownload() );
          }
        /*
          if ( cm->GetEnableRemoteCacheOverwriting () != app->GetEnableRemoteCacheOverwriting() )
          {
          cm->SetEnableRemoteCacheOverwriting (app->GetEnableRemoteCacheOverwriting() );
          }
        */
        if ( cm->GetRemoteCacheLimit() != app->GetRemoteCacheLimit() )
          {
          cm->SetRemoteCacheLimit (app->GetRemoteCacheLimit() );
          }
        if ( cm->GetRemoteCacheFreeBufferSize() != app->GetRemoteCacheFreeBufferSize() )
          {
          cm->SetRemoteCacheFreeBufferSize (app->GetRemoteCacheFreeBufferSize() );
          }
        }
      //---- update DataIOManager
      vtkDataIOManager *dm = scene->GetDataIOManager();
      if ( dm != NULL )
        {
        if ( dm->GetEnableAsynchronousIO() != app->GetEnableAsynchronousIO() )
          {
          dm->SetEnableAsynchronousIO (app->GetEnableAsynchronousIO() );
          }
        }

      //---- update application settings interface if required...
      if ( this->GetMainSlicerWindow() != NULL )
        {
        if ( this->GetMainSlicerWindow()->GetApplicationSettingsInterface() != NULL )
          {
          (vtkSlicerApplicationSettingsInterface::SafeDownCast (this->GetMainSlicerWindow()->GetApplicationSettingsInterface()))->UpdateRemoteIOSettings();
          }
        }
      }
    }
}




//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::UpdateRemoteIOConfigurationForRegistry()
{

  vtkMRMLScene *scene = this->GetMRMLScene();
  if ( this->GetApplication() != NULL )
    {
    vtkSlicerApplication *app = vtkSlicerApplication::SafeDownCast(this->GetApplication() );
    if ( scene != NULL )
      {
      vtkCacheManager *cm = scene->GetCacheManager();
      if ( cm != NULL )
        {
        app->SetRemoteCacheDirectory (cm->GetRemoteCacheDirectory() );
        app->SetEnableForceRedownload (cm->GetEnableForceRedownload() );
        //app->SetEnableRemoteCacheOverwriting (cm->GetEnableRemoteCacheOverwriting() );
        app->SetRemoteCacheLimit (cm->GetRemoteCacheLimit() );
        app->SetRemoteCacheFreeBufferSize (cm->GetRemoteCacheFreeBufferSize() );
        }
      vtkDataIOManager *dm = scene->GetDataIOManager();
      if ( dm != NULL )
        {
        app->SetEnableAsynchronousIO (dm->GetEnableAsynchronousIO() );
        }
      }
    }
}
