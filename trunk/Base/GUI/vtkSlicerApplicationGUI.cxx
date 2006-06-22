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

#include "vtkImplicitPlaneWidget.h"

#include "vtkKWApplication.h"
#include "vtkKWTclInteractor.h"
#include "vtkKWFrame.h"
#include "vtkKWMenu.h"
#include "vtkKWMenuButtonWithLabel.h"
#include "vtkKWLabel.h"
#include "vtkKWNotebook.h"
#include "vtkKWPushButton.h"
#include "vtkKWRenderWidget.h"
#include "vtkKWScale.h"
#include "vtkKWUserInterfacePanel.h"
#include "vtkKWWidget.h"
#include "vtkKWCheckButton.h"
#include "vtkKWEntry.h"
#include "vtkKWEntryWithLabel.h"
#include "vtkKWResourceUtilities.h"
#include "vtkKWToolbarSet.h"
#include "vtkKWToolbar.h"

#include "vtkKWSplitFrame.h"
#include "vtkKWUserInterfaceManagerNotebook.h"

#include "vtkSlicerWindow.h"
#include "vtkSlicerApplication.h"
#include "vtkSlicerApplicationGUI.h"
#include "vtkSlicerApplicationGUI.h"
#include "vtkSlicerApplicationLogic.h"
#include "vtkSlicerModuleGUI.h"
#include "vtkSlicerGUILayout.h"
#include "vtkSlicerTheme.h"
#include "vtkSlicerToolbarIcons.h"
#include "vtkSlicerLogoIcons.h"
#include "vtkSlicerModuleNavigationIcons.h"
#include "vtkSlicerViewControlIcons.h"

//---------------------------------------------------------------------------
vtkStandardNewMacro (vtkSlicerApplicationGUI);
vtkCxxRevisionMacro(vtkSlicerApplicationGUI, "$Revision: 1.0 $");


//---------------------------------------------------------------------------
vtkSlicerApplicationGUI::vtkSlicerApplicationGUI (  )
{
    //---  
    // widgets used in the Slice module
    //---

    //--- slicer main window
    this->MainSlicerWin = vtkSlicerWindow::New ( );

    //--- slicer toolbars
    this->ModulesToolbar = vtkKWToolbar::New ( );
    this->LoadSaveToolbar = vtkKWToolbar::New ( );
    this->ViewToolbar = vtkKWToolbar::New ( );
    
    //--- slicer icons
    this->SlicerLogoIcons = vtkSlicerLogoIcons::New ();
    this->SlicerViewControlIcons = vtkSlicerViewControlIcons::New();
    this->SlicerToolbarIcons = vtkSlicerToolbarIcons::New ();
    this->SlicerModuleNavigationIcons = vtkSlicerModuleNavigationIcons::New ();

    //--- logo widgets to which icons are assigned.
    this->SlicerLogoLabel = vtkKWLabel::New();

    //--- toolbar widgets to which icons are assigned.
    this->HomeIconButton = vtkKWPushButton::New ( );
    this->DataIconButton = vtkKWPushButton::New ( );
    this->VolumeIconButton = vtkKWPushButton::New ( );
    this->ModelIconButton = vtkKWPushButton::New ( );
    this->EditorIconButton = vtkKWPushButton::New ( );
    this->EditorToolboxIconButton = vtkKWPushButton::New ( );
    this->ColorIconButton = vtkKWPushButton::New ( );
    this->FiducialsIconButton = vtkKWPushButton::New ( );
    this->TransformIconButton = vtkKWPushButton::New ( );
    this->SaveSceneIconButton = vtkKWPushButton::New ( );
    this->LoadSceneIconButton = vtkKWPushButton::New ( );
    this->ConventionalViewIconButton = vtkKWPushButton::New ( );
    this->OneUp3DViewIconButton = vtkKWPushButton::New ( );
    this->OneUpSliceViewIconButton = vtkKWPushButton::New ( );
    this->FourUpViewIconButton = vtkKWPushButton::New ( );
    this->TabbedViewIconButton = vtkKWPushButton::New ( );
    this->LightBoxViewIconButton = vtkKWPushButton::New ( );
    
    // Control frames that comprise the Main Slicer GUI
    this->LogoFrame = vtkKWFrame::New();
    this->ModuleChooseFrame = vtkKWFrame::New();
    this->SliceControlFrame = vtkKWFrame::New();    
    this->ViewControlFrame = vtkKWFrame::New();    
    this->DefaultSlice0Frame = vtkKWFrame::New ();
    this->DefaultSlice1Frame = vtkKWFrame::New ();
    this->DefaultSlice2Frame = vtkKWFrame::New ();

    //--- ui for the ModuleChooseFrame,
    this->ModulesMenuButton = vtkKWMenuButton::New();
    this->ModulesLabel = vtkKWLabel::New();
    this->ModulesPrev = vtkKWPushButton::New ( );
    this->ModulesNext = vtkKWPushButton::New ( );
    this->ModulesHistory = vtkKWPushButton::New ( );
    
    //--- ui for the SliceControlframe.
    this->ToggleAnnotationButton = vtkKWPushButton::New ( );
    this->ToggleFgBgButton = vtkKWPushButton::New ( );
    this->SliceFadeScale = vtkKWScale::New ( );
    this->SliceOpacityScale = vtkKWScale::New ( );
    
    //--- ui for the ViewControlFrame
    this->SpinButton = vtkKWCheckButton::New ( );
    this->RockButton = vtkKWCheckButton::New ( );
    this->OrthoButton = vtkKWCheckButton::New ( );
    this->CenterButton = vtkKWPushButton::New ( );
    this->SelectButton = vtkKWMenuButton::New ( );
    this->FOVEntry = vtkKWEntryWithLabel::New ( );

    //--- ui for the ViewControlFrame
    this->RotateAroundAIconButton = vtkKWLabel::New ( );
    this->RotateAroundPIconButton = vtkKWLabel::New ( );
    this->RotateAroundRIconButton = vtkKWLabel::New ( );
    this->RotateAroundLIconButton = vtkKWLabel::New ( );
    this->RotateAroundSIconButton = vtkKWLabel::New ( );
    this->RotateAroundIIconButton = vtkKWLabel::New ( );
    this->RotateAroundMiddleIconButton = vtkKWLabel::New ( );
    this->RotateAroundTopCornerIconButton = vtkKWLabel::New ( );
    this->RotateAroundBottomCornerIconButton = vtkKWLabel::New ( );

    this->LookFromAIconButton = vtkKWLabel::New ( );
    this->LookFromPIconButton = vtkKWLabel::New ( );
    this->LookFromRIconButton = vtkKWLabel::New ( );
    this->LookFromLIconButton = vtkKWLabel::New ( );
    this->LookFromSIconButton = vtkKWLabel::New ( );
    this->LookFromIIconButton = vtkKWLabel::New ( );
    this->LookFromMiddleIconButton = vtkKWLabel::New ( );
    this->LookFromTopCornerIconButton = vtkKWLabel::New ( );
    this->LookFromBottomCornerIconButton = vtkKWLabel::New ( );

    this->NavZoomInIconButton = vtkKWPushButton::New ( );
    this->NavZoomOutIconButton = vtkKWPushButton::New ( );
    this->NavZoomScale = vtkKWScale::New ( );
    
    //--- main viewer 
    this->MainViewer = vtkKWRenderWidget::New ( );
    this->PlaneWidget = NULL;

    this->LoadSceneDialog = vtkKWLoadSaveDialog::New();
    this->SaveSceneDialog = vtkKWLoadSaveDialog::New();   
}



//---------------------------------------------------------------------------
vtkSlicerApplicationGUI::~vtkSlicerApplicationGUI ( )
{


    if ( this->SlicerLogoIcons ) {
        this->SlicerLogoIcons->Delete ( );
        this->SlicerLogoIcons = NULL;
    }
    if ( this->SlicerViewControlIcons ) {
        this->SlicerViewControlIcons->Delete ( );
        this->SlicerViewControlIcons = NULL;
    }
    if ( this->SlicerToolbarIcons ) {
        this->SlicerToolbarIcons->Delete ( );
        this->SlicerToolbarIcons = NULL;
    }
    if ( this->SlicerModuleNavigationIcons ) {
        this->SlicerModuleNavigationIcons->Delete ( );
        this->SlicerModuleNavigationIcons = NULL;
    }

    this->DeleteGUIPanelWidgets ( );
    this->DeleteToolbarWidgets ( );

    vtkSlicerWindow *win = this->MainSlicerWin;
    if ( win ) {
        vtkKWToolbarSet *tbs = win->GetMainToolbarSet();
        if (tbs ) {
            tbs->RemoveAllToolbars () ;
        }
    }
    if ( this->ModulesToolbar ) {
        this->ModulesToolbar->Delete ( );
        this->ModulesToolbar = NULL;
    }
    if ( this->LoadSaveToolbar ) {
        this->LoadSaveToolbar->Delete ( );
        this->LoadSaveToolbar = NULL;
    }
    if ( this->ViewToolbar ) {
        this->ViewToolbar->Delete ( );
        this->ViewToolbar = NULL;
    }

    this->DeleteFrames ( );

    if ( this->MainViewer ) {
        this->MainViewer->RemoveAllViewProps ( );
        this->MainViewer->Delete ( );
        this->MainViewer = NULL;
    }
    if ( this->PlaneWidget ) {
        this->PlaneWidget->Delete ( );
        this->PlaneWidget = NULL;
    }
    if ( this->LoadSceneDialog ) {
        this->LoadSceneDialog->Delete();
        this->LoadSceneDialog = NULL;
    }
    if ( this->SaveSceneDialog ) {
        this->SaveSceneDialog->Delete();
        this->SaveSceneDialog = NULL;
    }
    if ( this->MainSlicerWin ) {
        this->MainSlicerWin->Delete ( );
        this->MainSlicerWin = NULL;
    }
}


//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::PrintSelf ( ostream& os, vtkIndent indent )
{
    this->vtkObject::PrintSelf ( os, indent );

    os << indent << "SlicerApplicationGUI: " << this->GetClassName ( ) << "\n";
    os << indent << "MainViewer: " << this->GetMainViewer ( ) << "\n";
    os << indent << "MainSlicerWin: " << this->GetMainSlicerWin ( ) << "\n";
    // print widgets?
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::ProcessLoadSceneCommand()
{
    this->LoadSceneDialog->RetrieveLastPathFromRegistry(
      "OpenPath");

    this->LoadSceneDialog->Invoke();
    // If a file has been selected for loading...
    char *fileName = this->LoadSceneDialog->GetFileName();
    if ( fileName ) 
      {
        if (this->GetMRMLScene()) 
          {
          this->GetMRMLScene()->SetURL(fileName);
          this->GetMRMLScene()->Connect();
          this->LoadSceneDialog->SaveLastPathToRegistry("OpenPath");
          }
      }
    return;
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::ProcessSaveSceneCommand()
{
    this->SaveSceneDialog->RetrieveLastPathFromRegistry(
      "OpenPath");

     this->SaveSceneDialog->Invoke();

    // If a file has been selected for loading...
    char *fileName = this->SaveSceneDialog->GetFileName();
    if ( fileName ) 
      {
        if (this->GetMRMLScene()) 
          {
          this->GetMRMLScene()->SetURL(fileName);
          this->GetMRMLScene()->Commit();  
          this->SaveSceneDialog->SaveLastPathToRegistry("OpenPath");
          }
      }
    return;
}    

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::AddGUIObservers ( )
{

    // add observers onto the buttons and menubutton in the SlicerControl frame
    this->HomeIconButton->AddObserver (vtkKWPushButton::InvokedEvent, (vtkCommand *)this->GUICallbackCommand );
    this->DataIconButton->AddObserver (vtkKWPushButton::InvokedEvent, (vtkCommand *)this->GUICallbackCommand );
    this->VolumeIconButton->AddObserver (vtkKWPushButton::InvokedEvent, (vtkCommand *)this->GUICallbackCommand );
    this->ModelIconButton->AddObserver (vtkKWPushButton::InvokedEvent, (vtkCommand *)this->GUICallbackCommand );
    this->ModulesMenuButton->AddObserver (vtkCommand::ModifiedEvent, (vtkCommand *)this->GUICallbackCommand );

    this->GetMainSlicerWin()->GetFileMenu()->AddObserver (vtkKWMenu::MenuItemInvokedEvent, (vtkCommand *)this->GUICallbackCommand );
    
    this->LoadSceneDialog->AddObserver ( vtkCommand::ModifiedEvent, (vtkCommand *)this->GUICallbackCommand );
    this->SaveSceneDialog->AddObserver ( vtkCommand::ModifiedEvent, (vtkCommand *)this->GUICallbackCommand );

    this->SliceFadeScale->AddObserver ( vtkKWScale::ScaleValueStartChangingEvent, (vtkCommand *)this->GUICallbackCommand );
    this->SliceFadeScale->AddObserver ( vtkKWScale::ScaleValueChangingEvent, (vtkCommand *)this->GUICallbackCommand );

    this->ToggleFgBgButton->AddObserver ( vtkKWPushButton::InvokedEvent, (vtkCommand *)this->GUICallbackCommand );
}


//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::RemoveGUIObservers ( )
{
    this->HomeIconButton->RemoveObservers (vtkKWPushButton::InvokedEvent, (vtkCommand *)this->GUICallbackCommand );
    this->DataIconButton->RemoveObservers (vtkKWPushButton::InvokedEvent, (vtkCommand *)this->GUICallbackCommand );
    this->VolumeIconButton->RemoveObservers (vtkKWPushButton::InvokedEvent, (vtkCommand *)this->GUICallbackCommand );
    this->ModelIconButton->RemoveObservers (vtkKWPushButton::InvokedEvent, (vtkCommand *)this->GUICallbackCommand );
    this->ModulesMenuButton->RemoveObservers (vtkCommand::ModifiedEvent, (vtkCommand *)this->GUICallbackCommand );

    this->LoadSceneDialog->RemoveObservers ( vtkCommand::ModifiedEvent, (vtkCommand *) this->GUICallbackCommand );
    this->SaveSceneDialog->RemoveObservers ( vtkCommand::ModifiedEvent, (vtkCommand *) this->GUICallbackCommand );
    this->GetMainSlicerWin()->GetFileMenu()->RemoveObservers ( vtkKWMenu::MenuItemInvokedEvent, (vtkCommand *)this->GUICallbackCommand );
    this->SliceFadeScale->RemoveObservers ( vtkKWScale::ScaleValueStartChangingEvent, (vtkCommand *)this->GUICallbackCommand );
    this->SliceFadeScale->RemoveObservers ( vtkKWScale::ScaleValueChangingEvent, (vtkCommand *)this->GUICallbackCommand );    
    this->ToggleFgBgButton->RemoveObservers ( vtkKWPushButton::InvokedEvent, (vtkCommand *)this->GUICallbackCommand );
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
    vtkSlicerModuleGUI * m;
    const char *mName;
    vtkKWPushButton *pushb = vtkKWPushButton::SafeDownCast (caller );
    vtkKWMenuButton *menub = vtkKWMenuButton::SafeDownCast (caller );
    vtkKWMenu *menu = vtkKWMenu::SafeDownCast (caller );
    vtkKWLoadSaveDialog *filebrowse = vtkKWLoadSaveDialog::SafeDownCast(caller);
    vtkKWScale *scale = vtkKWScale::SafeDownCast(caller);

    vtkSlicerApplication *app = vtkSlicerApplication::SafeDownCast( this->GetApplication() );
        
    // Process events from top row of buttons
    // For now, Home button takes us to the Volumes module.
    if ( pushb == this->HomeIconButton && event == vtkKWPushButton::InvokedEvent ) {
        vtkSlicerModuleGUI *m = vtkSlicerApplication::SafeDownCast(
          this->GetApplication())->GetModuleGUIByName("Volumes");
        if ( m != NULL ) { m->GetUIPanel()->Raise(); }
        this->ModulesMenuButton->SetValue ( "Volumes" );
    }
    else if (pushb == this->DataIconButton && event == vtkKWPushButton::InvokedEvent ) {
        vtkSlicerModuleGUI *m = vtkSlicerApplication::SafeDownCast(this->GetApplication())->GetModuleGUIByName("Data");
        if ( m != NULL ) { m->GetUIPanel()->Raise(); }
        this->ModulesMenuButton->SetValue ( "Data" );
    }
    else if (pushb == this->VolumeIconButton && event == vtkKWPushButton::InvokedEvent ) {
        vtkSlicerModuleGUI *m = vtkSlicerApplication::SafeDownCast(this->GetApplication())->GetModuleGUIByName("Volumes");
        if ( m != NULL ) { m->GetUIPanel()->Raise(); }
        this->ModulesMenuButton->SetValue ( "Volumes" );
    }
    else if (pushb == this->ModelIconButton && event == vtkKWPushButton::InvokedEvent ) {
        vtkSlicerModuleGUI *m = vtkSlicerApplication::SafeDownCast(this->GetApplication())->GetModuleGUIByName("Models");
        if ( m != NULL ) { m->GetUIPanel()->Raise(); }
        this->ModulesMenuButton->SetValue ( "Models" );
    }
    else if (pushb == this->TransformIconButton && event == vtkKWPushButton::InvokedEvent ) {
        //vtkSlicerModuleGUI *m = vtkSlicerApplication::SafeDownCast(this->GetApplication())->GetModuleGUIByName("Transformments");
        //if ( m != NULL ) { m->GetUIPanel()->Raise(); }
        this->ModulesMenuButton->SetValue ( "Transform" );
    }
    else if (menu == this->GetMainSlicerWin()->GetFileMenu() && event == vtkKWMenu::MenuItemInvokedEvent)
    {
      int index = (int) (*((int *)callData));
      if (index == 2)
        {
          // use command directly instead of this
          //this->ProcessLoadSceneCommand()
        }
      else if (index == 3)
        {
          // use command directly instead of this
          //this->ProcessSaveSceneCommand()
        }
    }

    //--- Process events from menubutton
    //--- TODO: change the Logic's "active module" and raise the appropriate UIPanel.
    if ( menub == this->ModulesMenuButton && event == vtkCommand::ModifiedEvent )
        {
            if ( app->GetModuleGUICollection ( ) != NULL )
                {
                    app->GetModuleGUICollection( )->InitTraversal( );
                    m = vtkSlicerModuleGUI::SafeDownCast( app->GetModuleGUICollection( )->GetNextItemAsObject( ) );
                    while (m != NULL )
                        {
                            mName = m->GetUIPanel()->GetName();
                            if ( !strcmp (this->ModulesMenuButton->GetValue(), mName) ) {
                                m->GetUIPanel()->Raise();
                                break;
                            }
                            m = vtkSlicerModuleGUI::SafeDownCast( app->GetModuleGUICollection( )->GetNextItemAsObject( ) );
                        }
                    //this->ModulesMenuButton->SetValue ( "Modules" );
                }
        }

    // Process the Fade scale and button
    // -- set save state when manipulation starts
    // -- toggle the value if needed
    // -- adjust the Opacity of every composite node on every event
    if ( scale == this->SliceFadeScale && event == vtkKWScale::ScaleValueStartChangingEvent ||
         pushb == this->ToggleFgBgButton && event == vtkKWPushButton::InvokedEvent )
      {
      if (this->GetMRMLScene()) 
        {
        this->GetMRMLScene()->SaveStateForUndo();
        }
      }

    if ( scale == this->SliceFadeScale && event == vtkKWScale::ScaleValueChangingEvent ||
         pushb == this->ToggleFgBgButton && event == vtkKWPushButton::InvokedEvent )
      {

      if ( pushb == this->ToggleFgBgButton && event == vtkKWPushButton::InvokedEvent ) 
        {
        this->SliceFadeScale->SetValue( 1.0 - this->SliceFadeScale->GetValue() );
        }

      int i, nnodes = this->MRMLScene->GetNumberOfNodesByClass("vtkMRMLSliceCompositeNode");
      vtkMRMLSliceCompositeNode *cnode;
      for (i = 0; i < nnodes; i++)
        {
        cnode = vtkMRMLSliceCompositeNode::SafeDownCast (
                this->MRMLScene->GetNthNodeByClass( i, "vtkMRMLSliceCompositeNode" ) );
        cnode->SetOpacity( this->SliceFadeScale->GetValue() );
        }
      }

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
void vtkSlicerApplicationGUI::BuildGUI ( )
{
    int i;
    
    // Set up the conventional window: 3Dviewer, slice widgets, UI panel for now.
    if ( this->GetApplication() != NULL ) {
        vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();

        app->GetMainLayout()->InitializeLayout ( );

        if ( this->MainSlicerWin != NULL ) {

            // set up Slicer's main window
            this->MainSlicerWin->SecondaryPanelVisibilityOn ( );
            this->MainSlicerWin->MainPanelVisibilityOn ( );
            app->AddWindow ( this->MainSlicerWin );

            // Create the console before the window
            // - this will make the console independent of the main window
            //   so it can be raised/lowered independently
            this->MainSlicerWin->GetTclInteractor()->SetApplication(app);
            this->MainSlicerWin->GetTclInteractor()->Create();

            this->MainSlicerWin->Create ( );        

            // Construct menu bar and set up global key bindings


            this->GetMainSlicerWin()->GetFileMenu()->InsertCommand (this->GetMainSlicerWin()->GetFileMenuInsertPosition(),
                                              "Load Scene", this, "ProcessLoadSceneCommand");
            this->GetMainSlicerWin()->GetFileMenu()->InsertCommand (this->GetMainSlicerWin()->GetFileMenuInsertPosition(),
                                               "Save Scene", this, "ProcessSaveSceneCommand");

            this->GetMainSlicerWin()->GetFileMenu()->InsertSeparator (
                this->GetMainSlicerWin()->GetFileMenuInsertPosition());

            i = this->MainSlicerWin->GetEditMenu()->AddCommand ("Set Home", NULL, NULL);
            this->MainSlicerWin->GetEditMenu()->SetItemAccelerator ( i, "Ctrl+H");
            i = this->MainSlicerWin->GetEditMenu()->AddCommand ( "Undo", NULL, "$::slicer3::MRMLScene Undo" );
            this->MainSlicerWin->GetEditMenu()->SetItemAccelerator ( i, "Ctrl+Z");
            i = this->MainSlicerWin->GetEditMenu()->AddCommand ( "Redo", NULL, "$::slicer3::MRMLScene Redo" );
            this->MainSlicerWin->GetEditMenu()->SetItemAccelerator ( i, "Ctrl+Y");
            //i = this->MainSlicerWin->GetViewMenu()->AddCommand ( ? );
            //i = this->MainSlicerWin->GetWindowMenu()->AddCommand ( ? );
            //i = this->MainSlicerWin->GetHelpMenu()->AddCommand ( ? );

            // configure default size of GUI
            this->ConfigureMainSlicerWindow ( );
            this->ConfigureMainViewerPanel ( );
            this->ConfigureSliceViewersPanel ( );
            this->ConfigureGUIPanel ( );

            // Populate toolbar
            this->BuildToolBar();

            // Build 3DViewer
            this->BuildMainViewer ( );

            // Build main GUI panel
            this->BuildLogoGUIPanel ( );
            this->BuildSlicerControlGUIPanel ( );

            // Turn off the tabs for pages in the ModuleControlGUI
            this->MainSlicerWin->GetMainNotebook()->ShowIconsOff ( );
            //this->MainSlicerWin->GetMainNotebook()->SetAlwaysShowTabs ( 0 );
            this->MainSlicerWin->GetMainNotebook()->SetUseFrameWithScrollbars ( 1 );
            
            this->BuildSliceControlGUIPanel ( );
            this->BuildViewControlGUIPanel ( );

            this->LoadSceneDialog->SetParent ( this->MainSlicerWin );
            this->LoadSceneDialog->Create ( );
            this->LoadSceneDialog->SetFileTypes("{ {MRML Scene} {*.mrml} }");
            this->LoadSceneDialog->RetrieveLastPathFromRegistry("OpenPath");

            this->SaveSceneDialog->SetParent ( this->MainSlicerWin );
            this->SaveSceneDialog->Create ( );
            this->SaveSceneDialog->SetFileTypes("{ {MRML Scene} {*.mrml} }");
            this->SaveSceneDialog->SaveDialogOn();
            this->SaveSceneDialog->RetrieveLastPathFromRegistry("OpenPath");
        }
    }
}


//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::DisplayMainSlicerWindow ( )
{

    this->MainSlicerWin->Display ( );
}





//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::DeleteToolbarWidgets ( )
{

    if ( this->ModulesToolbar ) {
        this->ModulesToolbar->RemoveAllWidgets( );
    }
    if ( this->LoadSaveToolbar ) {
        this->LoadSaveToolbar->RemoveAllWidgets ( );
    }
    if ( this->ViewToolbar ) {
        this->ViewToolbar->RemoveAllWidgets ( );
    }

    if ( this->HomeIconButton ) {
        this->HomeIconButton->Delete ( );
        this->HomeIconButton = NULL;
    }
    if ( this->DataIconButton ) {
        this->DataIconButton->Delete ( );
        this->DataIconButton = NULL;
    }
    if ( this->VolumeIconButton ) {
        this->VolumeIconButton->Delete ( );
        this->VolumeIconButton = NULL;
    }
    if ( this->ModelIconButton ) {
        this->ModelIconButton->Delete ( );
        this->ModelIconButton = NULL;
    }
    if ( this->EditorIconButton ) {
        this->EditorIconButton->Delete ( );
        this->EditorIconButton = NULL;
    }
    if ( this->EditorToolboxIconButton ) {
        this->EditorToolboxIconButton->Delete ( );
        this->EditorToolboxIconButton = NULL;
    }
    if ( this->TransformIconButton ) {
        this->TransformIconButton->Delete ( );
        this->TransformIconButton = NULL;
    }
    if ( this->ColorIconButton ) {
        this->ColorIconButton->Delete ( );
        this->ColorIconButton = NULL;
    }
    if ( this->FiducialsIconButton ) {
        this->FiducialsIconButton->Delete ( );
        this->FiducialsIconButton = NULL;
    }
    if ( this->SaveSceneIconButton ) {
        this->SaveSceneIconButton->Delete ( );
        this->SaveSceneIconButton = NULL;
    }
    if ( this->LoadSceneIconButton ) {
        this->LoadSceneIconButton->Delete ( );
        this->LoadSceneIconButton = NULL;
    }
    if ( this->ConventionalViewIconButton ) {
        this->ConventionalViewIconButton->Delete ( );
        this->ConventionalViewIconButton = NULL;
    }
    if ( this->OneUp3DViewIconButton ) {
        this->OneUp3DViewIconButton->Delete ( );
        this->OneUp3DViewIconButton = NULL;
    }
    if ( this->OneUpSliceViewIconButton ) {
        this->OneUpSliceViewIconButton->Delete ( );
        this->OneUpSliceViewIconButton = NULL;
    }
    if ( this->FourUpViewIconButton ) {
        this->FourUpViewIconButton->Delete ( );
        this->FourUpViewIconButton = NULL;
    }
    if ( this->TabbedViewIconButton ) {
        this->TabbedViewIconButton->Delete ( );
        this->TabbedViewIconButton = NULL;
    }
    if ( this->LightBoxViewIconButton ) {
        this->LightBoxViewIconButton->Delete ( );
        this->LightBoxViewIconButton = NULL;
    }
    
}


//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::DeleteGUIPanelWidgets ( )
{
    //--- widgets from the ModuleChooseFrame
    if ( this->ModulesMenuButton ) {
        this->ModulesMenuButton->Delete();
        this->ModulesMenuButton = NULL;
    }
    if ( this->ModulesLabel ) {
        this->ModulesLabel->Delete ( );
        this->ModulesLabel = NULL;
    }
    if ( this->ModulesPrev ) {
        this->ModulesPrev->Delete ( );
        this->ModulesPrev = NULL;
    }
    if ( this->ModulesNext ) {
        this->ModulesNext->Delete ( );
        this->ModulesNext = NULL;
    }
    if ( this->ModulesHistory) {
        this->ModulesHistory->Delete ( );
        this->ModulesHistory = NULL;
    }

    //--- widgets from ViewControlFrame
    if ( this->SpinButton ) {
        this->SpinButton->Delete();
        this->SpinButton = NULL;
    }
    if ( this->RockButton) {
        this->RockButton->Delete();
        this->RockButton = NULL;
    }
    if ( this->OrthoButton ) {
        this->OrthoButton->Delete();
        this->OrthoButton = NULL;
    }
    if ( this->CenterButton ) {
        this->CenterButton->Delete();
        this->CenterButton = NULL;
    }
    if ( this->SelectButton ) {
        this->SelectButton->Delete();
        this->SelectButton = NULL;
    }
    if ( this->FOVEntry ) {
        this->FOVEntry->Delete();
        this->FOVEntry= NULL;
    }

    //--- widgets from LogoFrame
    if (this->SlicerLogoLabel ) {
        this->SlicerLogoLabel->Delete();
        this->SlicerLogoLabel = NULL;
    }

    //--- widgets from the SliceControlFrame
    if ( this->ToggleAnnotationButton ) {
        this->ToggleAnnotationButton->Delete ( );
        this->ToggleAnnotationButton = NULL;
    }
    if ( this->ToggleFgBgButton ) {
        this->ToggleFgBgButton->Delete ( );
        this->ToggleFgBgButton = NULL;
    }
    if ( this->SliceFadeScale ) {
        this->SliceFadeScale->Delete ( );
        this->SliceFadeScale = NULL;
    }
    if ( this->SliceOpacityScale ) {
        this->SliceOpacityScale->Delete ( );
        this->SliceOpacityScale = NULL;
    }

    //--- widgets from the ViewControlFrame
    if ( this->RotateAroundAIconButton ) {
        this->RotateAroundAIconButton->Delete ( );
        this->RotateAroundAIconButton = NULL;
    }
    if ( this->RotateAroundPIconButton ) {
        this->RotateAroundPIconButton->Delete ( );
        this->RotateAroundPIconButton = NULL;
    }
    if ( this->RotateAroundRIconButton ) {
        this->RotateAroundRIconButton->Delete ( );
        this->RotateAroundRIconButton = NULL;
    }
    if ( this->RotateAroundLIconButton ) {
        this->RotateAroundLIconButton->Delete ( );
        this->RotateAroundLIconButton = NULL;
    }
    if ( this->RotateAroundSIconButton ) {
        this->RotateAroundSIconButton->Delete ( );
        this->RotateAroundSIconButton = NULL;
    }
    if ( this->RotateAroundIIconButton ) {
        this->RotateAroundIIconButton->Delete ( );
        this->RotateAroundIIconButton = NULL;
    }
    if ( this->RotateAroundMiddleIconButton ) {
        this->RotateAroundMiddleIconButton->Delete ( );
        this->RotateAroundMiddleIconButton = NULL;
    }
    if ( this->RotateAroundTopCornerIconButton ) {
        this->RotateAroundTopCornerIconButton->Delete ( );
        this->RotateAroundTopCornerIconButton = NULL;
    }
    if ( this->RotateAroundBottomCornerIconButton ) {
        this->RotateAroundBottomCornerIconButton->Delete ( );
        this->RotateAroundBottomCornerIconButton = NULL;
    }
    if ( this->LookFromAIconButton ) {
        this->LookFromAIconButton->Delete ( );
        this->LookFromAIconButton = NULL;
    }
    if ( this->LookFromPIconButton ) {
        this->LookFromPIconButton->Delete ( );
        this->LookFromPIconButton = NULL;
    }
    if ( this->LookFromRIconButton ) {
        this->LookFromRIconButton->Delete ( );
        this->LookFromRIconButton = NULL;
    }
    if ( this->LookFromLIconButton ) {
        this->LookFromLIconButton->Delete ( );
        this->LookFromLIconButton = NULL;
    }
    if ( this->LookFromSIconButton ) {
        this->LookFromSIconButton->Delete ( );
        this->LookFromSIconButton = NULL;
    }
    if ( this->LookFromIIconButton ) {
        this->LookFromIIconButton->Delete ( );
        this->LookFromIIconButton = NULL;
    }
    if ( this->LookFromMiddleIconButton ) {
        this->LookFromMiddleIconButton->Delete ( );
        this->LookFromMiddleIconButton = NULL;
    }
    if ( this->LookFromTopCornerIconButton ) {
        this->LookFromTopCornerIconButton->Delete ( );
        this->LookFromTopCornerIconButton = NULL;
    }
    if ( this->LookFromBottomCornerIconButton ) {
        this->LookFromBottomCornerIconButton->Delete ( );
        this->LookFromBottomCornerIconButton = NULL;
    }
    if ( this->NavZoomInIconButton ) {
        this->NavZoomInIconButton->Delete ( );
        this->NavZoomInIconButton = NULL;
    }
    if ( this->NavZoomOutIconButton ) {
        this->NavZoomOutIconButton->Delete ( );
        this->NavZoomOutIconButton = NULL;
    }
    if ( this->NavZoomScale ) {
        this->NavZoomScale->Delete ( );
        this->NavZoomScale = NULL;
    }
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::DeleteFrames ( )
{
    if ( this->LogoFrame ) {
        this->LogoFrame->Delete ();
        this->LogoFrame = NULL;
    }
    if ( this->ModuleChooseFrame ) {
        this->ModuleChooseFrame->Delete ();
        this->ModuleChooseFrame = NULL;
    }
    if ( this->SliceControlFrame ) {
        this->SliceControlFrame->Delete ( );
        this->SliceControlFrame = NULL;
    }
    if ( this->ViewControlFrame ) {
        this->ViewControlFrame->Delete ( );
        this->ViewControlFrame = NULL;
    }
    if ( this->DefaultSlice0Frame ) {
        this->DefaultSlice0Frame->Delete ();
        this->DefaultSlice0Frame = NULL;
    }
    if ( this->DefaultSlice1Frame ) {
        this->DefaultSlice1Frame->Delete ();
        this->DefaultSlice1Frame = NULL;
    }
    if ( this->DefaultSlice2Frame ) {
        this->DefaultSlice2Frame->Delete ();
        this->DefaultSlice2Frame = NULL;
    }
}




//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::BuildMainViewer ( )
{

    if ( this->GetApplication() != NULL ) {
        vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();

        vtkSlicerWindow *win = this->MainSlicerWin;
        if ( this->MainViewer != NULL ) {
            this->MainViewer->SetParent (win->GetViewPanelFrame ( ) );
            this->MainViewer->Create ( );
            app->Script  ("pack %s -side top -fill both -expand y -padx 0 -pady 0",
                          this->MainViewer->GetWidgetName ( ) );
            vtkSlicerApplication *app = vtkSlicerApplication::SafeDownCast(this->GetApplication()); 
            this->MainViewer->SetRendererBackgroundColor ( 
                app->GetSlicerTheme()->GetSlicerColors()->ViewerBlue );
            this->MainViewer->GetRenderer()->GetActiveCamera()->ParallelProjectionOff();

            // set up antialiasing
            this->MainViewer->GetRenderWindow()->LineSmoothingOn();
            this->MainViewer->GetRenderWindow()->PolygonSmoothingOn ( );
            this->MainViewer->GetRenderWindow()->PointSmoothingOn();
            // this->MainViewer->SetMultiSamples ( 4 );
            
            // put in a plane interactor to test
            vtkCubeSource *cubeSource = vtkCubeSource::New();
            vtkPolyDataMapper *cubeMapper = vtkPolyDataMapper::New ();
            cubeMapper->SetInputConnection ( cubeSource->GetOutputPort() );
            vtkActor *cubeActor = vtkActor::New ( );
            cubeActor->SetMapper ( cubeMapper );
            // don't add the actor, so we can see the interactor
            MainViewer->AddViewProp ( cubeActor );

            // TODO: this requires a change to KWWidgets
            this->PlaneWidget = vtkImplicitPlaneWidget::New();
            // this->PlaneWidget->SetInteractor( this->GetRenderWindowInteractor() );
            this->PlaneWidget->SetInteractor( this->GetRenderWindowInteractor() );
            this->PlaneWidget->PlaceWidget();
            this->PlaneWidget->On();

            MainViewer->ResetCamera ( );
        
            cubeSource->Delete ();
            cubeActor->Delete ();
            cubeMapper->Delete ();
        }
    }
}


//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::BuildToolBar()
{
    if ( this->GetApplication() != NULL ) {
        vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();

        //--- configure the window's main toolbarset.
        vtkSlicerWindow *win = this->MainSlicerWin;
        vtkKWToolbarSet *tbs = win->GetMainToolbarSet();
        tbs->SetToolbarsWidgetsFlatAspect ( 1 );
        tbs->BottomSeparatorVisibilityOn ( );
        tbs->TopSeparatorVisibilityOn ( );

        //--- configure toolbars
        vtkKWToolbar *mtb = this->GetModulesToolbar ( );
        mtb->SetParent ( tbs->GetToolbarsFrame ( ) );
        mtb->Create();
        mtb->SetReliefToGroove ( );
        mtb->SetWidgetsPadX ( 3 );
        mtb->SetWidgetsPadY ( 2 );

        vtkKWToolbar *ltb = this->GetLoadSaveToolbar ( );
        ltb->SetParent ( tbs->GetToolbarsFrame ( ) );
        ltb->Create();
        ltb->SetReliefToGroove ( );
        ltb->SetWidgetsPadX ( 3 );
        ltb->SetWidgetsPadY ( 2 );

        vtkKWToolbar *vtb = this->GetViewToolbar ( );
        vtb->SetParent ( tbs->GetToolbarsFrame ( ) );
        vtb->Create();
        vtb->SetReliefToGroove ( );
        vtb->SetWidgetsPadX ( 3 );
        vtb->SetWidgetsPadY ( 2 );

        //--- and add toolbars to the window's main toolbar set.        
        tbs->AddToolbar ( this->GetModulesToolbar() );
        tbs->AddToolbar ( this->GetViewToolbar() );
        tbs->AddToolbar ( this->GetLoadSaveToolbar() );

        //--- create icons and the labels that display them and add to toolbar

        // home icon
        this->HomeIconButton->SetParent ( mtb->GetFrame ( ));
        this->HomeIconButton->Create ( );
        this->HomeIconButton->SetReliefToFlat ( );
        this->HomeIconButton->SetBorderWidth ( 0 );
        this->HomeIconButton->SetOverReliefToNone ( );
        this->HomeIconButton->SetImageToIcon ( this->SlicerToolbarIcons->GetHomeIcon( ) );
        this->HomeIconButton->SetBalloonHelpString ( "Home" );
        mtb->AddWidget ( this->HomeIconButton );

        // data module icon
        this->DataIconButton->SetParent ( mtb->GetFrame ( ));
        this->DataIconButton->Create ( );
        this->DataIconButton->SetReliefToFlat ( );
        this->DataIconButton->SetBorderWidth ( 0 );
        this->DataIconButton->SetOverReliefToNone ( );
        this->DataIconButton->SetImageToIcon ( this->SlicerToolbarIcons->GetDataIcon ( ) );
        this->DataIconButton->SetBalloonHelpString ( "Data");
        mtb->AddWidget ( this->DataIconButton );

        // volume module icon
        this->VolumeIconButton->SetParent ( mtb->GetFrame ( ));
        this->VolumeIconButton->Create ( );
        this->VolumeIconButton->SetReliefToFlat ( );
        this->VolumeIconButton->SetBorderWidth ( 0 );
        this->VolumeIconButton->SetOverReliefToNone ( );
        this->VolumeIconButton->SetImageToIcon ( this->SlicerToolbarIcons->GetVolumeIcon ( ));
        this->VolumeIconButton->SetBalloonHelpString ( "Volumes");
        mtb->AddWidget ( this->VolumeIconButton );

        // models module icon
        this->ModelIconButton->SetParent (mtb->GetFrame ( ) );
        this->ModelIconButton->Create ( );
        this->ModelIconButton->SetReliefToFlat ( );
        this->ModelIconButton->SetBorderWidth ( 0 );
        this->ModelIconButton->SetOverReliefToNone ( );
        this->ModelIconButton->SetImageToIcon ( this->SlicerToolbarIcons->GetModelIcon ( ) );
        this->ModelIconButton->SetBalloonHelpString ( "Models");
        mtb->AddWidget ( this->ModelIconButton );

        // transforms module icon
        this->TransformIconButton->SetParent ( mtb->GetFrame ( ) );
        this->TransformIconButton->Create ( );
        this->TransformIconButton->SetReliefToFlat ( );
        this->TransformIconButton->SetBorderWidth ( 0 );
        this->TransformIconButton->SetOverReliefToNone ( );
        this->TransformIconButton->SetImageToIcon ( this->SlicerToolbarIcons->GetTransformIcon ( ) );
        this->TransformIconButton->SetBalloonHelpString ( "Transforms");
        mtb->AddWidget ( this->TransformIconButton );

        // fiducial utility icon
        this->FiducialsIconButton->SetParent ( mtb->GetFrame ( ) );
        this->FiducialsIconButton->Create ( );
        this->FiducialsIconButton->SetReliefToFlat ( );
        this->FiducialsIconButton->SetBorderWidth ( 0 );
        this->FiducialsIconButton->SetOverReliefToNone ( );
        this->FiducialsIconButton->SetImageToIcon ( this->SlicerToolbarIcons->GetFiducialsIcon ( ) );
        this->FiducialsIconButton->SetBalloonHelpString ( "Fiducials");
        mtb->AddWidget ( this->FiducialsIconButton );

        // editor module icon
        this->EditorToolboxIconButton->SetParent ( mtb->GetFrame ( ) );
        this->EditorToolboxIconButton->Create ( );
        this->EditorToolboxIconButton->SetReliefToFlat ( );
        this->EditorToolboxIconButton->SetBorderWidth ( 0 );
        this->EditorToolboxIconButton->SetOverReliefToNone ( );
        this->EditorToolboxIconButton->SetImageToIcon ( this->SlicerToolbarIcons->GetEditorToolboxIcon ( ) );
        this->EditorToolboxIconButton->SetBalloonHelpString ( "Editor Toolbox");        
        mtb->AddWidget ( this->EditorToolboxIconButton );
        // editor module icon
        this->EditorIconButton->SetParent ( mtb->GetFrame ( ) );
        this->EditorIconButton->Create ( );
        this->EditorIconButton->SetReliefToFlat ( );
        this->EditorIconButton->SetBorderWidth ( 0 );
        this->EditorIconButton->SetOverReliefToNone ( );
        this->EditorIconButton->SetImageToIcon ( this->SlicerToolbarIcons->GetEditorIcon ( ) );
        this->EditorIconButton->SetBalloonHelpString ( "Editor");        
        mtb->AddWidget ( this->EditorIconButton );

        // color utility icon
        this->ColorIconButton->SetParent ( mtb->GetFrame ( ) );
        this->ColorIconButton->Create ( );
        this->ColorIconButton->SetReliefToFlat ( );
        this->ColorIconButton->SetBorderWidth ( 0 );
        this->ColorIconButton->SetOverReliefToNone ( );
        this->ColorIconButton->SetImageToIcon ( this->SlicerToolbarIcons->GetColorIcon ( ) );
        this->ColorIconButton->SetBalloonHelpString ( "Colors");
        mtb->AddWidget ( this->ColorIconButton );

        // save scene icon
        this->SaveSceneIconButton->SetParent ( ltb->GetFrame ( ));
        this->SaveSceneIconButton->Create ( );
        this->SaveSceneIconButton->SetReliefToFlat ( );
        this->SaveSceneIconButton->SetBorderWidth ( 0 );
        this->SaveSceneIconButton->SetOverReliefToNone ( );
        this->SaveSceneIconButton->SetImageToIcon ( this->SlicerToolbarIcons->GetSaveSceneIcon( ) );
        this->SaveSceneIconButton->SetBalloonHelpString ( "Save a MRML scene to a file.");
        ltb->AddWidget ( this->SaveSceneIconButton );

        // load scene icon
        this->LoadSceneIconButton->SetParent ( ltb->GetFrame ( ) );
        this->LoadSceneIconButton->Create();
        this->LoadSceneIconButton->SetReliefToFlat ( );
        this->LoadSceneIconButton->SetBorderWidth ( 0 );
        this->LoadSceneIconButton->SetOverReliefToNone ( );
        this->LoadSceneIconButton->SetImageToIcon ( this->SlicerToolbarIcons->GetLoadSceneIcon( ) );
        this->LoadSceneIconButton->SetBalloonHelpString ( "Load a MRML scene.");
        ltb->AddWidget ( this->LoadSceneIconButton );

        // conventional view icon
        this->ConventionalViewIconButton->SetParent (vtb->GetFrame ( ) );
        this->ConventionalViewIconButton->Create ( );
        this->ConventionalViewIconButton->SetReliefToFlat ( );
        this->ConventionalViewIconButton->SetBorderWidth ( 0 );
        this->ConventionalViewIconButton->SetOverReliefToNone ( );
        this->ConventionalViewIconButton->SetImageToIcon ( this->SlicerToolbarIcons->GetConventionalViewIcon ( ) );        
        this->ConventionalViewIconButton->SetBalloonHelpString ("Display the 3D viewer over 3 slice windows");
        vtb->AddWidget ( this->ConventionalViewIconButton );
        // 3Dview-only icon
        this->OneUp3DViewIconButton->SetParent ( vtb->GetFrame ( ) );
        this->OneUp3DViewIconButton->Create ( );
        this->OneUp3DViewIconButton->SetReliefToFlat ( );
        this->OneUp3DViewIconButton->SetBorderWidth ( 0 );
        this->OneUp3DViewIconButton->SetOverReliefToNone ( );
        this->OneUp3DViewIconButton->SetImageToIcon ( this->SlicerToolbarIcons->GetOneUp3DViewIcon ( ) );
        this->OneUp3DViewIconButton->SetBalloonHelpString ( "Display the 3D viewer without any slice windows" );
        vtb->AddWidget (this->OneUp3DViewIconButton );

        // Slice view-only icon
        this->OneUpSliceViewIconButton->SetParent ( vtb->GetFrame ( ) );
        this->OneUpSliceViewIconButton->Create ( );
        this->OneUpSliceViewIconButton->SetReliefToFlat ( );
        this->OneUpSliceViewIconButton->SetBorderWidth ( 0 );
        this->OneUpSliceViewIconButton->SetOverReliefToNone ( );
        this->OneUpSliceViewIconButton->SetImageToIcon ( this->SlicerToolbarIcons->GetOneUpSliceViewIcon ( ) );
        this->OneUpSliceViewIconButton->SetBalloonHelpString ( "Display one slice window with no 3D viewer" );
        vtb->AddWidget (this->OneUpSliceViewIconButton );

        // 4 equal windows icon
        this->FourUpViewIconButton->SetParent ( vtb->GetFrame ( ) );
        this->FourUpViewIconButton->Create ( );
        this->FourUpViewIconButton->SetReliefToFlat ( );
        this->FourUpViewIconButton->SetBorderWidth ( 0 );
        this->FourUpViewIconButton->SetOverReliefToNone ( );
        this->FourUpViewIconButton->SetImageToIcon ( this->SlicerToolbarIcons->GetFourUpViewIcon ( ) );
        this->FourUpViewIconButton->SetBalloonHelpString ( "Display the 3D viewer and 3 slice windows in a matrix" );
        vtb->AddWidget ( this->FourUpViewIconButton );

        // tabbed view icon
        this->TabbedViewIconButton->SetParent ( vtb->GetFrame ( ) );
        this->TabbedViewIconButton->Create ( );
        this->TabbedViewIconButton->SetReliefToFlat ( );
        this->TabbedViewIconButton->SetBorderWidth ( 0 );
        this->TabbedViewIconButton->SetOverReliefToNone ( );
        this->TabbedViewIconButton->SetImageToIcon ( this->SlicerToolbarIcons->GetTabbedViewIcon ( ) );
        this->TabbedViewIconButton->SetBalloonHelpString ( "Display a collection of scenes in a notebook" );
        vtb->AddWidget ( this->TabbedViewIconButton );

        // lightbox view icon
        this->LightBoxViewIconButton->SetParent ( vtb->GetFrame ( ));
        this->LightBoxViewIconButton->Create ( );
        this->LightBoxViewIconButton->SetReliefToFlat ( );
        this->LightBoxViewIconButton->SetBorderWidth ( 0 );
        this->LightBoxViewIconButton->SetOverReliefToNone ( );
        this->LightBoxViewIconButton->SetImageToIcon ( this->SlicerToolbarIcons->GetLightBoxViewIcon( ) );
        this->LightBoxViewIconButton->SetBalloonHelpString ( "Display a slice-matrix and no 3D view" );
        vtb->AddWidget ( this->LightBoxViewIconButton );

        tbs->ShowToolbar ( this->GetModulesToolbar ( ));
        tbs->ShowToolbar ( this->GetLoadSaveToolbar ( ));
        tbs->ShowToolbar ( this->GetViewToolbar ( ));
    }

}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::BuildLogoGUIPanel ( )
{
    if ( this->GetApplication( )  != NULL ) {
        vtkSlicerApplication *app = vtkSlicerApplication::SafeDownCast ( this->GetApplication () );
        this->SlicerLogoLabel->SetParent ( this->LogoFrame );
        this->SlicerLogoLabel->Create();
        this->SlicerLogoLabel->SetImageToIcon ( this->SlicerLogoIcons->GetSlicerLogo() );
        this->SlicerLogoLabel->SetBalloonHelpString ("placeholder logo");
        app->Script ( "pack %s -side top -anchor w -padx 2 -pady 0", this->SlicerLogoLabel->GetWidgetName( ) );        
    }
    
}



//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::BuildSlicerControlGUIPanel ( )
{
    const char* mName;
    vtkSlicerModuleGUI *m;
    
    if ( this->GetApplication( )  != NULL ) {
        vtkSlicerApplication *app = vtkSlicerApplication::SafeDownCast( this->GetApplication() );
        
        //--- ALL modules menu button label
        this->ModulesLabel->SetParent ( this->ModuleChooseFrame );
        this->ModulesLabel->Create ( );
        this->ModulesLabel->SetText ( "Modules:");
        this->ModulesLabel->SetAnchorToWest ( );
        this->ModulesLabel->SetWidth ( 7 );

        //--- All modules menu button
        this->ModulesMenuButton->SetParent ( this->ModuleChooseFrame );
        this->ModulesMenuButton->Create ( );
        this->ModulesMenuButton->SetWidth ( 28 );
        this->ModulesMenuButton->IndicatorVisibilityOn ( );
        this->ModulesMenuButton->SetBalloonHelpString ("Select a Slicer module.");
        //--- ALL modules pull-down menu 
        if ( app->GetModuleGUICollection ( ) != NULL ) {
            app->GetModuleGUICollection( )->InitTraversal( );
            m = vtkSlicerModuleGUI::SafeDownCast( app->GetModuleGUICollection( )->GetNextItemAsObject( ));
            while ( m != NULL ) {
                mName = m->GetUIPanel( )->GetName( );
                this->ModulesMenuButton->GetMenu( )->AddRadioButton( mName );
                m = vtkSlicerModuleGUI::SafeDownCast( app->GetModuleGUICollection( )->GetNextItemAsObject( ));
            }
        }
        //--- TODO: make the initial value be module user sets as "home"
        this->ModulesMenuButton->SetValue ("Volumes");
        
        //--- Next and previous module button
        this->ModulesNext->SetParent ( this->ModuleChooseFrame );
        this->ModulesNext->Create ( );
        this->ModulesNext->SetBorderWidth ( 0 );
        this->ModulesNext->SetImageToIcon ( this->SlicerModuleNavigationIcons->GetModuleNextIcon() );
        this->ModulesNext->SetBalloonHelpString ("Navigate to the next module in your use history.");

        this->ModulesPrev->SetParent ( this->ModuleChooseFrame );
        this->ModulesPrev->Create ( );
        this->ModulesPrev->SetBorderWidth ( 0 );
        this->ModulesPrev->SetImageToIcon ( this->SlicerModuleNavigationIcons->GetModulePrevIcon() );
        this->ModulesPrev->SetBalloonHelpString ("Navigate to the previous module in your use history.");
        
        this->ModulesHistory->SetParent ( this->ModuleChooseFrame );
        this->ModulesHistory->Create ( );
        this->ModulesHistory->SetBorderWidth ( 0 );
        this->ModulesHistory->SetImageToIcon ( this->SlicerModuleNavigationIcons->GetModuleHistoryIcon() );
        this->ModulesHistory->SetBalloonHelpString ("Pop up a window showing your module use history.");
        
        //--- pack everything up.
        app->Script ( "pack %s -side left -anchor n -padx 1 -ipadx 1 -pady 3", this->ModulesLabel->GetWidgetName( ) );
        app->Script ( "pack %s -side left -anchor n -padx 1 -ipady 0 -pady 2", this->ModulesMenuButton->GetWidgetName( ) );
        app->Script ( "pack %s -side left -anchor c -padx 2 -pady 2", this->ModulesPrev->GetWidgetName( ) );
        app->Script ( "pack %s -side left -anchor c -padx 2 -pady 2", this->ModulesNext->GetWidgetName( ) );
        app->Script ( "pack %s -side left -anchor c -padx 2 -pady 2", this->ModulesHistory->GetWidgetName( ) );
    }
}



//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::BuildSliceControlGUIPanel ( )
{

    //--- Populate the Slice Control Frame

    if ( this->GetApplication( )  != NULL ) {
        vtkSlicerApplication *app = vtkSlicerApplication::SafeDownCast( this->GetApplication() );
        this->SliceControlFrame->SetReliefToGroove();
        
        //--- create frames
        vtkKWFrame *f1 = vtkKWFrame::New ( );
        f1->SetParent ( this->SliceControlFrame );
        f1->Create ( );
        vtkKWFrame *f2 = vtkKWFrame::New ( );
        f2->SetParent ( this->SliceControlFrame );
        f2->Create ( );
        vtkKWFrame *f3 = vtkKWFrame::New ( );
        f3->SetParent ( this->SliceControlFrame );
        f3->Create ( );
        
        //--- pack everything up: buttons, labels, scales
        app->Script ( "pack %s -side left -anchor n -padx 0 -pady 5", f1->GetWidgetName( ) );
        app->Script ( "pack %s -side left -anchor n -padx 0 -pady 5", f2->GetWidgetName( ) );
        app->Script ( "pack %s -side left -anchor n -padx 0 -pady 5", f3->GetWidgetName( ) );

        //--- make buttons for toggling FG/BG and annotations
        this->ToggleFgBgButton->SetParent ( f1 );
        this->ToggleFgBgButton->Create ( );
        this->ToggleFgBgButton->SetWidth ( 16 );
        this->ToggleFgBgButton->SetText ( "Toggle FG/BG" );
        this->ToggleAnnotationButton->SetParent ( f1 );
        this->ToggleAnnotationButton->Create ( );
        this->ToggleAnnotationButton->SetWidth ( 16 );
        this->ToggleAnnotationButton->SetText ( "Toggle Annotation" );
    
        app->Script ( "pack %s -side top -anchor w -padx 1 -pady 1", this->ToggleFgBgButton->GetWidgetName( ) );
        app->Script ( "pack %s -side top -anchor w -padx 1 -pady 1", this->ToggleAnnotationButton->GetWidgetName( ) );

        //--- make labels (can't reposition the Scale's labels, so
        //--- supressing those and using a new set.)
        vtkKWLabel *fadeLabel = vtkKWLabel::New ( );
        vtkKWLabel *opacityLabel = vtkKWLabel::New ( );
        fadeLabel->SetParent ( f2 );
        fadeLabel->Create ( );
        fadeLabel->SetWidth ( 14 );
        fadeLabel->SetAnchorToEast ( );
        fadeLabel->SetText ( "Fade (FG/BG):");
        opacityLabel->SetParent ( f2 );
        opacityLabel->Create ( );
        opacityLabel->SetWidth ( 14 );
        opacityLabel->SetAnchorToEast ( );
        opacityLabel->SetText ( "Opacity (0,1):");
        app->Script ( "pack %s -side top -anchor e -padx 1 -pady 1", fadeLabel->GetWidgetName( ) );
        app->Script ( "pack %s -side top -anchor e -padx 1 -pady 2", opacityLabel->GetWidgetName( ) );
        
        //--- make scales for sliding slice visibility in the SliceViewers
        //--- and for sliding slice opacity in the 3D Viewer.
        this->SliceFadeScale->SetParent ( f3 );
        this->SliceFadeScale->Create ( );
        this->SliceFadeScale->SetRange (0.0, 1.0);
        this->SliceFadeScale->SetResolution ( 0.01 );
        this->SliceFadeScale->SetValue ( 0.0 );
        this->SliceFadeScale->SetLength ( 120 );
        this->SliceFadeScale->SetOrientationToHorizontal ( );
        this->SliceFadeScale->ValueVisibilityOff ( );
        this->SliceFadeScale->SetBalloonHelpString ( "Scale fades between FG and BG Slice Layers" );

        this->SliceOpacityScale->SetParent ( f3 );
        this->SliceOpacityScale->Create ( );
        this->SliceOpacityScale->SetRange ( 0.0, 1.0 );
        this->SliceOpacityScale->SetResolution ( 0.01 );
        this->SliceOpacityScale->SetValue ( 0.0 );
        this->SliceOpacityScale->SetLength ( 120 );
        this->SliceOpacityScale->SetOrientationToHorizontal ( );
        this->SliceOpacityScale->ValueVisibilityOff ( );
        this->SliceOpacityScale->SetBalloonHelpString ( "Scale sets the opacity of the slice plane in the 3D Viewer" );

        app->Script ( "pack %s -side top -anchor w -padx 0 -pady 1", this->SliceFadeScale->GetWidgetName( ) );
        app->Script ( "pack %s -side top -anchor w -padx 0 -pady 0", this->SliceOpacityScale->GetWidgetName( ) );

        fadeLabel->Delete ( );
        opacityLabel->Delete ( );
        f1->Delete ( );
        f2->Delete ( );
        f3->Delete ( );
    }
}


//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::AssignViewControlIcons ( )
{
        //--- assign image data to each label
        this->RotateAroundAIconButton->SetImageToIcon ( this->SlicerViewControlIcons->GetRotateAroundAIconLO() );
        this->RotateAroundPIconButton->SetImageToIcon ( this->SlicerViewControlIcons->GetRotateAroundPIconLO( ) );
        this->RotateAroundRIconButton->SetImageToIcon ( this->SlicerViewControlIcons->GetRotateAroundRIconLO ( ));
        this->RotateAroundLIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundLIconLO ( ));        
        this->RotateAroundSIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundSIconLO ( ));
        this->RotateAroundIIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundIIconLO ( ) );
        this->RotateAroundMiddleIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundMiddleIcon ( ) );
        this->RotateAroundTopCornerIconButton->SetImageToIcon ( this->SlicerViewControlIcons->GetRotateAroundTopCornerIcon ( ));
        this->RotateAroundBottomCornerIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundBottomCornerIcon ( ));
        this->LookFromAIconButton->SetImageToIcon ( this->SlicerViewControlIcons->GetLookFromAIconLO() );
        this->LookFromPIconButton->SetImageToIcon ( this->SlicerViewControlIcons->GetLookFromPIconLO( ) );
        this->LookFromRIconButton->SetImageToIcon ( this->SlicerViewControlIcons->GetLookFromRIconLO ( ));
        this->LookFromLIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetLookFromLIconLO ( ));        
        this->LookFromSIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetLookFromSIconLO ( ));
        this->LookFromIIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetLookFromIIconLO ( ) );
        this->LookFromMiddleIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetLookFromMiddleIcon ( ) );
        this->LookFromTopCornerIconButton->SetImageToIcon ( this->SlicerViewControlIcons->GetLookFromTopCornerIcon ( ));
        this->LookFromBottomCornerIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetLookFromBottomCornerIcon ( ));
}



//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::MakeViewControlRolloverBehavior ( )
{

        //--- configure and bind for rollover interaction
        this->RotateAroundAIconButton->SetBorderWidth (0);
        this->RotateAroundAIconButton->SetBinding ( "<Enter>",  this, "EnterRotateAroundACallback");
        this->RotateAroundAIconButton->SetBinding ( "<Leave>",  this, "LeaveRotateAroundACallback");
        this->Script ( "%s ListMethods", this->GetTclName() );

        this->RotateAroundPIconButton->SetBorderWidth (0);
        this->RotateAroundPIconButton->SetBinding ( "<Enter>", this, "EnterRotateAroundPCallback");
        this->RotateAroundPIconButton->SetBinding ( "<Leave>", this, "LeaveRotateAroundPCallback");

        this->RotateAroundRIconButton->SetBorderWidth (0);
        this->RotateAroundRIconButton->SetBinding ( "<Enter>", this, "EnterRotateAroundRCallback");
        this->RotateAroundRIconButton->SetBinding ( "<Leave>", this, "LeaveRotateAroundRCallback");

        this->RotateAroundLIconButton->SetBorderWidth (0);
        this->RotateAroundLIconButton->SetBinding ( "<Enter>", this, "EnterRotateAroundLCallback");
        this->RotateAroundLIconButton->SetBinding ( "<Leave>", this, "LeaveRotateAroundLCallback");

        this->RotateAroundSIconButton->SetBorderWidth (0);
        this->RotateAroundSIconButton->SetBinding ( "<Enter>", this, "EnterRotateAroundSCallback");
        this->RotateAroundSIconButton->SetBinding ( "<Leave>", this, "LeaveRotateAroundSCallback");
        
        this->RotateAroundIIconButton->SetBorderWidth (0);
        this->RotateAroundIIconButton->SetBinding ( "<Enter>", this, "EnterRotateAroundICallback");
        this->RotateAroundIIconButton->SetBinding ( "<Leave>", this, "LeaveRotateAroundICallback");
        
        this->RotateAroundMiddleIconButton->SetBorderWidth (0);
        this->RotateAroundTopCornerIconButton->SetBorderWidth (0);
        this->RotateAroundBottomCornerIconButton->SetBorderWidth (0);

        this->LookFromAIconButton->SetBorderWidth (0);
        this->LookFromAIconButton->SetBinding ( "<Enter>", this, "EnterLookFromACallback");
        this->LookFromAIconButton->SetBinding ( "<Leave>", this, "LeaveLookFromACallback");
        
        this->LookFromPIconButton->SetBorderWidth (0);
        this->LookFromPIconButton->SetBinding ( "<Enter>", this, "EnterLookFromPCallback");
        this->LookFromPIconButton->SetBinding ( "<Leave>", this, "LeaveLookFromPCallback");
        
        this->LookFromRIconButton->SetBorderWidth (0);
        this->LookFromRIconButton->SetBinding ( "<Enter>", this, "EnterLookFromRCallback");
        this->LookFromRIconButton->SetBinding ( "<Leave>", this, "LeaveLookFromRCallback");
        
        this->LookFromLIconButton->SetBorderWidth (0);
        this->LookFromLIconButton->SetBinding ( "<Enter>", this, "EnterLookFromLCallback");
        this->LookFromLIconButton->SetBinding ( "<Leave>", this, "LeaveLookFromLCallback");
        
        this->LookFromSIconButton->SetBorderWidth (0);
        this->LookFromSIconButton->SetBinding ( "<Enter>", this, "EnterLookFromSCallback");
        this->LookFromSIconButton->SetBinding ( "<Leave>", this, "LeaveLookFromSCallback");
        
        this->LookFromIIconButton->SetBorderWidth (0);
        this->LookFromIIconButton->SetBinding ( "<Enter>", this, "EnterLookFromICallback");
        this->LookFromIIconButton->SetBinding ( "<Leave>", this, "LeaveLookFromICallback");
        
        this->LookFromMiddleIconButton->SetBorderWidth (0);
        this->LookFromTopCornerIconButton->SetBorderWidth (0);
        this->LookFromBottomCornerIconButton->SetBorderWidth (0);
}


//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::BuildViewControlGUIPanel ( )
{
    if ( this->GetApplication( )  != NULL ) {

        vtkSlicerApplication *app = vtkSlicerApplication::SafeDownCast( this->GetApplication() );
        vtkSlicerGUILayout *layout = app->GetMainLayout ( );

        vtkKWFrame *f1 = vtkKWFrame::New ( );
        vtkKWFrame *f1a = vtkKWFrame::New ( );    
        vtkKWFrame *f1b = vtkKWFrame::New ( );    
        vtkKWFrame *f2 = vtkKWFrame::New ( );

        // divide the GUI panel into two frames of identical wid.
        int wid = layout->GetDefaultGUIPanelWidth() ;
        int buf = 4;
        int thirdwid = wid/3 - buf;
        
        // create frames and set their widths.
        f1->SetParent ( this->ViewControlFrame );
        f1->Create ( );
        f1->SetWidth ( thirdwid );
        f1->SetHeight (layout->GetDefaultViewControlFrameHeight ( ) );
        f1->SetReliefToGroove();

        f2->SetParent ( this->ViewControlFrame );
        f2->Create ( );
        f2->SetWidth ( 2 * thirdwid );
        f2->SetHeight (layout->GetDefaultViewControlFrameHeight ( ) );

        f1a->SetParent ( f1 );
        f1a->Create ( );
        f1a->SetWidth (thirdwid );

        f1b->SetParent ( f1 );
        f1b->Create ( );
        f1b->SetWidth ( thirdwid );
        
        //--- create rotate-around and look-from image mosaics from vtkKWLabels
        this->RotateAroundAIconButton->SetParent ( f1b );
        this->RotateAroundAIconButton->Create ( );
        this->RotateAroundAIconButton->SetBalloonHelpString ("Rotate camera in 3D view around A-P axis.");
        this->RotateAroundPIconButton->SetParent ( f1b );
        this->RotateAroundPIconButton->Create ( );
        this->RotateAroundPIconButton->SetBalloonHelpString ("Rotate camera in 3D view around A-P axis.");
        this->RotateAroundRIconButton->SetParent ( f1b );
        this->RotateAroundRIconButton->Create ( );
        this->RotateAroundRIconButton->SetBalloonHelpString ("Rotate camera in 3D view around R-L axis.");
        this->RotateAroundLIconButton->SetParent ( f1b );
        this->RotateAroundLIconButton->Create ( );
        this->RotateAroundLIconButton->SetBalloonHelpString ("Rotate camera in 3D view around R-L axis.");
        this->RotateAroundSIconButton->SetParent ( f1b );
        this->RotateAroundSIconButton->Create ( );
        this->RotateAroundSIconButton->SetBalloonHelpString ("Rotate camera in 3D view around S-I axis.");
        this->RotateAroundIIconButton->SetParent ( f1b );
        this->RotateAroundIIconButton->Create ( );
        this->RotateAroundIIconButton->SetBalloonHelpString ("Rotate camera in 3D view around S-I axis.");
        this->RotateAroundMiddleIconButton->SetParent ( f1b );
        this->RotateAroundMiddleIconButton->Create ( );
        this->RotateAroundTopCornerIconButton->SetParent ( f1b );
        this->RotateAroundTopCornerIconButton->Create ( );
        this->RotateAroundBottomCornerIconButton->SetParent ( f1b );
        this->RotateAroundBottomCornerIconButton->Create ( );
        this->LookFromAIconButton->SetParent ( f1b );
        this->LookFromAIconButton->Create ( );
        this->LookFromAIconButton->SetBalloonHelpString ("Position 3D view camera down the A-axis looking toward center.");
        this->LookFromPIconButton->SetParent ( f1b );
        this->LookFromPIconButton->Create ( );
        this->LookFromPIconButton->SetBalloonHelpString ("Position 3D view camera down the P-axis looking toward center.");
        this->LookFromRIconButton->SetParent ( f1b );
        this->LookFromRIconButton->Create ( );
        this->LookFromRIconButton->SetBalloonHelpString ("Position 3D view camera down the R-axis looking toward center.");
        this->LookFromLIconButton->SetParent ( f1b );
        this->LookFromLIconButton->Create ( );
        this->LookFromLIconButton->SetBalloonHelpString ("Position 3D view camera down the L-axis looking toward center.");
        this->LookFromSIconButton->SetParent ( f1b );
        this->LookFromSIconButton->Create ( );
        this->LookFromSIconButton->SetBalloonHelpString ("Position 3D view camera down the S-axis looking toward center.");
        this->LookFromIIconButton->SetParent ( f1b );
        this->LookFromIIconButton->Create ( );
        this->LookFromIIconButton->SetBalloonHelpString ("Position 3D view camera down the I-axis looking toward center.");
        this->LookFromMiddleIconButton->SetParent ( f1b );
        this->LookFromMiddleIconButton->Create ( );
        this->LookFromTopCornerIconButton->SetParent ( f1b );
        this->LookFromTopCornerIconButton->Create ( );
        this->LookFromBottomCornerIconButton->SetParent ( f1b );
        this->LookFromBottomCornerIconButton->Create ( );

        this->AssignViewControlIcons ( );
        this->MakeViewControlRolloverBehavior ( );
        
        //--- create the nav/zoom widgets
        this->NavZoomInIconButton->SetParent ( f2 );
        this->NavZoomInIconButton->Create ( );
        this->NavZoomInIconButton->SetReliefToFlat ( );        
        this->NavZoomOutIconButton->SetParent ( f2 );        
        this->NavZoomOutIconButton->Create ( );
        this->NavZoomOutIconButton->SetReliefToFlat ( );
        this->NavZoomScale->SetParent ( f2 );
        this->NavZoomScale->Create ( );
        this->NavZoomScale->SetRange (0.0, 1.0);
        this->NavZoomScale->SetResolution ( 0.01 );
        this->NavZoomScale->SetBorderWidth ( 1 );
        this->NavZoomScale->SetValue ( 0.0 );
        // make scale long enough to fill the frame,
        // leaving room for the zoomin, zoomout buttons.
        this->NavZoomScale->SetLength ( 120 );
        this->NavZoomScale->SetOrientationToHorizontal ( );
        this->NavZoomScale->ValueVisibilityOff ( );
        this->NavZoomScale->SetBalloonHelpString ( "Use scale to zoom the navigation window in/out" );
        //--- assign image data to the zoom buttons
        this->NavZoomInIconButton->SetImageToIcon ( this->SlicerViewControlIcons->GetNavZoomInIcon() );
        this->NavZoomOutIconButton->SetImageToIcon ( this->SlicerViewControlIcons->GetNavZoomOutIcon() );

        // temporary thing until navzoom window is built.
        vtkKWLabel *tmpNavZoom = vtkKWLabel::New ( );
        tmpNavZoom->SetParent (f2);
        tmpNavZoom->SetWidth ( 20);
        tmpNavZoom->SetHeight (10 );
        tmpNavZoom->Create();
        tmpNavZoom->SetText ( "3DNav / SliceZoom" );
        tmpNavZoom->SetBackgroundColor ( app->GetSlicerTheme()->GetSlicerColors()->ViewerBlue );

        //--- other camera control widgets
        this->SpinButton->SetParent ( f1a);
        this->SpinButton->Create ( );
        this->SpinButton->SetText ( "Spin" );

        this->RockButton->SetParent ( f1a );
        this->RockButton->Create ( );
        this->RockButton->SetText ( "Rock" );

        this->OrthoButton->SetParent ( f1a );
        this->OrthoButton->Create ( );
        this->OrthoButton->SetText ( "Ortho" );

        this->CenterButton->SetParent ( f1a );
        this->CenterButton->Create ( );
        this->CenterButton->SetText ( "Center");

        this->SelectButton->SetParent ( f1a );
        this->SelectButton->Create ( );
        this->SelectButton->SetValue ( "Select");

        this->FOVEntry->SetParent ( f1a );
        this->FOVEntry->Create ( );
        this->FOVEntry->SetLabelText ( "FOV: ");
        this->FOVEntry->GetWidget()->SetWidth (4);

        
        this->Script ( "pack %s -side left -anchor n -padx 2 -pady 2 -expand n", f1->GetWidgetName ( ) );
        this->Script ( "pack %s -side left -anchor n -fill x -padx 5 -pady 2 -expand n", f2->GetWidgetName( ) );    

        this->Script ( "pack %s -side top -padx 0 -pady 0 -anchor n -expand n ", f1a->GetWidgetName( ) );
        this->Script ( "pack %s -side top -padx 0 -pady 0 -anchor n -expand n ", f1b->GetWidgetName() );
        
        this->Script ("grid %s -row 0 -column 0 -sticky w -padx 3 -pady 2", this->SpinButton->GetWidgetName ( ) );
        this->Script ("grid %s -row 1 -column 0 -sticky w -padx 3 -pady 2", this->RockButton->GetWidgetName ( ) );
        this->Script ("grid %s -row 2 -column 0 -sticky w -padx 3 -pady 2", this->OrthoButton->GetWidgetName ( ) );
        this->Script ("grid %s -row 0 -column 1 -sticky ew -padx 0 -pady 2", this->CenterButton->GetWidgetName ( ) );
        this->Script ("grid %s -row 1 -column 1 -sticky ew -padx 0 -pady 2", this->SelectButton->GetWidgetName ( ) );
        this->Script ("grid %s -row 2 -column 1 -sticky ew -padx 0 -pady 2", this->FOVEntry->GetWidgetName ( ) );
        
        this->Script ("grid %s -row 0 -column 0 -sticky w -padx 0 -pady 0 -ipadx 0 -ipady 0", this->RotateAroundPIconButton->GetWidgetName ( ));
        this->Script ("grid %s -row 0 -column 1 -sticky w -padx 0 -pady 0 -ipadx 0 -ipady 0", this->RotateAroundSIconButton->GetWidgetName ( ));
        this->Script ("grid %s -row 0 -column 2 -sticky w -padx 0 -pady 0 -ipadx 0 -ipady 0",this->RotateAroundTopCornerIconButton->GetWidgetName ( ) );
        this->Script ("grid %s -row 0 -column 3 -sticky w -padx 0 -pady 0 -ipadx 0 -ipady 0", this->LookFromPIconButton->GetWidgetName ( ));
        this->Script ("grid %s -row 0 -column 4 -sticky w -padx 0 -pady 0 -ipadx 0 -ipady 0",  this->LookFromSIconButton->GetWidgetName ( ));
        this->Script ("grid %s -row 0 -column 5 -sticky w -padx 0 -pady 0 -ipadx 0 -ipady 0", this->LookFromTopCornerIconButton->GetWidgetName ( ));        
                      
        this->Script ("grid %s -row 1 -column 0 -sticky w -padx 0 -pady 0 -ipadx 0 -ipady 0", this->RotateAroundRIconButton->GetWidgetName ( ) );
        this->Script ("grid %s -row 1 -column 1 -sticky w -padx 0 -pady 0 -ipadx 0 -ipady 0", this->RotateAroundMiddleIconButton->GetWidgetName ( ));
        this->Script ("grid %s -row 1 -column 2 -sticky w -padx 0 -pady 0 -ipadx 0 -ipady 0", this->RotateAroundLIconButton->GetWidgetName ( ) );
        this->Script ("grid %s -row 1 -column 3 -sticky w -padx 0 -pady 0 -ipadx 0 -ipady 0", this->LookFromRIconButton->GetWidgetName ( ) );
        this->Script ("grid %s -row 1 -column 4 -sticky w -padx 0 -pady 0 -ipadx 0 -ipady 0", this->LookFromMiddleIconButton->GetWidgetName ( ) );
        this->Script ("grid %s -row 1 -column 5 -sticky w -padx 0 -pady 0 -ipadx 0 -ipady 0", this->LookFromLIconButton->GetWidgetName ( ));        

        this->Script ("grid %s -row 2 -column 0 -sticky w -padx 0 -pady 0 -ipadx 0 -ipady 0", this->RotateAroundBottomCornerIconButton->GetWidgetName ( ));
        this->Script ("grid %s -row 2 -column 1 -sticky w -padx 0 -pady 0 -ipadx 0 -ipady 0", this->RotateAroundIIconButton->GetWidgetName ( ) );
        this->Script ("grid %s -row 2 -column 2 -sticky w -padx 0 -pady 0 -ipadx 0 -ipady 0", this->RotateAroundAIconButton->GetWidgetName ( ) );
        this->Script ("grid %s -row 2 -column 3 -sticky w -padx 0 -pady 0 -ipadx 0 -ipady 0", this->LookFromBottomCornerIconButton->GetWidgetName ( ) );
        this->Script ("grid %s -row 2 -column 4 -sticky w -padx 0 -pady 0 -ipadx 0 -ipady 0", this->LookFromIIconButton->GetWidgetName ( ) );
        this->Script ("grid %s -row 2 -column 5 -sticky w -padx 0 -pady 0 -ipadx 0 -ipady 0", this->LookFromAIconButton->GetWidgetName ( ));        

        this->Script ( "grid %s -row 1 -column 0 -padx 0 -pady 0 -sticky ew", this->NavZoomOutIconButton->GetWidgetName() );
        this->Script ( "grid %s -row 1 -column 1 -padx 0 -pady 0 -sticky ew", this->NavZoomScale->GetWidgetName() );
        this->Script ( "grid %s -row 1 -column 2 -padx 0 -pady 0 -sticky ew", this->NavZoomInIconButton->GetWidgetName() );
        this->Script ("grid %s -row 0 -columnspan 3 -ipadx 40 -ipady 60 -padx 0 -pady 0 -sticky nsew", tmpNavZoom->GetWidgetName ( ) );
        
        f1a->Delete();
        f1b->Delete();
        f1->Delete();
        f2->Delete();
    }
}



//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::ConfigureMainSlicerWindow ( )
{

    if ( this->GetApplication() != NULL ) {
        vtkSlicerApplication *app = vtkSlicerApplication::SafeDownCast(this->GetApplication() );
        if ( this->MainSlicerWin != NULL ) {
            this->MainSlicerWin->MainPanelVisibilityOn ();
            this->MainSlicerWin->SecondaryPanelVisibilityOn ();
            this->MainSlicerWin->SetSize ( app->GetMainLayout()->GetDefaultSlicerWindowWidth ( ),
                           app->GetMainLayout()->GetDefaultSlicerWindowHeight () );
            // Configure the minimum width of Slicer's GUI panel.
            // Panel can be expanded and collapsed entirely, but
            // can't be resized by hand to a value smaller than what's set.
            this->MainSlicerWin->GetMainSplitFrame()->SetFrame1Size (app->GetMainLayout()->GetDefaultGUIPanelWidth() );
            this->MainSlicerWin->GetMainSplitFrame()->SetFrame1MinimumSize (app->GetMainLayout()->GetDefaultGUIPanelWidth ( ) );
        }
    }

}


//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::ConfigureMainViewerPanel ( )
{
    if ( this->GetApplication() != NULL ) {
        // pointers for convenience
        vtkSlicerApplication *app = vtkSlicerApplication::SafeDownCast(this->GetApplication() );


        if ( this->MainSlicerWin != NULL ) {
            this->MainSlicerWin->GetViewPanelFrame()->SetWidth ( app->GetMainLayout()->GetDefaultMainViewerWidth() );
        }
    }

}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::ConfigureSliceViewersPanel ( )
{
    if ( this->GetApplication() != NULL ) {
        // pointers for convenience
        vtkSlicerApplication *app = vtkSlicerApplication::SafeDownCast( this->GetApplication() );

        this->MainSlicerWin->GetSecondarySplitFrame()->SetFrame2Size (app->GetMainLayout()->GetDefaultSliceGUIFrameWidth ( ) );
        this->MainSlicerWin->GetSecondarySplitFrame()->SetFrame2MinimumSize (app->GetMainLayout()->GetDefaultSliceGUIFrameWidth ( ) );
        
        if ( this->MainSlicerWin != NULL ) {
            this->MainSlicerWin->GetSecondaryPanelFrame()->SetWidth ( 3 * app->GetMainLayout()->GetDefaultSliceGUIFrameWidth () );
            this->MainSlicerWin->GetSecondaryPanelFrame()->SetHeight ( app->GetMainLayout()->GetDefaultSliceGUIFrameHeight () );

            // Parent and configure Slice0 frame
            this->DefaultSlice0Frame->SetParent ( this->MainSlicerWin->GetSecondaryPanelFrame ( ) );
            this->DefaultSlice0Frame->Create ( );

            // Parent and configure Slice1 frame
            this->DefaultSlice1Frame->SetParent ( this->MainSlicerWin->GetSecondaryPanelFrame ( ) );
            this->DefaultSlice1Frame->Create ( );

            // Parent and configure Slice2 frame
            this->DefaultSlice2Frame->SetParent ( this->MainSlicerWin->GetSecondaryPanelFrame ( ) );
            this->DefaultSlice2Frame->Create ( );
            
            // pack them.
            app->Script ("pack %s -side left  -expand y -fill both -padx 0 -pady 0", this->DefaultSlice0Frame->GetWidgetName( ) );
            app->Script ("pack %s -side left  -expand y -fill both -padx 0 -pady 0", this->DefaultSlice1Frame->GetWidgetName( ) );
            app->Script ("pack %s -side left  -expand y -fill both -padx 0 -pady 0", this->DefaultSlice2Frame->GetWidgetName( ) );


        }
    }

}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::ConfigureGUIPanel ( )
{

    if ( this->GetApplication() != NULL ) {
        // pointers for convenience
        vtkSlicerApplication *app = vtkSlicerApplication::SafeDownCast(this->GetApplication() );

        if ( this->MainSlicerWin != NULL ) {

            this->MainSlicerWin->GetMainPanelFrame()->SetWidth ( app->GetMainLayout()->GetDefaultGUIPanelWidth() );
            this->MainSlicerWin->GetMainPanelFrame()->SetHeight ( app->GetMainLayout()->GetDefaultGUIPanelHeight() );
            this->MainSlicerWin->GetMainPanelFrame()->SetReliefToSunken();

            this->LogoFrame->SetParent ( this->MainSlicerWin->GetMainPanelFrame ( ) );
            this->LogoFrame->Create( );
            this->LogoFrame->SetHeight ( app->GetMainLayout()->GetDefaultLogoFrameHeight ( ) );

            this->ModuleChooseFrame->SetParent ( this->MainSlicerWin->GetMainPanelFrame ( ) );
            this->ModuleChooseFrame->Create( );
            this->ModuleChooseFrame->SetHeight ( app->GetMainLayout()->GetDefaultModuleChooseFrameHeight ( ) );

            // pack logo and slicer control frames
            app->Script ( "pack %s -side top -fill x -padx 1 -pady 1", this->LogoFrame->GetWidgetName() );
            app->Script ( "pack %s -side top -fill x -padx 1 -pady 10", this->ModuleChooseFrame->GetWidgetName() );

            this->SliceControlFrame->SetParent ( this->MainSlicerWin->GetMainPanelFrame ( ) );
            this->SliceControlFrame->Create( );
            this->SliceControlFrame->SetHeight ( app->GetMainLayout()->GetDefaultSliceControlFrameHeight ( ) );
            
            this->ViewControlFrame->SetParent ( this->MainSlicerWin->GetMainPanelFrame ( ) );
            this->ViewControlFrame->Create( );
            this->ViewControlFrame->SetHeight ( app->GetMainLayout()->GetDefaultViewControlFrameHeight ( ) );
            
            // pack slice and view control frames
            app->Script ( "pack %s -side bottom -fill x -padx 1 -pady 10", this->ViewControlFrame->GetWidgetName() );
            app->Script ( "pack %s -side bottom -fill x -padx 1 -pady 10", this->SliceControlFrame->GetWidgetName() );

        }
    }

}




//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::EnterRotateAroundACallback ( ) {
    this->RotateAroundPIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundPIconHI() );
    this->RotateAroundAIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundAIconHI() );
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::LeaveRotateAroundACallback ( ) {
    this->RotateAroundPIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundPIconLO() );
    this->RotateAroundAIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundAIconLO() );
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::EnterRotateAroundPCallback ( ) {
    this->RotateAroundPIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundPIconHI() );
    this->RotateAroundAIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundAIconHI() );
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::LeaveRotateAroundPCallback ( ) {
    this->RotateAroundPIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundPIconLO() );
    this->RotateAroundAIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundAIconLO() );
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::EnterRotateAroundRCallback ( ) {
    this->RotateAroundRIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundRIconHI() );
    this->RotateAroundLIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundLIconHI() );
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::LeaveRotateAroundRCallback ( ) {
    this->RotateAroundRIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundRIconLO() );
    this->RotateAroundLIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundLIconLO() );
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::EnterRotateAroundLCallback ( ) {
    this->RotateAroundRIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundRIconHI() );
    this->RotateAroundLIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundLIconHI() );
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::LeaveRotateAroundLCallback ( ) {
    this->RotateAroundRIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundRIconLO() );
    this->RotateAroundLIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundLIconLO() );
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::EnterRotateAroundSCallback ( ) {
    this->RotateAroundSIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundSIconHI() );
    this->RotateAroundIIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundIIconHI() );
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::LeaveRotateAroundSCallback ( ) {
    this->RotateAroundSIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundSIconLO() );
    this->RotateAroundIIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundIIconLO() );
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::EnterRotateAroundICallback ( ) {
    this->RotateAroundIIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundIIconHI() );
    this->RotateAroundSIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundSIconHI() );
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::LeaveRotateAroundICallback ( ) {
    this->RotateAroundIIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundIIconLO() );
    this->RotateAroundSIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetRotateAroundSIconLO() );

}



//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::EnterLookFromACallback ( ) {
    this->LookFromAIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetLookFromAIconHI() );
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::LeaveLookFromACallback ( ) {
    this->LookFromAIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetLookFromAIconLO() );
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::EnterLookFromPCallback ( ) {
    this->LookFromPIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetLookFromPIconHI() );
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::LeaveLookFromPCallback ( ) {
    this->LookFromPIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetLookFromPIconLO() );
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::EnterLookFromRCallback ( ) {
    this->LookFromRIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetLookFromRIconHI() );
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::LeaveLookFromRCallback ( ) {
    this->LookFromRIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetLookFromRIconLO() );
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::EnterLookFromLCallback ( ) {
    this->LookFromLIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetLookFromLIconHI() );
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::LeaveLookFromLCallback ( ) {
    this->LookFromLIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetLookFromLIconLO() );
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::EnterLookFromSCallback ( ) {
    this->LookFromSIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetLookFromSIconHI() );
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::LeaveLookFromSCallback ( ) {
    this->LookFromSIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetLookFromSIconLO() );
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::EnterLookFromICallback ( ) {
    this->LookFromIIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetLookFromIIconHI() );
}

//---------------------------------------------------------------------------
void vtkSlicerApplicationGUI::LeaveLookFromICallback ( ) {
    this->LookFromIIconButton->SetImageToIcon (this->SlicerViewControlIcons->GetLookFromIIconLO() );
}
