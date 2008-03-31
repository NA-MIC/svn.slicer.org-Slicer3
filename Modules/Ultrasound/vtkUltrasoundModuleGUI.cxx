#include "vtkUltrasoundModuleGUI.h"
#include "vtkObjectFactory.h"

#include "vtkCallbackCommand.h"

#include "vtkKWEvent.h"
#include "vtkKWCheckButton.h"
#include "vtkKWMenuButton.h"
#include "vtkKWMenu.h"
#include "vtkKWScaleWithLabel.h"
#include "vtkKWLabel.h"

#include "vtkSlicerApplication.h"
#include "vtkImageData.h"
#include "vtkMRMLVolumeNode.h"
#include "vtkMRMLScalarVolumeDisplayNode.h"
#include "vtkMRMLScalarVolumeNode.h"
#include "vtkSlicerColorLogic.h"
#include "vtkMRMLScene.h"

vtkCxxRevisionMacro(vtkUltrasoundModuleGUI, "$ Revision 1.0$");
vtkStandardNewMacro(vtkUltrasoundModuleGUI);


vtkUltrasoundModuleGUI::vtkUltrasoundModuleGUI()
{
    // Initialize Values
    this->cb_Enabled = NULL;
    this->sc_RefreshRate = NULL;
    this->cb_Scanning = NULL;
    this->Logic = NULL;
    this->Node = NULL;

}

vtkUltrasoundModuleGUI::~vtkUltrasoundModuleGUI()
{
    if (this->cb_Enabled != NULL)
    this->cb_Enabled->Delete();
    if (this->sc_RefreshRate != NULL)
        this->sc_RefreshRate->Delete();

}




// Description:
// Process events generated by Logic
void vtkUltrasoundModuleGUI::ProcessLogicEvents ( vtkObject *caller, unsigned long event,
                                                 void *callData )
{
}

/// GUI part
void vtkUltrasoundModuleGUI::BuildGUI ( )
{
    vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();
    this->GetUIPanel()->AddPage("Ultrasound", "Ultrasound",NULL);

    // Define your help text and build the help frame here.
    const char *help = "Ultrasound. This is currently a prototype and will be under active development throughout 3DSlicer's Beta release.";
    const char *about = "This work is supported by NA-MIC, NAC, BIRN, NCIGT, and the Slicer Community. See http://www.slicer.org for details.";
    vtkKWWidget *page = this->UIPanel->GetPageWidget ( "Ultrasound" );
    this->BuildHelpAndAboutFrame ( page, help, about );
    //
    //Ultrasound 
    //
    vtkSlicerModuleCollapsibleFrame *usFrame = vtkSlicerModuleCollapsibleFrame::New ( );
    usFrame->SetParent (page);
    usFrame->Create();
    usFrame->ExpandFrame();
    usFrame->SetLabelText("Ultrasound");
    app->Script ( "pack %s -side top -anchor nw -fill x -padx 2 -pady 2",
        usFrame->GetWidgetName());


    this->cb_Enabled = vtkKWCheckButton::New();
    this->cb_Enabled->SetParent(usFrame->GetFrame());
    this->cb_Enabled->Create();
    this->cb_Enabled->SetText("Enabled");
    app->Script("pack %s -side top -anchor nw -fill x -padx 2 -pady 2",
        this->cb_Enabled->GetWidgetName());

    this->sc_RefreshRate = vtkKWScaleWithLabel::New();
    this->sc_RefreshRate->SetParent(usFrame->GetFrame());
    this->sc_RefreshRate->Create();

    this->sc_RefreshRate->GetLabel()->SetText("Refresh Rate");
    app->Script("pack %s -side top -anchor nw -fill x -padx 2 -pady 2",
        this->sc_RefreshRate->GetWidgetName());

}

// This method releases references and key-bindings,
// and optionally removes observers.
void vtkUltrasoundModuleGUI::TearDownGUI ( )
{
}


// Description:
// Methods for adding module-specific key bindings and
// removing them.
void vtkUltrasoundModuleGUI::CreateModuleEventBindings ( )
{
}

void vtkUltrasoundModuleGUI::ReleaseModuleEventBindings ( )
{
}


// Description:
// Add obsereves to GUI widgets
void vtkUltrasoundModuleGUI::AddGUIObservers ( )
{
}


// Description:
// Remove obsereves to GUI widgets
void vtkUltrasoundModuleGUI::RemoveGUIObservers ( )
{
}

void vtkUltrasoundModuleGUI::RemoveMRMLNodeObservers ( )
{
}

void vtkUltrasoundModuleGUI::RemoveLogicObservers ( )
{
}



// Description:
// Process events generated by GUI widgets
void vtkUltrasoundModuleGUI::ProcessGUIEvents ( vtkObject *caller, unsigned long event,
                                               void *callData )
{
}


// Description:
// Process events generated by MRML
void vtkUltrasoundModuleGUI::ProcessMRMLEvents ( vtkObject *caller, unsigned long event,
                                                void *callData)
{
}


void vtkUltrasoundModuleGUI::CreateWidget()
{
}


// Description:
// Methods describe behavior at module enter and exit.
void vtkUltrasoundModuleGUI::Enter ( )
{
}

void vtkUltrasoundModuleGUI::Exit ( )
{
}



void vtkUltrasoundModuleGUI::PrintSelf(ostream& os, vtkIndent indent)
{
}

vtkMRMLVolumeNode* vtkUltrasoundModuleGUI::AddVolumeNode(const char* volumeNodeName)
{

    std::cerr << "AddVolumeNode(): called." << std::endl;

    vtkMRMLVolumeNode *volumeNode = NULL;

    if (volumeNode == NULL)  // if real-time volume node has not been created
    {

        //vtkMRMLVolumeDisplayNode *displayNode = NULL;
        vtkMRMLScalarVolumeDisplayNode *displayNode = NULL;
        vtkMRMLScalarVolumeNode *scalarNode = vtkMRMLScalarVolumeNode::New();
        vtkImageData* image = vtkImageData::New();

        float fov = 300.0;
        image->SetDimensions(256, 256, 1);
        image->SetExtent(0, 255, 0, 255, 0, 0 );
        image->SetSpacing( fov/256, fov/256, 10 );
        image->SetOrigin( -fov/2, -fov/2, -0.0 );
        image->SetScalarTypeToShort();
        image->AllocateScalars();

        short* dest = (short*) image->GetScalarPointer();
        if (dest)
        {
            memset(dest, 0x00, 256*256*sizeof(short));
            image->Update();
        }

        /*
        vtkSlicerSliceLayerLogic *reslice = vtkSlicerSliceLayerLogic::New();
        reslice->SetUseReslice(0);
        */
        scalarNode->SetAndObserveImageData(image);


        /* Based on the code in vtkSlicerVolumeLogic::AddHeaderVolume() */
        //displayNode = vtkMRMLVolumeDisplayNode::New();
        displayNode = vtkMRMLScalarVolumeDisplayNode::New();
        scalarNode->SetLabelMap(0);
        volumeNode = scalarNode;

        if (volumeNode != NULL)
        {
            volumeNode->SetName(volumeNodeName);
            this->GetMRMLScene()->SaveStateForUndo();

            vtkDebugMacro("Setting scene info");
            volumeNode->SetScene(this->GetMRMLScene());
            displayNode->SetScene(this->GetMRMLScene());


            double range[2];
            vtkDebugMacro("Set basic display info");
            volumeNode->GetImageData()->GetScalarRange(range);
            range[0] = 0.0;
            range[1] = 256.0;
            displayNode->SetLowerThreshold(range[0]);
            displayNode->SetUpperThreshold(range[1]);
            displayNode->SetWindow(range[1] - range[0]);
            displayNode->SetLevel(0.5 * (range[1] - range[0]) );

            vtkDebugMacro("Adding node..");
            this->GetMRMLScene()->AddNode(displayNode);

            //displayNode->SetDefaultColorMap();
            vtkSlicerColorLogic *colorLogic = vtkSlicerColorLogic::New();
            displayNode->SetAndObserveColorNodeID(colorLogic->GetDefaultVolumeColorNodeID());
            //colorLogic->Delete();

            volumeNode->SetAndObserveDisplayNodeID(displayNode->GetID());

            vtkDebugMacro("Name vol node "<<volumeNode->GetClassName());
            vtkDebugMacro("Display node "<<displayNode->GetClassName());

            this->GetMRMLScene()->AddNode(volumeNode);
            vtkDebugMacro("Node added to scene");
        }

        //scalarNode->Delete();
        /*
        if (displayNode)
        {
        displayNode->Delete();
        }
        */

    }
    return volumeNode;
}
