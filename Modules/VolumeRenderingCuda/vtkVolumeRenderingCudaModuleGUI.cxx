#include "vtkVolumeRenderingCudaModuleGUI.h"
#include "vtkVolumeRenderingCudaModuleLogic.h"
#include "vtkSlicerApplication.h"
#include "vtkKWWidget.h"
#include "vtkKWPushButton.h"
#include "vtkSlicerNodeSelectorWidget.h"
#include "vtkSlicerModuleCollapsibleFrame.h"
#include "vtkMRMLScene.h"
#include "vtkVolume.h"

#include "vtkImageData.h"

#include "vtkVolumeCudaMapper.h"
#include "vtkKWTypeChooserBox.h"
#include "vtkKWMatrixWidget.h"
#include "vtkKWLabel.h"

#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkTexture.h"

#include "vtkRenderer.h"

#include "vtkCudaMemory.h"

extern "C" {
#include "CUDA_renderAlgo.h"
}
/// TEMPORARY
#include <stdio.h>
#include <stdlib.h>
#include <cutil.h>


vtkVolumeRenderingCudaModuleGUI::vtkVolumeRenderingCudaModuleGUI()
{
    this->LoadButton = NULL;
    this->CreatePiplineTestButton = NULL;
    this->CudaMapper = NULL;
    this->CudaActor = NULL;

    this->ImageData = NULL;
    this->InputTypeChooser = NULL;
    this->InputResolutionMatrix = NULL;
    this->Color = NULL;

    this->CudaInputMemoryCache = NULL;
}


vtkVolumeRenderingCudaModuleGUI::~vtkVolumeRenderingCudaModuleGUI()
{
    if (this->LoadButton != NULL)
    {
        this->LoadButton->SetParent(NULL);
        this->LoadButton->Delete();
        this->LoadButton = NULL; 
    }

    if (this->CreatePiplineTestButton != NULL)
    {
        this->CreatePiplineTestButton->SetParent(NULL);
        this->CreatePiplineTestButton->Delete();
        this->CreatePiplineTestButton = NULL;  
    }

    if (this->CudaMapper != NULL)
    {
        this->CudaMapper->Delete();
    }
    if (this->CudaActor != NULL)
    {
        this->CudaActor->Delete();  
    }

    DeleteWidget(this->InputTypeChooser);
    DeleteWidget(this->InputResolutionMatrix);
    DeleteWidget(this->Color);

    if (this->CudaInputMemoryCache != NULL)
        this->CudaInputMemoryCache->Delete();
}

void vtkVolumeRenderingCudaModuleGUI::DeleteWidget(vtkKWWidget* widget)
{
    if (widget != NULL)
    {
        widget->SetParent(NULL);
        widget->Delete();
    }
}

vtkVolumeRenderingCudaModuleGUI* vtkVolumeRenderingCudaModuleGUI::New()
{
    vtkObject* ret = vtkObjectFactory::CreateInstance("vtkVolumeRenderingCudaModuleGUI");
    if (ret)
        return (vtkVolumeRenderingCudaModuleGUI*)ret;
    // If the Factory was unable to create the object, we do it ourselfes.
    return new vtkVolumeRenderingCudaModuleGUI();
}


void vtkVolumeRenderingCudaModuleGUI::BuildGUI ( )
{
    vtkSlicerApplication *app = (vtkSlicerApplication *)this->GetApplication();
    this->GetUIPanel()->AddPage("VolumeRenderingCuda","VolumeRenderingCuda",NULL);

    // Define your help text and build the help frame here.
    const char *help = "VolumeRenderingCuda. 3D Segmentation This module is currently a prototype and will be under active development throughout 3DSlicer's Beta release.";
    const char *about = "This work is supported by NA-MIC, NAC, BIRN, NCIGT, and the Slicer Community. See http://www.slicer.org for details.";
    vtkKWWidget *page = this->UIPanel->GetPageWidget ( "VolumeRenderingCuda" );
    this->BuildHelpAndAboutFrame ( page, help, about );
    //
    //Load and save
    //
    vtkSlicerModuleCollapsibleFrame *loadSaveDataFrame = vtkSlicerModuleCollapsibleFrame::New ( );
    loadSaveDataFrame->SetParent (page);
    loadSaveDataFrame->Create();
    loadSaveDataFrame->ExpandFrame();
    loadSaveDataFrame->SetLabelText("Load and Save");
    app->Script ( "pack %s -side top -anchor nw -fill x -padx 2 -pady 2 -in %s",
        loadSaveDataFrame->GetWidgetName(), page->GetWidgetName());

    this->LoadButton = vtkKWPushButton::New();
    this->LoadButton->SetParent(loadSaveDataFrame->GetFrame());
    this->LoadButton->Create();
    this->LoadButton->SetText("Load new Model");
    app->Script( "pack %s -side top -anchor nw -fill x -padx 2 -pady 2 -in %s",
        this->LoadButton->GetWidgetName(), loadSaveDataFrame->GetFrame()->GetWidgetName());


    this->CreatePiplineTestButton = vtkKWPushButton::New();
    this->CreatePiplineTestButton->SetParent(loadSaveDataFrame->GetFrame());
    this->CreatePiplineTestButton->Create();
    this->CreatePiplineTestButton->SetText("Test Creating the pipeline");
    app->Script( "pack %s -side top -anchor nw -fill x -padx 2 -pady 2 -in %s",
        this->CreatePiplineTestButton->GetWidgetName(), loadSaveDataFrame->GetFrame()->GetWidgetName());


    this->InputTypeChooser = vtkKWTypeChooserBox::New();
    this->InputTypeChooser->SetParent(loadSaveDataFrame->GetFrame());
    this->InputTypeChooser->Create();
    app->Script( "pack %s -side top -anchor nw -fill x -padx 2 -pady 2 -in %s",
        this->InputTypeChooser->GetWidgetName(), loadSaveDataFrame->GetFrame()->GetWidgetName());

    this->InputResolutionMatrix = vtkKWMatrixWidget::New();
    this->InputResolutionMatrix->SetParent(loadSaveDataFrame->GetFrame());
    this->InputResolutionMatrix->Create();
    this->InputResolutionMatrix->SetRestrictElementValueToInteger();
    this->InputResolutionMatrix->SetNumberOfColumns(4);
    this->InputResolutionMatrix->SetNumberOfRows(1);
    this->InputResolutionMatrix->SetElementValueAsInt(0,0,128);
    this->InputResolutionMatrix->SetElementValueAsInt(0,1,128);
    this->InputResolutionMatrix->SetElementValueAsInt(0,2,1);    
    this->InputResolutionMatrix->SetElementValueAsInt(0,3,4);
    app->Script( "pack %s -side top -anchor nw -fill x -padx 2 -pady 2 -in %s",
        this->InputResolutionMatrix->GetWidgetName(), loadSaveDataFrame->GetFrame()->GetWidgetName());


    /// CameraPosition
    this->CameraPosition = vtkKWMatrixWidget::New();
    this->CameraPosition->SetParent(loadSaveDataFrame->GetFrame());
    this->CameraPosition->Create();
    this->CameraPosition->SetRestrictElementValueToDouble();
    this->CameraPosition->SetNumberOfColumns(3);
    this->CameraPosition->SetNumberOfRows(3);
    this->CameraPosition->SetElementValueAsDouble(0,0, 10);
    this->CameraPosition->SetElementValueAsDouble(0,1, 10);
    this->CameraPosition->SetElementValueAsDouble(0,2, 10);

    this->CameraPosition->SetElementValueAsDouble(1,0, 0);
    this->CameraPosition->SetElementValueAsDouble(1,1, 0);
    this->CameraPosition->SetElementValueAsDouble(1,2, 0);

    this->CameraPosition->SetElementValueAsDouble(2,0, 0);
    this->CameraPosition->SetElementValueAsDouble(2,1, 1);
    this->CameraPosition->SetElementValueAsDouble(2,2, 0);
    app->Script( "pack %s -side top -anchor nw -fill x -padx 2 -pady 2 -in %s",
        this->CameraPosition->GetWidgetName(), loadSaveDataFrame->GetFrame()->GetWidgetName());

    vtkKWLabel* label = vtkKWLabel::New();
    label->SetParent(loadSaveDataFrame->GetFrame());
    label->Create();
    label->SetText("Color:");
    app->Script( "pack %s -side top -anchor nw -fill x -padx 2 -pady 2 -in %s",
        label->GetWidgetName(), loadSaveDataFrame->GetFrame()->GetWidgetName());


    this->Color = vtkKWMatrixWidget::New();
    this->Color->SetParent(loadSaveDataFrame->GetFrame());
    this->Color->Create();
    this->Color->SetRestrictElementValueToInteger();
    this->Color->SetNumberOfColumns(3);
    this->Color->SetNumberOfRows(1);
    this->Color->SetElementValueAsInt(0,0,255);
    this->Color->SetElementValueAsInt(0,1,255);
    this->Color->SetElementValueAsInt(0,2,255);
    app->Script( "pack %s -side top -anchor nw -fill x -padx 2 -pady 2 -in %s",
        this->Color->GetWidgetName(), loadSaveDataFrame->GetFrame()->GetWidgetName());  


    this->UpdateButton = vtkKWPushButton::New();
    this->UpdateButton->SetParent(loadSaveDataFrame->GetFrame());
    this->UpdateButton->Create();
    this->UpdateButton->SetText("Update Renderer");
    app->Script( "pack %s -side top -anchor nw -fill x -padx 2 -pady 2 -in %s",
        this->UpdateButton->GetWidgetName(), loadSaveDataFrame->GetFrame()->GetWidgetName());

    this->Built=true;
}

void vtkVolumeRenderingCudaModuleGUI::TearDownGUI ( )
{
    this->Exit();
    if ( this->Built )
    {
        this->RemoveGUIObservers();
    }
}

void vtkVolumeRenderingCudaModuleGUI::CreateModuleEventBindings ( )
{
    vtkDebugMacro("VolumeRenderingCudaModule: CreateModuleEventBindings: No ModuleEventBindings yet");
}
void vtkVolumeRenderingCudaModuleGUI::ReleaseModuleEventBindings ( )
{
    vtkDebugMacro("VolumeRenderingCudaModule: ReleaseModuleEventBindings: No ModuleEventBindings to remove yet");
}

void vtkVolumeRenderingCudaModuleGUI::AddGUIObservers ( )
{
    this->LoadButton->AddObserver(vtkKWPushButton::InvokedEvent, (vtkCommand*)this->GUICallbackCommand);
    this->CreatePiplineTestButton->AddObserver(vtkKWPushButton::InvokedEvent, (vtkCommand*)this->GUICallbackCommand);

    this->InputTypeChooser->GetMenu()->AddObserver(vtkKWMenu::MenuItemInvokedEvent, (vtkCommand*)this->GUICallbackCommand);
    this->InputResolutionMatrix->AddObserver(vtkKWMatrixWidget::ElementChangedEvent, (vtkCommand*)this->GUICallbackCommand);
    this->CameraPosition->AddObserver(vtkKWMatrixWidget::ElementChangedEvent, (vtkCommand*)this->GUICallbackCommand);
    this->Color->AddObserver(vtkKWMatrixWidget::ElementChangedEvent, (vtkCommand*)this->GUICallbackCommand);
    this->UpdateButton->AddObserver(vtkKWPushButton::InvokedEvent, (vtkCommand*)this->GUICallbackCommand);
}

void vtkVolumeRenderingCudaModuleGUI::RemoveGUIObservers ( )
{
    this->LoadButton->RemoveObservers(vtkKWPushButton::InvokedEvent, (vtkCommand*)this->GUICallbackCommand);
    this->CreatePiplineTestButton->RemoveObservers(vtkKWPushButton::InvokedEvent, (vtkCommand*)this->GUICallbackCommand);

    this->InputTypeChooser->GetMenu()->RemoveObservers(vtkKWMenu::MenuItemInvokedEvent, (vtkCommand*)this->GUICallbackCommand);
    this->InputResolutionMatrix->RemoveObservers(vtkKWMatrixWidget::ElementChangedEvent, (vtkCommand*)this->GUICallbackCommand);
    this->CameraPosition->RemoveObservers(vtkKWMatrixWidget::ElementChangedEvent, (vtkCommand*)this->GUICallbackCommand);
    this->Color->RemoveObservers(vtkKWMatrixWidget::ElementChangedEvent, (vtkCommand*)this->GUICallbackCommand);
    this->UpdateButton->RemoveObservers(vtkKWPushButton::InvokedEvent, (vtkCommand*)this->GUICallbackCommand);
}
void vtkVolumeRenderingCudaModuleGUI::RemoveMRMLNodeObservers ( )
{
}
void vtkVolumeRenderingCudaModuleGUI::RemoveLogicObservers ( )
{
}

void vtkVolumeRenderingCudaModuleGUI::ProcessGUIEvents ( vtkObject *caller, unsigned long event,
                                                        void *callData )
{
    vtkDebugMacro("vtkVolumeRenderingModuleGUI::ProcessGUIEvents: event = " << event);

    if (caller == this->LoadButton)
    {
        if (this->CudaMapper == NULL)
            this->CudaMapper = vtkVolumeCudaMapper::New();
        if (this->CudaActor == NULL)
        {
            this->CudaActor = vtkVolume::New();
            this->CudaActor->SetMapper(this->CudaMapper);
        }

        this->TestCudaViewer();
        this->CudaMapper->Render(NULL, NULL);
    }

    if (caller == this->CreatePiplineTestButton)
    {
        this->CreatePipelineTest();
    }


    /// INPUT TYPE OR SIZE CHANGED CHANGED
    if (
        /*
        caller == this->InputTypeChooser->GetMenu() ||
        caller == this->InputResolutionMatrix ||
        caller == this->Color ||
        caller == this->CameraPosition
        */
        caller == this->UpdateButton
        )
    {
        cerr << "Type" << this->InputTypeChooser->GetSelectedName() << " " << this->InputTypeChooser->GetSelectedType() << 
            " X:" << this->InputResolutionMatrix->GetElementValueAsInt(0,0) << 
            " Y:" << this->InputResolutionMatrix->GetElementValueAsInt(0,1) <<
            " Z:" << this->InputResolutionMatrix->GetElementValueAsInt(0,2) <<
            " A:" << this->InputResolutionMatrix->GetElementValueAsInt(0,3) << endl;

        TestCudaViewer();

        try {
            cerr << "ALLOCATE" << endl;  
            //this->ImageData->SetScalarType(this->InputTypeChooser->GetSelectedType());
            this->ImageData->SetDimensions(this->InputResolutionMatrix->GetElementValueAsInt(0,0), this->InputResolutionMatrix->GetElementValueAsInt(0,0), 1);
            this->ImageData->SetExtent(0, this->InputResolutionMatrix->GetElementValueAsInt(0,0) - 1, 
                0, this->InputResolutionMatrix->GetElementValueAsInt(0,1) - 1, 
                0, this->InputResolutionMatrix->GetElementValueAsInt(0,2) - 1);
            this->ImageData->SetNumberOfScalarComponents(this->InputResolutionMatrix->GetElementValueAsInt(0,3));
            this->ImageData->SetScalarTypeToUnsignedChar();
            this->ImageData->AllocateScalars();

            this->RenderWithCUDA("C:\\Documents and Settings\\bensch\\Desktop\\svn\\orxonox\\subprojects\\volrenSample\\heart256.raw", 256, 256, 256);
            
            cerr << "FINISHED" << endl;
        }
        catch (...)
        {
            cerr << "ERROR READING" << endl;  
        }
    }
}



void vtkVolumeRenderingCudaModuleGUI::RenderWithCUDA(const char* inputFile, int inX, int inY, int inZ)
{  
    int outX = this->InputResolutionMatrix->GetElementValueAsInt(0,0);
    int outY = this->InputResolutionMatrix->GetElementValueAsInt(0,1);


    printf("Input '%s' rendered to 'output.raw' with resolution of %dx%d\n", inputFile, outX, outY);

    unsigned char* inputBuffer=(unsigned char*)malloc(inX*inY*inZ*sizeof(unsigned char));
    unsigned char* outputBuffer;


    /// LOAD THE DATA ONCE!
    if (this->CudaInputMemoryCache == NULL)
    {
        FILE *fp;
        fp=fopen(inputFile,"r");
        fread(inputBuffer, sizeof(unsigned char), inX*inY*inZ, fp);
        fclose(fp);

        // Initialization. Prepare and allocate GPU memory to accomodate 3D data and Result image.
        cerr << "CUDA Initialization.\n";
        CUDArenderAlgo_init(inX, inY, inZ, outX,outY);

        this->CudaInputMemoryCache = vtkCudaMemory::New();
        //Allocate and Copy
        this->CudaInputMemoryCache->CopyFromMemory(inputBuffer, sizeof(unsigned char) * inX * inY * inZ);
    }

    // Setting transformation matrix. This matrix will be used to do rotation and translation on ray tracing.

    float color[6]={this->Color->GetElementValueAsInt(0,0),
        this->Color->GetElementValueAsInt(0,1), 
        this->Color->GetElementValueAsInt(0,2),1,1,1};
    float minmax[6]={0,255,0,255,0,255};
    float lightVec[3]={0, 0, 1};

    vtkCamera* cam =
        this->GetApplicationGUI()->GetViewerWidget()->GetMainViewer()->GetRenderer()->GetActiveCamera();
    
    cerr << cam->GetViewTransformMatrix();

    double ax,ay,az;
    double bx,by,bz;
    double cx,cy,cz;

    ax = cam->GetFocalPoint()[0] - cam->GetPosition()[0];
    ay = cam->GetFocalPoint()[1] - cam->GetPosition()[1];
    az = cam->GetFocalPoint()[2] - cam->GetPosition()[2];
    cam->GetViewUp(bx, by, bz);
    cx = ay*bz-az*by;
    cy = az*bx-ax*bz;
    cz = ax*by-ay*bx;
    
    bx = cy*az-cz*ay;
    by = cz*ax-cx*az;
    bz = cx*ay-cy*ax;

    double len = sqrt(ax*ax + ay*ay + az*az);
    ax /= len; ay /= len; az /= len;

    len = sqrt(bx*bx + by*by + bz*bz);
    bx /= len; by /= len; bz /= len;

    len = sqrt(cx*cx + cy*cy + cz*cz);
    cx /= len; cy /= len; cz /= len;
    //cam = vtkCamera::New();
    //cam->SetPosition(  this->CameraPosition->GetElementValueAsDouble(0,0),
    //    this->CameraPosition->GetElementValueAsDouble(0,1),
    //    this->CameraPosition->GetElementValueAsDouble(0,2));
    //cam->SetViewUp(    this->CameraPosition->GetElementValueAsDouble(1,0),
    //    this->CameraPosition->GetElementValueAsDouble(1,1),
    //    this->CameraPosition->GetElementValueAsDouble(1,2));
    //cam->SetFocalPoint(this->CameraPosition->GetElementValueAsDouble(2,0),
    //    this->CameraPosition->GetElementValueAsDouble(2,1),
    //    this->CameraPosition->GetElementValueAsDouble(2,2));

    vtkMatrix4x4*  viewMat = vtkMatrix4x4::New();
    //cam->GetPerspectiveTransformMatrix(1.0,.1,1000);
    viewMat->DeepCopy(cam->GetViewTransformMatrix());//1.0, .1, 1000));
    //  viewMat->Invert();
    //viewMat->Invert();
    cerr << *viewMat;


    float rotationMatrix[4][4]=
   /* {{1,0,0,0},
    {0,1,0,0},
    {0,0,1,0},
    {0,0,0,1}};*/
    {{ax,bx,cx,2},
     {ay,by,cy,2},
     {az,bz,cz,2},
     {0,0,0,1}};
    /*
    for (unsigned int i = 0; i < 4; i++)
    {
    for (unsigned int j = 0; j < 4; j++)
    rotationMatrix[i][j] = viewMat->GetElement(i,j);
    }
    */


    cerr << "Volume rendering.\n";
    // Do rendering. 

    CUDArenderAlgo_doRender(this->CudaInputMemoryCache->GetMemPointerAs<unsigned char>(),
        (float*)rotationMatrix, color, minmax, lightVec, 
        inX, inY, inZ,    //3D data size
        outX, outY,     //result image size
        0,0,0,          //translation of data in x,y,z direction
        1, 1, 1,        //voxel dimension
        90, 255,        //min and max threshold
        -100);          //slicing distance from center of 3D data
    // Get the resulted image.

    cerr << "Copy result from GPU to CPU.\n";
    CUDArenderAlgo_getResult((unsigned char**)&outputBuffer, outX,outY);
    memcpy(this->ImageData->GetScalarPointer(), outputBuffer,
        sizeof(unsigned char) *
        this->InputResolutionMatrix->GetElementValueAsInt(0,0) *
        this->InputResolutionMatrix->GetElementValueAsInt(0,1) *
        this->InputResolutionMatrix->GetElementValueAsInt(0,2) * 
        this->InputResolutionMatrix->GetElementValueAsInt(0,3));

    this->RenderToScreen(this->ImageData);

    // Free allocated GPU memory.
    CUDArenderAlgo_delete();
    this->CudaInputMemoryCache->Delete();
    this->CudaInputMemoryCache = NULL;
    free(inputBuffer);
}


void vtkVolumeRenderingCudaModuleGUI::RenderToScreen(vtkImageData* imageData)
{
vtkImageExtractComponents *components=vtkImageExtractComponents::New();
        components->SetInput(imageData);
        components->SetComponents(0,1,2);

    vtkRenderer* renPlane = this->GetApplicationGUI()->GetViewerWidget()->GetMainViewer()->GetRenderer();

    vtkRenderWindow *renWin=this->GetApplicationGUI()->GetViewerWidget()->GetMainViewer()->GetRenderWindow();
        //Get current size of window
        int *size=renWin->GetSize();
    
    //renPlane->SetBackground(this->renViewport->GetBackground());
    //renPlane->SetActiveCamera(this->renViewport->GetActiveCamera());

    renPlane->SetDisplayPoint(0,0,0.5);
    renPlane->DisplayToWorld();
    double coordinatesA[4];
    renPlane->GetWorldPoint(coordinatesA);

    renPlane->SetDisplayPoint(size[0],0,0.5);
    renPlane->DisplayToWorld();
    double coordinatesB[4];
    renPlane->GetWorldPoint(coordinatesB);

    renPlane->SetDisplayPoint(size[0],size[1],0.5);
    renPlane->DisplayToWorld();
    double coordinatesC[4];
    renPlane->GetWorldPoint(coordinatesC);

    renPlane->SetDisplayPoint(0,size[1],0.5);
    renPlane->DisplayToWorld();
    double coordinatesD[4];
    renPlane->GetWorldPoint(coordinatesD);

    //Create the Polydata
    vtkPoints *points=vtkPoints::New();
    points->InsertPoint(0,coordinatesA);
    points->InsertPoint(1,coordinatesB);
    points->InsertPoint(2,coordinatesC);
    points->InsertPoint(3,coordinatesD);

    vtkCellArray *polygon=vtkCellArray::New();
    polygon->InsertNextCell(4);
    polygon->InsertCellPoint(0);
    polygon->InsertCellPoint(1);
    polygon->InsertCellPoint(2);
    polygon->InsertCellPoint(3);
    //Take care about Texture coordinates
    vtkFloatArray *textCoords=vtkFloatArray::New();
    textCoords->SetNumberOfComponents(2);
    textCoords->Allocate(8);
    float tc[2];
    tc[0]=0;
    tc[1]=0;
    textCoords->InsertNextTuple(tc);
    tc[0]=1;
    tc[1]=0;
    textCoords->InsertNextTuple(tc);
    tc[0]=1;
    tc[1]=1;
    textCoords->InsertNextTuple(tc);
    tc[0]=0;
    tc[1]=1;
    textCoords->InsertNextTuple(tc);

    vtkPolyData *polydata=vtkPolyData::New();
    polydata->SetPoints(points);
    polydata->SetPolys(polygon);
    polydata->GetPointData()->SetTCoords(textCoords);

    vtkPolyDataMapper *polyMapper=vtkPolyDataMapper::New();
    polyMapper->SetInput(polydata);

    vtkActor *actor=vtkActor::New(); 
    actor->SetMapper(polyMapper);

    //Take care about the texture
    vtkTexture *atext=vtkTexture::New();
    atext->SetInput(components->GetOutput());
    atext->SetInterpolate(1);
    actor->SetTexture(atext);

    //Remove all old Actors
    renPlane->RemoveAllViewProps();

    renPlane->AddActor(actor);
    //Remove the old Renderer

    renWin->SwapBuffersOn();


    //Delete everything we have done
    components->Delete();
    points->Delete();
    polygon->Delete();
    textCoords->Delete();
    polydata->Delete();
    polyMapper->Delete();
    actor->Delete();
    atext->Delete();
    this->GetApplicationGUI()->GetViewerWidget()->GetMainViewer()->GetRenderWindow()->Render();

}


#include "vtkImageReader.h"

void vtkVolumeRenderingCudaModuleGUI::TestCudaViewer()
{
    if (this->ImageData == NULL) { 
    printf("START CREATING WINDOW\n");

        /// ImageData (from File) in Memory.
        this->ImageData = vtkImageData::New();
        //this->ImageData->SetDimensions(256, 256, 1);
        this->ImageData->SetScalarTypeToUnsignedChar();
        this->ImageData->SetNumberOfScalarComponents(4);

        printf("FINISHED CREATING WINDOW\n");
    }
}

void vtkVolumeRenderingCudaModuleGUI::CreatePipelineTest()
{
    vtkImageReader* reader = vtkImageReader::New();
    reader->SetFileName("/projects/igtdev/bensch/svn/volrenSample/heart256.raw");




    reader->Delete();
}

void vtkVolumeRenderingCudaModuleGUI::ProcessMRMLEvents ( vtkObject *caller, unsigned long event,
                                                         void *callData)
{
}


void vtkVolumeRenderingCudaModuleGUI::SetViewerWidget(vtkSlicerViewerWidget *viewerWidget)
{
}
void vtkVolumeRenderingCudaModuleGUI::SetInteractorStyle(vtkSlicerViewerInteractorStyle *interactorStyle)
{
}


void vtkVolumeRenderingCudaModuleGUI::Enter ( )
{
    vtkDebugMacro("Enter Volume Rendering Cuda Module");

    if ( this->Built == false )
    {
        this->BuildGUI();
        this->AddGUIObservers();
    }
    this->CreateModuleEventBindings();
    //this->UpdateGUI();
}
void vtkVolumeRenderingCudaModuleGUI::Exit ( )
{
    vtkDebugMacro("Exit: removeObservers for VolumeRenderingModule");
    this->ReleaseModuleEventBindings();
}


void vtkVolumeRenderingCudaModuleGUI::PrintSelf(ostream& os, vtkIndent indent)
{
    this->SuperClass::PrintSelf(os, indent);

    os<<indent<<"vtkVolumeRenderingCudaModuleGUI"<<endl;
    os<<indent<<"vtkVolumeRenderingCudaModuleLogic"<<endl;
    if(this->GetLogic())
    {
        this->GetLogic()->PrintSelf(os,indent.GetNextIndent());
    }
}
