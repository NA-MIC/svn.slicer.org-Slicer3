#include "vtkVolumeCudaMapper.h"
#include "vtkVolumeRenderingCudaFactory.h"

#include "vtkActor.h"
#include "vtkPolyDataMapper.h"

#include "vtkPolyData.h"
#include "vtkTexture.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkImageExtractComponents.h"

#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkCamera.h"

#include "vtkImageData.h"
#include "vtkCudaMemory.h"
#include "vtkCudaHostMemory.h"
#include "vtkCudaMemoryArray.h"

#include "vtkObjectFactory.h"

#include "vtkCudaMemory.h"
#include <vector_types.h>

extern "C" {
#include "CUDA_renderAlgo.h"
}



vtkCxxRevisionMacro(vtkVolumeCudaMapper, "$Revision: 1.6 $");
vtkStandardNewMacro(vtkVolumeCudaMapper);

vtkVolumeCudaMapper::vtkVolumeCudaMapper()
{
    this->LocalOutputImage = vtkImageData::New();

    this->CudaInputBuffer = vtkCudaMemory::New();
    this->CudaOutputBuffer = vtkCudaMemory::New();
    
    this->Test = vtkCudaMemoryArray::New();

    this->OutputDataSize[0] = this->OutputDataSize[1] = 0;
    this->SetColor(255, 255, 255);
    this->UpdateOutputResolution(128, 128, 4);
}  

vtkVolumeCudaMapper::~vtkVolumeCudaMapper()
{
    this->CudaInputBuffer->Delete();
    this->CudaOutputBuffer->Delete();
    this->LocalOutputImage->Delete();
    this->Test->Delete();
}

void vtkVolumeCudaMapper::SetInput(vtkImageData * input)
{
    this->Superclass::SetInput(input);

    if (input != NULL)
    {
        this->CudaInputBuffer->CopyFrom(input);
    }
    else
    {
        this->CudaInputBuffer->Free();
    }
}


void vtkVolumeCudaMapper::UpdateOutputResolution(unsigned int width, unsigned int height, unsigned int colors)
{
    if (this->OutputDataSize[0] == width &&
        this->OutputDataSize[1] == height )
        return;
    // Set the data Size
    this->OutputDataSize[0] = width;
    this->OutputDataSize[1] = height;

    // Re-allocate the memory
    this->CudaOutputBuffer->Allocate<uchar4>(width * height);
    this->Test->SetFormat<uchar4>();
    this->Test->Allocate(width, height);

    // Allocate the Image Data
    this->LocalOutputImage->SetScalarTypeToUnsignedChar();
    this->LocalOutputImage->SetNumberOfScalarComponents(4);
    this->LocalOutputImage->SetDimensions(width, height, 1);
    this->LocalOutputImage->SetExtent(0, width - 1, 
        0, height - 1, 
        0, 1 - 1);
    this->LocalOutputImage->SetNumberOfScalarComponents(colors);
    this->LocalOutputImage->SetScalarTypeToUnsignedChar();
    this->LocalOutputImage->AllocateScalars();
}

void vtkVolumeCudaMapper::Update()
{
    cerr << "TEST\n";
    cout << "TEST2\n";
}


void vtkVolumeCudaMapper::Render(vtkRenderer *renderer, vtkVolume *volume)
{
    cerr << "TEST\n";
    cout << "TEST2\n";
    this->UpdateRenderPlane(renderer, volume);
}


#include "vtkTimerLog.h"
#include "texture_types.h"

void vtkVolumeCudaMapper::UpdateRenderPlane(vtkRenderer* renderer, vtkVolume* volume)
{
    float color[6]={this->Color[0],this->Color[1],this->Color[2], 1,1,1};
    float minmax[6]={0,255,0,255,0,255};
    float lightVec[3]={0, 0, 1};

    vtkRenderWindow *renWin= renderer->GetRenderWindow();
    //Get current size of window
    int *size=renWin->GetSize();

    int width = size[0], height = size[1];
    this->UpdateOutputResolution(width, height, 4);

    vtkCamera* cam =
        renderer->GetActiveCamera();

    // Build the Rotation Matrix
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

    double distance = cam->GetDistance();
    ax /= distance; ay /= distance; az /= distance;

    double len = sqrt(bx*bx + by*by + bz*bz);
    bx /= len; by /= len; bz /= len;

    len = sqrt(cx*cx + cy*cy + cz*cz);
    cx /= len; cy /= len; cz /= len;

    float rotationMatrix[4][4]=
    {{ax,bx,cx,0},
    {ay,by,cy,0},
    {az,bz,cz,0},
    {0,0,0,1}};
    /*  {{1,0,0,0},
    {0,1,0,0},
    {0,0,1,0},
    {0,0,0,1}};*/

    vtkErrorMacro( << "Volume rendering.\n");
    // Do rendering.

    int* dims = this->GetInput()->GetDimensions();

    vtkTimerLog* log = vtkTimerLog::New();
    log->StartTimer();
    CUDArenderAlgo_doRender(this->CudaOutputBuffer->GetMemPointerAs<uchar4>(),
        this->CudaInputBuffer->GetMemPointerAs<unsigned char>(),
        (float*)rotationMatrix, color, minmax, lightVec, 
        dims[0], dims[1], dims[2],    //3D data size
        this->OutputDataSize[0], this->OutputDataSize[1],     //result image size
        0,0,0,          //translation of data in x,y,z direction
        1, 1, 1,        //voxel dimension
        90, 255,        //min and max threshold
        -100);          //slicing distance from center of 3D data
    // Get the resulted image.
    log->StopTimer();

    vtkErrorMacro(<< "Elapsed Time to Render:: " << log->GetElapsedTime());
    log->StartTimer();
    this->CudaOutputBuffer->CopyTo(this->LocalOutputImage);
    log->StopTimer();
    vtkErrorMacro(<< "Elapsed Time to Copy Memory:: " << log->GetElapsedTime());

    log->StartTimer();

    /// EXPERIMENTAL BEGIN

    //textureReference tex;

    //tex.channelDesc = cudaCreateChannelDesc<uchar4>();
    //tex.addressMode[0] = cudaAddressModeWrap;
    //tex.addressMode[1] = cudaAddressModeWrap;
    //tex.filterMode = cudaFilterModeLinear;
    //tex.normalized = true;

    //log->StartTimer();
    //this->CudaOutputBuffer->CopyTo(this->Test);
    //cudaBindTextureToArray(&tex, this->Test->GetArray(), &tex.channelDesc);
    //log->StopTimer();
    //vtkErrorMacro(<< "Elapsed Time to Copy Memory on Device:: " << log->GetElapsedTime());
    //



//      
//// TEXTURE CODE
//        printf("Creating GL texture...\n");
//        glEnable(GL_TEXTURE_2D);
//        glGenTextures(1, &gl_Tex);
//        glBindTexture(GL_TEXTURE_2D, gl_Tex);
//        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
//        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
//        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
//        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
//        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, this->LocalOutputBuffer);
//    printf("Texture created.\n");
//
//    printf("Creating PBO...\n");
//        glGenBuffers(1, &gl_PBO);
//        glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, gl_PBO);
//        glBufferData(GL_PIXEL_UNPACK_BUFFER_ARB, width * height * 4, h_Src, GL_STREAM_COPY);
//        While a PBO is registered to CUDA, it can't be used 
//        as the destination for OpenGL drawing calls.
//        But in our particular case OpenGL is only used 
//        to display the content of the PBO, specified by CUDA kernels,
//        so we need to register/unregister it only once.
//        CUDA_SAFE_CALL( cudaGLRegisterBufferObject(gl_PBO) );
//    printf("PBO created.\n");

    /// EXPERIMENTAL END

    vtkImageExtractComponents *components = vtkImageExtractComponents::New();
    components->SetInput(this->LocalOutputImage);
    components->SetComponents(0,1,2);

    

    //renderer->SetBackground(this->renViewport->GetBackground());
    //renderer->SetActiveCamera(this->renViewport->GetActiveCamera());

    renderer->SetDisplayPoint(0,0,0.5);
    renderer->DisplayToWorld();
    double coordinatesA[4];
    renderer->GetWorldPoint(coordinatesA);

    renderer->SetDisplayPoint(size[0],0,0.5);
    renderer->DisplayToWorld();
    double coordinatesB[4];
    renderer->GetWorldPoint(coordinatesB);

    renderer->SetDisplayPoint(size[0],size[1],0.5);
    renderer->DisplayToWorld();
    double coordinatesC[4];
    renderer->GetWorldPoint(coordinatesC);

    renderer->SetDisplayPoint(0,size[1],0.5);
    renderer->DisplayToWorld();
    double coordinatesD[4];
    renderer->GetWorldPoint(coordinatesD);

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
    vtkFloatArray *textCoords = vtkFloatArray::New();
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
    renderer->RemoveAllViewProps();

    renderer->AddActor(actor);
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

       log->StopTimer();
    vtkErrorMacro(<< "FINISH:: " << log->GetElapsedTime());

}

void vtkVolumeCudaMapper::PrintSelf(ostream& os, vtkIndent indent)
{
    vtkVolumeMapper::PrintSelf(os, indent);
}
