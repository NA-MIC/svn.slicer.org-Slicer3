// Type
#include "vtkVolumeCudaMapper.h"
#include "vtkVolumeRenderingCudaFactory.h"
#include "vtkObjectFactory.h"
// Extended Type
#include "vtkVolume.h"
#include "vtkPolyDataMapper.h"
#include "vtkVolumeProperty.h"
#include "vtkColorTransferFunction.h"
#include "vtkPiecewiseFunction.h"

//Data Types
#include "vtkPolyData.h"
#include "vtkTexture.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkImageExtractComponents.h"

// Rendering
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkCamera.h"


// CUDA
#include "vtkCudaSupport.h"
#include "vtkImageData.h"
#include "vtkCudaMemory.h"
#include "vtkCudaHostMemory.h"
#include "vtkCudaMemoryArray.h"

#include <vector_types.h>

// openGL
#include "vtkOpenGLExtensionManager.h"
#include "vtkgl.h"
#include "cuda_gl_interop.h"

extern "C" {
#include "CUDA_renderAlgo.h"
}
#include <cutil.h>


vtkCxxRevisionMacro(vtkVolumeCudaMapper, "$Revision: 1.6 $");

vtkVolumeCudaMapper* vtkVolumeCudaMapper::New()
{
    vtkCudaSupport* support = vtkCudaSupport::New();
    bool cudaIsSupported = support->IsSupported();
    support->Delete();
    if (cudaIsSupported)
        return new vtkVolumeCudaMapper;
    else
        return new vtkVolumeCudaMapper; // HACK

}

vtkVolumeCudaMapper::vtkVolumeCudaMapper()
{
    this->LocalOutputImage = vtkImageData::New();

    this->CudaInputBuffer = vtkCudaMemory::New();
    this->CudaOutputBuffer = vtkCudaMemory::New();

    this->Texture = 0;
    this->BufferObject = 0;

    this->SetThreshold(90,255);

    this->OutputDataSize[0] = this->OutputDataSize[1] = 0;

    // check for the RenderMode
    vtkOpenGLExtensionManager *extensions = vtkOpenGLExtensionManager::New();
    extensions->SetRenderWindow(NULL);
    if (extensions->ExtensionSupported("GL_ARB_vertex_buffer_object"))
    {
        extensions->LoadExtension("GL_ARB_vertex_buffer_object");
        this->GLBufferObjectsAvailiable = true;
        this->CurrentRenderMode = RenderToTexture;
    }
    else
    {
        this->GLBufferObjectsAvailiable = false;
        this->CurrentRenderMode = RenderToMemory;
    }
    extensions->Delete();

    this->CudaColorTransferFunction = vtkCudaMemory::New();
    this->CudaColorTransferFunction->Allocate<float3>(256);
    this->LocalColorTransferFunction = vtkCudaHostMemory::New();
    this->LocalColorTransferFunction->Allocate<float3>(256);

    this->CudaAlphaTransferFunction = vtkCudaMemory::New();
    this->CudaAlphaTransferFunction->Allocate<float>(256);
    this->LocalAlphaTransferFunction = vtkCudaHostMemory::New();
    this->LocalAlphaTransferFunction->Allocate<float>(256);

}  

vtkVolumeCudaMapper::~vtkVolumeCudaMapper()
{
    this->CudaInputBuffer->Delete();
    this->CudaOutputBuffer->Delete();
    this->LocalOutputImage->Delete();

    this->CudaColorTransferFunction->Delete();
    this->LocalColorTransferFunction->Delete();
    this->CudaAlphaTransferFunction->Delete();
    this->LocalAlphaTransferFunction->Delete();

    if (this->Texture == 0 || !glIsTexture(this->Texture))
        glGenTextures(1, &this->Texture);

    if (this->GLBufferObjectsAvailiable == true)
        if (this->BufferObject != 0 && vtkgl::IsBufferARB(this->BufferObject))
            vtkgl::DeleteBuffersARB(1, &this->BufferObject);
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

void vtkVolumeCudaMapper::SetRenderMode(vtkVolumeCudaMapper::RenderMode mode)
{
    if (mode == RenderToTexture && this->GLBufferObjectsAvailiable)
    {
        this->CurrentRenderMode = mode;
    }
    else
    {
        this->CurrentRenderMode = RenderToMemory;
    }
    this->UpdateOutputResolution(this->OutputDataSize[0], this->OutputDataSize[1], true);
}

void vtkVolumeCudaMapper::UpdateOutputResolution(unsigned int width, unsigned int height, bool TypeChanged)
{
    if (this->OutputDataSize[0] == width &&
        this->OutputDataSize[1] == height && !TypeChanged)
        return;

    // Set the data Size
    this->OutputDataSize[0] = width;
    this->OutputDataSize[1] = height;

    // Re-allocate the memory
    this->CudaOutputBuffer->Allocate<uchar4>(width * height);

    {
        // Allocate the Image Data
        this->LocalOutputImage->SetScalarTypeToUnsignedChar();
        this->LocalOutputImage->SetNumberOfScalarComponents(4);
        this->LocalOutputImage->SetDimensions(width, height, 1);
        this->LocalOutputImage->SetExtent(0, width - 1, 
            0, height - 1, 
            0, 1 - 1);
        this->LocalOutputImage->SetNumberOfScalarComponents(4);
        this->LocalOutputImage->AllocateScalars();
    }
    // TEXTURE CODE
    glEnable(GL_TEXTURE_2D);
    if (this->Texture != 0 && glIsTexture(this->Texture))
        glDeleteTextures(1, &this->Texture);
    glGenTextures(1, &this->Texture);
    glBindTexture(GL_TEXTURE_2D, this->Texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, this->LocalOutputImage->GetScalarPointer());
    glBindTexture(GL_TEXTURE_2D, 0);

    if (this->CurrentRenderMode == RenderToTexture)
    {
        // OpenGL Buffer Code
        if (this->BufferObject != 0 && vtkgl::IsBufferARB(this->BufferObject))
            vtkgl::DeleteBuffersARB(1, &this->BufferObject);
        vtkgl::GenBuffersARB(1, &this->BufferObject);
        vtkgl::BindBufferARB(vtkgl::PIXEL_UNPACK_BUFFER_ARB, this->BufferObject);
        vtkgl::BufferDataARB(vtkgl::PIXEL_UNPACK_BUFFER_ARB, width * height * 4, this->LocalOutputImage->GetScalarPointer(), vtkgl::STREAM_COPY);
        vtkgl::BindBufferARB(vtkgl::PIXEL_UNPACK_BUFFER_ARB, 0);
    }
}

void vtkVolumeCudaMapper::UpdateVolumeProperties(vtkVolumeProperty *property)
{
    double range[2];
    property->GetRGBTransferFunction()->GetRange(range);
    property->GetRGBTransferFunction()->GetTable(range[0], range[1], 256, this->LocalColorTransferFunction->GetMemPointerAs<float>());
    LocalColorTransferFunction->CopyTo(CudaColorTransferFunction);

    property->GetScalarOpacity()->GetTable(range[0], range[1], 256, this->LocalAlphaTransferFunction->GetMemPointerAs<float>());
    LocalAlphaTransferFunction->CopyTo(CudaAlphaTransferFunction);

}

#include "vtkTimerLog.h"
#include "texture_types.h"

#include "vtkPNGWriter.h"
void vtkVolumeCudaMapper::Render(vtkRenderer *renderer, vtkVolume *volume)
{
    // Renice!!
    int* minMax = this->GetInput()->GetExtent();
    float minmax[6]; //={minMax[0],minMax[1],minMax[2],minMax[3],minMax[4],minMax[5]};
    for (unsigned int i = 0; i < 6; i++)
        minmax[i] = minMax[i];
    float lightVec[3]={0, 0, 1};


//    this->CudaInputBuffer->CopyFrom(this->GetInput());


    vtkRenderWindow *renWin= renderer->GetRenderWindow();
    //Get current size of window
    int *size=renWin->GetSize();

    int width = size[0], height = size[1];
    this->UpdateOutputResolution(width, height);

    this->UpdateVolumeProperties(volume->GetProperty());

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

    int* dims = this->GetInput()->GetDimensions();

    // Do rendering.
    uchar4* RenderDestination = NULL;

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, this->Texture);

    if (this->CurrentRenderMode == RenderToTexture)
    {
        CUDA_SAFE_CALL( cudaGLRegisterBufferObject(this->BufferObject) );

        vtkgl::BindBufferARB(vtkgl::PIXEL_UNPACK_BUFFER_ARB, this->BufferObject);
        CUDA_SAFE_CALL(cudaGLMapBufferObject((void**)&RenderDestination, this->BufferObject));
    }
    else
    {
        RenderDestination = this->CudaOutputBuffer->GetMemPointerAs<uchar4>();
    }
    vtkTimerLog* log = vtkTimerLog::New();
    log->StartTimer();

    CUDArenderAlgo_doRender(RenderDestination,
        this->CudaInputBuffer->GetMemPointer(),
        this->GetInput()->GetScalarType(),
        (float*)rotationMatrix,
        this->CudaColorTransferFunction->GetMemPointerAs<float>(),
        this->CudaAlphaTransferFunction->GetMemPointerAs<float>(),
        minmax, lightVec, 
        dims[0], dims[1], dims[2],                            //3D data size
        this->OutputDataSize[0], this->OutputDataSize[1],     //result image size
        0,0,0,                                                //translation of data in x,y,z direction
        1, 1, 1,                                              //voxel dimension
        this->Threshold[0], this->Threshold[1],               //min and max threshold
        -500                                                  //slicing distance from center of 3D data
        );         

    // Get the resulted image.
    log->StopTimer();
    //vtkErrorMacro(<< "Elapsed Time to Render:: " << log->GetElapsedTime());

    log->StartTimer();
    if (this->CurrentRenderMode == RenderToTexture)
    {
        CUDA_SAFE_CALL(cudaGLUnmapBufferObject(this->BufferObject));
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, (0));
        CUDA_SAFE_CALL( cudaGLUnregisterBufferObject(this->BufferObject) );
        vtkgl::BindBufferARB(vtkgl::PIXEL_UNPACK_BUFFER_ARB, 0);
    }
    else // (this->CurrentRenderMode == RenderToMemory)
    {
        this->CudaOutputBuffer->CopyTo(this->LocalOutputImage);
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, this->LocalOutputImage->GetScalarPointer());
    }
    log->StopTimer();
    //vtkErrorMacro(<< "Elapsed Time to Copy Memory:: " << log->GetElapsedTime());

    log->StartTimer();

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


    glPushAttrib(GL_BLEND);
    glEnable(GL_BLEND);
    glPushAttrib(GL_LIGHTING);
    glDisable(GL_LIGHTING);

    glBegin(GL_QUADS);
    glTexCoord2i(0,1);
    glVertex4dv(coordinatesA);
    glTexCoord2i(1,1);
    glVertex4dv(coordinatesB);
    glTexCoord2i(1,0);
    glVertex4dv(coordinatesC);
    glTexCoord2i(0,0);
    glVertex4dv(coordinatesD);
    glEnd();

    glPopAttrib();
    glPopAttrib();

    log->Delete();

    return;
}

void vtkVolumeCudaMapper::PrintSelf(ostream& os, vtkIndent indent)
{
    vtkVolumeMapper::PrintSelf(os, indent);
}
