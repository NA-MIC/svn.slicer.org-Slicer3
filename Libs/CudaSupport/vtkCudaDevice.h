#ifndef VTKCUDADEVICE_H_
#define VTKCUDADEVICE_H_

#include "vtkObject.h"
#include "vtkCudaSupportModule.h"

#include "driver_types.h"

class VTK_CUDASUPPORT_EXPORT vtkCudaDevice : public vtkObject
{
  vtkTypeRevisionMacro(vtkCudaDevice, vtkObject);
  
  static vtkCudaDevice* New();

/// Device Information
  vtkGetMacro(Initialized, bool);
  
  void SetDeviceNumber(int deviceNumber);
  vtkGetMacro(DeviceNumber, int);
  
  /** @returns the name of the CUDA device */
  const char* GetName() const { return DeviceProp.name; }
  /** @returns the global memory size of the CUDA device */
  size_t GetTotalGlobalMem() const { return DeviceProp.totalGlobalMem; }
  /** @returns the shared memory per block of the CUDA device */
  size_t GetSharedMemPerBlock() const { return DeviceProp. sharedMemPerBlock; }
  /** @reuturns the registers per block of the CUDA device */
  int GetRegsPerBlock() const { return DeviceProp.regsPerBlock; }
  /** @returns the wrap size of the CUDA device */
  int GetWrapSize() const { return DeviceProp.warpSize; }
  /** @returns the memory pitch of the CUDA device */
  size_t GetMemPitch() const { return DeviceProp.memPitch; }
  /** @returns the maximum numbers of threads that can be run in parallel on this CUDA device */
  int GetMaxThreadsPerBlock() const { return DeviceProp.maxThreadsPerBlock; }
  /** @returns the maximum sizes of each dimension of a block for the CUDA device as int[3] */
  const int* GetMaxThreadsDim() const { return DeviceProp.maxThreadsDim; } // [3]
  /** @returns the maximum sizes of each dimension of a grid for this CUDA device as int[3] */
  const int* GetMaxGridSize() const { return DeviceProp.maxGridSize; }  // [3]
  /** @returns the tatal amount (in bytes) of constant memory of this CUDA device */
  size_t GetTotalConstMem() const { return DeviceProp.totalConstMem; }
  /** @reuturns The major revision number of the CUDA device */
  int GetMajor() const { return DeviceProp.major; }
  /** @returns the minor revision number of the CUDA device */
  int GetMinor() const { return DeviceProp.minor; }
  /** @reuturns the clock rate in kiloherz of the CUDA device */
  int GetClockRate() const { return DeviceProp.clockRate; }
  /** @returns the texture alignment of the CUDA device */
  size_t GetTextureAlignment() const { return DeviceProp.textureAlignment; }
  /** @reuturns the entire device properties of this CUDA device as cudaGetDeviceProperties(this->GetDeviceNumber()) would return it */
  const cudaDeviceProp& GetCudaDeviceProperty() const { return this->DeviceProp; }
  
  
  
  void MakeActive();
  void SynchronizeThread();
  void ExitThread();

  /// Memory Management
  bool AllocateMemory();
        
  void PrintSelf(ostream& os, vtkIndent indent);
        
  protected:
    vtkCudaDevice();
    virtual ~vtkCudaDevice();

  void LoadDeviceProperties();

    bool Initialized;
    int DeviceNumber;
    cudaDeviceProp DeviceProp;    
};

#endif /*VTKCUDADEVICE_H_*/
