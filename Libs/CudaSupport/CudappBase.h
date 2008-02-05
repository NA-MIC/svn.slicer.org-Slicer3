#ifndef CUDAPPBASE_H_
#define CUDAPPBASE_H_

#include <ostream>
#include "driver_types.h"
#include "CudappSupportModule.h"

/// THIS IS A STATIC CLASS USED FOR BASIC CUDA FUNCTIONALITY!!
class CUDA_SUPPORT_EXPORT CudappBase 
{
public:
    //BTX
    typedef enum {
        Success,
        NotReadyError,
        InvalidValueError,
    } State;
    //ETX

    static CudappBase* New();

    static cudaError_t GetLastError();
    static const char* GetLastErrorString();
    static const char* GetErrorString(cudaError_t error);
    static void PrintError(cudaError_t error);

private:
    virtual ~CudappBase();
    CudappBase();
    CudappBase(const CudaBase&) {}
};

#endif /*CUDAPPBASE_H_*/
