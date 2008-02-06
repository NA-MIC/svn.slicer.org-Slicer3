#ifndef CUDAPPSUPPORT_H_
#define CUDAPPSUPPORT_H_

#include "CudappSupportModule.h"
#include <vector>

namespace Cudapp
{
    class Device;
    class CUDA_SUPPORT_EXPORT Support
    {
    public:
        typedef std::vector<Device*> DeviceList;

        Support();
        virtual ~Support();

        bool IsSupported() { return (this->GetDeviceCount() > 0); }
        bool IsSupported(const char* cudaVersion);

        //BTX
        int GetDeviceCount() const { return this->Devices.size(); }        
        const DeviceList GetDevices() { return this->Devices; }
        Device* operator[](int deviceNumber) const { return this->Devices[deviceNumber]; }
        //ETX

        virtual void PrintSelf(std::ostream&  os) const;

    protected:

        int CheckSupportedCudaVersion();
        //BTX
        DeviceList Devices;
        //ETX
    };
    inline std::ostream& operator<<(std::ostream& os, const Support& in){
        in.PrintSelf(os);
        return os; 
    }

}
#endif /*CUDAPPSUPPORT_H_*/
