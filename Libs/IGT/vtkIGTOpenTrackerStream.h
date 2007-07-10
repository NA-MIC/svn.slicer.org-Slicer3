 // .NAME vtkIGTOpenTrackerStream - Central registry to provide control and I/O for
//  trackers and imagers
// .SECTION Description
// vtkIGTOpenTrackerStream registers arbitary number of trackers and imagers, created MRML nodes in the MRML secene. Designed and Coded by Nobuhiko Hata and Haiying Liu, Jan 12, 2007 @ NA-MIC All Hands Meeting, Salt Lake City, UT

#ifndef IGTOPENTRACKERSTREAM_H
#define IGTOPENTRACKERSTREAM_H


#include <string>
#include <vector>

#include "vtkIGTWin32Header.h" 
#include "vtkObject.h"

#include "vtkMatrix4x4.h"
#include "vtkTransform.h"

#include "vtkImageData.h"


#ifdef OT_VERSION_20
#include "OpenTracker/OpenTracker.h"
#include "OpenTracker/common/CallbackModule.h"
#include <OpenTracker/types/Image.h>
#endif
#ifdef OT_VERSION_13
#include "OpenTracker.h"
#include "common/CallbackModule.h"
#endif

using namespace ot;


class VTK_IGT_EXPORT vtkIGTOpenTrackerStream : public vtkObject
{
public:


    static vtkIGTOpenTrackerStream *New();
    vtkTypeRevisionMacro(vtkIGTOpenTrackerStream,vtkObject);
    void PrintSelf(ostream& os, vtkIndent indent);

    vtkSetMacro(Speed,int);
    vtkSetMacro(MultiFactor,float);

    vtkSetMacro(StartTimer,int);

    vtkSetObjectMacro(RegMatrix,vtkMatrix4x4);
    vtkGetObjectMacro(RegMatrix,vtkMatrix4x4);

   

     vtkGetObjectMacro(LocatorMatrix,vtkMatrix4x4);
   
     vtkSetMacro(position_cb2_FS0,float);
     vtkSetMacro(position_cb2_FS1,float);
     vtkSetMacro(position_cb2_FS2,float);

     vtkSetMacro(orientation_cb2_FS0,float);
     vtkSetMacro(orientation_cb2_FS1,float);
     vtkSetMacro(orientation_cb2_FS2,float);
     vtkSetMacro(orientation_cb2_FS3,float);


    vtkSetMacro(RealtimeXsize,int);
    vtkSetMacro(RealtimeYsize,int);
    // vtkSetMacro(RealtimeImageData,Image);
    
   
    vtkGetObjectMacro(LocatorNormalTransform,vtkTransform);
       
    /**
     * Constructor
     **/
    vtkIGTOpenTrackerStream();


    //Description:
    //Destructor

    virtual ~vtkIGTOpenTrackerStream ( );
  

    void Init(const char *configFile);
    void StopPolling();
    void PollRealtime();
    void SetLocatorTransforms();
   

    void ProcessTimerEvents();

    
     
    static void callbackF(Node&, Event&, void *data);
     static void callbackF_cb2(Node&, Event&, void *data_cb2);
    //BTX
    void SetTracker(std::vector<float> pos,std::vector<float> quat);
    //ETX
    //BTX
    void SetOpenTrackerforScannerControll(std::vector<std::string> scancommandkeys, std::vector<std::string> scancommandvalue);
    //ETX
    //BTX
    void SetOpenTrackerforBRPDataFlowValveFilter(std::vector<std::string> filtercommandkeys, std::vector<std::string> filtercommandvalue);
    //ETX
    //BTX
    void SetOrientationforRobot(float xsendrobotcoords, float ysendrobotcoords, float zsendrobotcoords, std::vector<float> sendrobotcoordsvector, std::string robotcommandvalue,std::string robotcommandkey);
    //ETX    
    //BTX
      void GetRealtimeImage(int*, vtkImageData* image);
      //ETX

      //BTX
      //  void GetCoordsOrientforScanner(std::vector<float> OrientationForScanner, std::vector<float> PositionForScanner);
      void GetCoordsOrientforScanner(float* OrientationForScanner,float* PositionForScanner);
     
 //ETX

private:

    int Speed;
    int StartTimer;
    float MultiFactor;

    vtkMatrix4x4 *LocatorMatrix;
    vtkMatrix4x4 *LocatorMatrix_cb2;
    vtkMatrix4x4 *RegMatrix;
    vtkMatrix4x4 *RegMatrix_cb2;
    vtkTransform *LocatorNormalTransform;
    vtkTransform *LocatorNormalTransform_cb2;
                                          
    Context *context;
    
    float position_cb2_FS0;
    float position_cb2_FS1;
    float position_cb2_FS2;

    float orientation_cb2_FS0;
    float orientation_cb2_FS1;
    float orientation_cb2_FS2;
    float orientation_cb2_FS3;
    // float orientation_cb2_FS[4];
    
    float position_cb2[3];
    float orientation_cb2[4];
    

    int RealtimeXsize;
    int RealtimeYsize;

  
    Image RealtimeImageData;
     int RealtimeImageSerial;


    void Normalize(float *a);
    void Cross(float *a, float *b, float *c);
    void ApplyTransform(float *position, float *norm, float *transnorm);
    void ApplyTransform_cb2(float *position_cb2, float *norm_cb2, float *transnorm_cb2);
    void CloseConnection();

    void quaternion2xyz(float* orientation, float *normal, float *transnormal); 
    void quaternion2xyz_cb2(float* orientation_cb2, float *normal_cb2, float *transnormal_cb2); 

    
    
};

#endif // IGTOPENTRACKERSTREAM_H

