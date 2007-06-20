
#include "vtkIGTOpenTrackerStream.h"
#include "vtkObjectFactory.h"

#include "vtkKWTkUtilities.h"
#include "vtkKWApplication.h"
#include "vtkCommand.h"
#include <OpenTracker/OpenTracker.h>
#include <OpenTracker/dllinclude.h>
#include <OpenTracker/input/SlicerNTModule.h>
#include <OpenTracker/core/Configurator.h>
//#include <OpenTracker/types/Image.h>


#include <vtksys/SystemTools.hxx>
#include "vtkCallbackCommand.h"


vtkStandardNewMacro(vtkIGTOpenTrackerStream);
vtkCxxRevisionMacro(vtkIGTOpenTrackerStream, "$Revision: 1.0 $");

vtkIGTOpenTrackerStream::vtkIGTOpenTrackerStream()
{
    this->Speed = 0;
    this->StartTimer = 0;
    this->LocatorNormalTransform = vtkTransform::New();
    this->LocatorNormalTransform_cb2 = vtkTransform::New();
    this->LocatorMatrix = vtkMatrix4x4::New();
    this->LocatorMatrix_cb2 = vtkMatrix4x4::New();// Identity
    /*
    this->LocatorNormalTransform = vtkTransform::New();
    this->LocatorMatrix = vtkMatrix4x4::New(); // Identity
    */
    this->RealtimeXsize = 0;
    this->RealtimeYsize = 0;
    this->RealtimeImageSerial = 0;

    //this->RealtimeImageData = Image::Image();
    this->RegMatrix = NULL;
    this->RegMatrix_cb2 = NULL;
    this->context = NULL;

}


vtkIGTOpenTrackerStream::~vtkIGTOpenTrackerStream()
{
    this->LocatorNormalTransform->Delete();
    this->LocatorMatrix->Delete();

    this->LocatorNormalTransform_cb2->Delete();
    this->LocatorMatrix_cb2->Delete();

    if (this->context)
    {
        delete this->context;
    }

}

void vtkIGTOpenTrackerStream::Init(const char *configFile)
{
    fprintf(stderr,"config file: %s\n",configFile);

    //OT_REGISTER_MODULE(SlicerNTModule,NULL);

    addSPLModules();

    this->context = new Context(1); 
    // get callback module from the context
    CallbackModule * callbackMod = (CallbackModule *)context->getModule("CallbackConfig");
    

    // parse the configuration file
    context->parseConfiguration(configFile);  
    context->start();

    // sets the callback function
    // if we use NaviTrack (not opentracker), use this function
    // callbackMod->setCallback( "cb1", (OTCallbackFunction*)&callbackF ,this);    
#ifdef OT_VERSION_20

    callbackMod->setCallback( "cb1", (OTCallbackFunction*)&callbackF ,this);
    callbackMod->setCallback( "cb2", (OTCallbackFunction*)&callbackF_cb2 ,this);

#endif
#ifdef OT_VERSION_13
    callbackMod->setCallback( "cb1", (CallbackFunction*)&callbackF ,this);
    callbackMod->setCallback( "cb2", (CallbackFunction*)&callbackF_cb2 ,this);
#endif

}
//this second callbackF is for receiving Orientation and Position data from a Robot with a Needle device, philip


void vtkIGTOpenTrackerStream::callbackF_cb2(Node&, Event &event, void *data_cb2)
{


    float position_cb2[3];
    float orientation_cb2[4];
    float norm_cb2[3];
    float transnorm_cb2[3];
    int j;    
    
    vtkIGTOpenTrackerStream *VOT_cb2 =(vtkIGTOpenTrackerStream *)data_cb2;

        // the original values are in the unit of meters

    position_cb2[0]=(float)(event.getPosition())[0] * VOT_cb2->MultiFactor; 
    position_cb2[1]=(float)(event.getPosition())[1] * VOT_cb2->MultiFactor;
    position_cb2[2]=(float)(event.getPosition())[2] * VOT_cb2->MultiFactor;

    orientation_cb2[0]=(float)(event.getOrientation())[0];
    orientation_cb2[1]=(float)(event.getOrientation())[1];
    orientation_cb2[2]=(float)(event.getOrientation())[2];
    orientation_cb2[3]=(float)(event.getOrientation())[3];
    
    VOT_cb2->quaternion2xyz_cb2(orientation_cb2, norm_cb2, transnorm_cb2);

     

    // Apply the transform matrix 
    // to the postion, norm and transnorm
    if (VOT_cb2->RegMatrix)
        VOT_cb2->ApplyTransform_cb2(position_cb2, norm_cb2, transnorm_cb2);

    for (j=0; j<3; j++) {
        VOT_cb2->LocatorMatrix->SetElement(j,0,position_cb2[j]);
    }


    for (j=0; j<3; j++) {
        VOT_cb2->LocatorMatrix->SetElement(j,1,norm_cb2[j]);
    }

    for (j=0; j<3; j++) {
        VOT_cb2->LocatorMatrix->SetElement(j,2,transnorm_cb2[j]);
    }

    for (j=0; j<3; j++) {
      VOT_cb2->LocatorMatrix->SetElement(j,3,0);
    }

    for (j=0; j<3; j++) {
        VOT_cb2->LocatorMatrix->SetElement(3,j,0);
    }

    VOT_cb2->LocatorMatrix->SetElement(3,3,1);


}

void vtkIGTOpenTrackerStream::callbackF(Node&, Event &event, void *data)
{

 
  
    float position[3];
    float orientation[4];
    float norm[3];
    float transnorm[3];
    int j;    
    
    vtkIGTOpenTrackerStream *VOT=(vtkIGTOpenTrackerStream *)data;

    // the original values are in the unit of meters

    position[0]=(float)(event.getPosition())[0] * VOT->MultiFactor; 
    position[1]=(float)(event.getPosition())[1] * VOT->MultiFactor;
    position[2]=(float)(event.getPosition())[2] * VOT->MultiFactor;

    orientation[0]=(float)(event.getOrientation())[0];
    orientation[1]=(float)(event.getOrientation())[1];
    orientation[2]=(float)(event.getOrientation())[2];
    orientation[3]=(float)(event.getOrientation())[3];

     VOT->quaternion2xyz(orientation, norm, transnorm);

        VOT->RealtimeXsize=(int)event.getAttribute(std::string("xdim"),0);
      VOT->RealtimeYsize=(int)event.getAttribute(std::string("ydim"),0);
      // int Xsize=(int)event.getAttribute(std::string("xdim"),0);
      // int Ysize=(int)event.getAttribute(std::string("ydim"),0);

      //  cout << "controll Callback x = " << Xsize << ", y = " << Ysize << endl;

     if(event.hasAttribute("image")) {
        VOT->RealtimeImageSerial = (VOT->RealtimeImageSerial + 1) % 32768;

     VOT->RealtimeImageData=(Image)event.getAttribute((Image*)NULL,"image");
       cout << "image size is " << VOT->RealtimeImageData.size() << endl;
     }

    // Apply the transform matrix 
    // to the postion, norm and transnorm
    if (VOT->RegMatrix)
        VOT->ApplyTransform(position, norm, transnorm);

    for (j=0; j<3; j++) {
        VOT->LocatorMatrix->SetElement(j,0,position[j]);
    }


    for (j=0; j<3; j++) {
        VOT->LocatorMatrix->SetElement(j,1,norm[j]);
    }

    for (j=0; j<3; j++) {
        VOT->LocatorMatrix->SetElement(j,2,transnorm[j]);
    }

    for (j=0; j<3; j++) {
        VOT->LocatorMatrix->SetElement(j,3,0);
    }

    for (j=0; j<3; j++) {
        VOT->LocatorMatrix->SetElement(3,j,0);
    }

    VOT->LocatorMatrix->SetElement(3,3,1);
}



void vtkIGTOpenTrackerStream::StopPolling()
{
    context->close();
}



void vtkIGTOpenTrackerStream::PollRealtime()
{
  if (context) {
    // cout <<"PollRealtime()"<< endlq;
    context->loopOnce();
  }
}



void vtkIGTOpenTrackerStream::PrintSelf(ostream& os, vtkIndent indent)
{


}


void vtkIGTOpenTrackerStream::quaternion2xyz(float* orientation, float *normal, float *transnormal) 
{
    float q0, qx, qy, qz;

    q0 = orientation[3];
    qx = orientation[0];
    qy = orientation[1];
    qz = orientation[2]; 

    transnormal[0] = 1-2*qy*qy-2*qz*qz;
    transnormal[1] = 2*qx*qy+2*qz*q0;
    transnormal[2] = 2*qx*qz-2*qy*q0;

    normal[0] = 2*qx*qz+2*qy*q0;
    normal[1] = 2*qy*qz-2*qx*q0;
    normal[2] = 1-2*qx*qx-2*qy*qy;
}


void vtkIGTOpenTrackerStream::quaternion2xyz_cb2(float* orientation_cb2, float *normal_cb2, float *transnormal_cb2) 
{
    float q0, qx, qy, qz;

    q0 = orientation_cb2[3];
    qx = orientation_cb2[0];
    qy = orientation_cb2[1];
    qz = orientation_cb2[2]; 

    transnormal_cb2[0] = 1-2*qy*qy-2*qz*qz;
    transnormal_cb2[1] = 2*qx*qy+2*qz*q0;
    transnormal_cb2[2] = 2*qx*qz-2*qy*q0;

    normal_cb2[0] = 2*qx*qz+2*qy*q0;
    normal_cb2[1] = 2*qy*qz-2*qx*q0;
    normal_cb2[2] = 1-2*qx*qx-2*qy*qy;
}




void vtkIGTOpenTrackerStream::SetLocatorTransforms()
{
    // Get locator matrix
    float p[3], n[3], t[3], c[3];
    p[0] = this->LocatorMatrix->GetElement(0, 0);
    p[1] = this->LocatorMatrix->GetElement(1, 0);
    p[2] = this->LocatorMatrix->GetElement(2, 0);
    n[0] = this->LocatorMatrix->GetElement(0, 1);
    n[1] = this->LocatorMatrix->GetElement(1, 1);
    n[2] = this->LocatorMatrix->GetElement(2, 1);
    t[0] = this->LocatorMatrix->GetElement(0, 2);
    t[1] = this->LocatorMatrix->GetElement(1, 2);
    t[2] = this->LocatorMatrix->GetElement(2, 2);


    // Ensure N, T orthogonal:
    //    C = N x T
    //    T = C x N
    this->Cross(c, n, t);
    this->Cross(t, c, n);

    // Ensure vectors are normalized
    this->Normalize(n);
    this->Normalize(t);
    this->Normalize(c); 


    /*
    # Find transform, N, that brings the locator coordinate frame 
    # into the scanner frame.  Then invert N to M and set it to the locator's
    # userMatrix to position the locator within the world space.
    #
    # 1.) Concatenate a translation, T, TO the origin which is (-x,-y,-z)
    #     where the locator's position is (x,y,z).
    # 2.) Concatenate the R matrix.  If the locator's reference frame has
    #     axis Ux, Uy, Uz, then Ux is the TOP ROW of R, Uy is the second, etc.
    # 3.) Translate the cylinder so its tip is at the origin instead
    #     of the center of its tube.  Call this matrix C.
    # Then: N = C*R*T, M = Inv(N)
    #
    # (See page 419 and 429 of "Computer Graphics", Hearn & Baker, 1997,
    #  ISBN 0-13-530924-7)
    # 
    # The alternative approach used here is to find the transform, M, that
    # moves the scanner coordinate frame to the locator's.  
    # 
    # 1.) Translate the cylinder so its tip is at the origin instead
    #     of the center of its tube.  Call this matrix C.
    # 2.) Concatenate the R matrix.  If the locator's reference frame has
    #     axis Ux, Uy, Uz, then Ux is the LEFT COL of R, Uy is the second,etc.
    # 3.) Concatenate a translation, T, FROM the origin which is (x,y,z)
    #     where the locator's position is (x,y,z).
    # Then: M = T*R*C
    */
    vtkMatrix4x4 *locator_matrix = vtkMatrix4x4::New();
    vtkTransform *locator_transform = vtkTransform::New();

    // Locator's offset: p[0], p[1], p[2]
    float x0 = p[0];
    float y0 = p[1];
    float z0 = p[2];


    // Locator's coordinate axis:
    // Ux = T
    float Uxx = t[0];
    float Uxy = t[1];
    float Uxz = t[2];

    // Uy = -N
    float Uyx = -n[0];
    float Uyy = -n[1];
    float Uyz = -n[2];

    // Uz = Ux x Uy
    float Uzx = Uxy*Uyz - Uyy*Uxz;
    float Uzy = Uyx*Uxz - Uxx*Uyz;
    float Uzz = Uxx*Uyy - Uyx*Uxy;

    // Ux
    locator_matrix->SetElement(0, 0, Uxx);
    locator_matrix->SetElement(1, 0, Uxy);
    locator_matrix->SetElement(2, 0, Uxz);
    locator_matrix->SetElement(3, 0, 0);
    // Uy
    locator_matrix->SetElement(0, 1, Uyx);
    locator_matrix->SetElement(1, 1, Uyy);
    locator_matrix->SetElement(2, 1, Uyz);
    locator_matrix->SetElement(3, 1, 0);
    // Uz
    locator_matrix->SetElement(0, 2, Uzx);
    locator_matrix->SetElement(1, 2, Uzy);
    locator_matrix->SetElement(2, 2, Uzz);
    locator_matrix->SetElement(3, 2, 0);
    // Bottom row
    locator_matrix->SetElement(0, 3, 0);
    locator_matrix->SetElement(1, 3, 0);
    locator_matrix->SetElement(2, 3, 0);
    locator_matrix->SetElement(3, 3, 1);

    // Set the vtkTransform to PostMultiply so a concatenated matrix, C,
    // is multiplied by the existing matrix, M: C*M (not M*C)
    locator_transform->PostMultiply();
    // M = T*R*C

    
    // NORMAL PART

    locator_transform->Identity();
    // C:
    locator_transform->Translate(0, (100 / 2.0), 0);
    // R:
    locator_transform->Concatenate(locator_matrix);
    // T:
    locator_transform->Translate(x0, y0, z0);

    this->LocatorNormalTransform->DeepCopy(locator_transform);
   

    locator_matrix->Delete();
    locator_transform->Delete();


}



void vtkIGTOpenTrackerStream::Normalize(float *a)
{
    float d;
    d = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);

    if (d == 0.0) return;

    a[0] = a[0] / d;
    a[1] = a[1] / d;
    a[2] = a[2] / d;
}



// a = b x c
void vtkIGTOpenTrackerStream::Cross(float *a, float *b, float *c)
{
    a[0] = b[1]*c[2] - c[1]*b[2];
    a[1] = c[0]*b[2] - b[0]*c[2];
    a[2] = b[0]*c[1] - c[0]*b[1];
}



void vtkIGTOpenTrackerStream::ApplyTransform(float *position, float *norm, float *transnorm)
{
    // Transform position, norm and transnorm
    // ---------------------------------------------------------
    float p[4];
    float n[4];
    float tn[4];

    for (int i = 0; i < 3; i++)
    {
        p[i] = position[i];
        n[i] = norm[i];
        tn[i] = transnorm[i];
    }
    p[3] = 1;     // translation affects a poistion
    n[3] = 0;     // translation doesn't affect an orientation
    tn[3] = 0;    // translation doesn't affect an orientation

    this->RegMatrix->MultiplyPoint(p, p);    // transform a position
    this->RegMatrix->MultiplyPoint(n, n);    // transform an orientation
    this->RegMatrix->MultiplyPoint(tn, tn);  // transform an orientation

    for (int i = 0; i < 3; i++)
    {
        position[i] = p[i];
        norm[i] = n[i];
        transnorm[i] = tn[i];
    }
}




void vtkIGTOpenTrackerStream::ApplyTransform_cb2(float *position_cb2, float *norm_cb2, float *transnorm_cb2)
{
    // Transform position, norm and transnorm
    // ---------------------------------------------------------
    float p[4];
    float n[4];
    float tn[4];

    for (int i = 0; i < 3; i++)
    {
        p[i] = position_cb2[i];
        n[i] = norm_cb2[i];
        tn[i] = transnorm_cb2[i];
    }
    p[3] = 1;     // translation affects a poistion
    n[3] = 0;     // translation doesn't affect an orientation
    tn[3] = 0;    // translation doesn't affect an orientation

    this->RegMatrix_cb2->MultiplyPoint(p, p);    // transform a position
    this->RegMatrix_cb2->MultiplyPoint(n, n);    // transform an orientation
    this->RegMatrix_cb2->MultiplyPoint(tn, tn);  // transform an orientation

    for (int i = 0; i < 3; i++)
    {
        position_cb2[i] = p[i];
        norm_cb2[i] = n[i];
        transnorm_cb2[i] = tn[i];
    }
}




void vtkIGTOpenTrackerStream::ProcessTimerEvents()
{
    if (this->StartTimer)
      {   
        //   cout << "vtkIGT=====================StartTimer " << endl;
   
      this->PollRealtime();
      // cout << "vtkIGT=====================PollRealtime() " << endl;
        this->InvokeEvent (vtkCommand::ModifiedEvent);
        // cout << "vtkIGT=====================InvokeEvent(vtkCommand...) " << endl;
        vtkKWTkUtilities::CreateTimerHandler(vtkKWApplication::GetMainInterp(), 
               100, this, "ProcessTimerEvents");  // RSierra 3/8/07 The integer defines the update rate. On my laptop there is no differenct in performance (i.e. the CPU load is minimal and approx. 10% for update of 2). Is the value equivalent to ms?        
    
        // cout << "vtkIGT=====================CreatTimerHandler " << endl;
   } 
   else
    {
        this->StopPolling();
    }
}

void vtkIGTOpenTrackerStream::SetTracker(std::vector<float> pos,std::vector<float> quat)
{
#if defined(OT_VERSION_20) || defined(OT_VERSION_13)
  SlicerNTModule * module = (SlicerNTModule *)context->getModule("SlicerConfig");
  module->SetTracker(pos,quat);
#endif

}

void vtkIGTOpenTrackerStream::SetOpenTrackerforScannerControll(std::vector<std::string> scancommandkeys, std::vector<std::string> scancommandvalue)
{
  #if defined(OT_VERSION_20) || defined(OT_VERSION_13)

  SlicerNTModule * module = (SlicerNTModule *)context->getModule("SlicerConfig");
  
  module->SetOpenTrackerforScannerControll(scancommandkeys, scancommandvalue);
 #endif
}


void vtkIGTOpenTrackerStream::SetOpenTrackerforBRPDataFlowValveFilter(std::vector<std::string> filtercommandkeys, std::vector<std::string> filtercommandvalue)
{
  #if defined(OT_VERSION_20) || defined(OT_VERSION_13)

  SlicerNTModule * module = (SlicerNTModule *)context->getModule("SlicerConfig");
  
  module->SetOpenTrackerforBRPDataFlowValveFilter(filtercommandkeys, filtercommandvalue);
 #endif
}



void vtkIGTOpenTrackerStream::SetOrientationforRobot(float xsendrobotcoords, float ysendrobotcoords, float zsendrobotcoords, std::vector<float> sendrobotcoordsvector, std::string robotcommandvalue,std::string robotcommandkey)
{
  #if defined(OT_VERSION_20) || defined(OT_VERSION_13)
  cout<<"opentrackerstream";
  SlicerNTModule * module = (SlicerNTModule *)context->getModule("SlicerConfig");
  
  module->SetOrientationforRobot(xsendrobotcoords, ysendrobotcoords, zsendrobotcoords, sendrobotcoordsvector,robotcommandvalue, robotcommandkey);
 #endif
}

/*
void vtkIGTOpenTrackerStream::GetSizeforRealtimeImaging(int* xsizevalueRI, int* ysizevalueRI)
{
 
  *xsizevalueRI = RealtimeXsize;
  *ysizevalueRI = RealtimeYsize;
 
}

void vtkIGTOpenTrackerStream::GetImageDataforRealtimeImaging(Image* ImageDataRI)
{
 
   *ImageDataRI = RealtimeImageData;
 
}
*/

void vtkIGTOpenTrackerStream::GetRealtimeImage(int* serial, vtkImageData* image)
{
  //std::cerr << "Serial : " << this->RealtimeImageSerial << ", " <<  *serial << std::endl;
  //std::cerr << "(xsize, ysize) = (" << RealtimeXsize << ", " << RealtimeYsize << ")" << std::endl;

    if (*serial != this->RealtimeImageSerial)
    {
      std::cerr << "Serial : " << this->RealtimeImageSerial << ", " <<  *serial << std::endl;
        *serial = this->RealtimeImageSerial;
        if (image && RealtimeImageData.size() > 0)
        {
            image->SetDimensions(RealtimeXsize, RealtimeYsize, 1);
            //image->SetExtent( xmin, xmax, ymin, ymax, zmin, zmax );
            image->SetExtent(0, RealtimeXsize-1, 0, RealtimeYsize-1, 0, 0 );
            image->SetNumberOfScalarComponents( 1 );
            image->SetOrigin( 0, 0, 0 );
            image->SetSpacing( 1, 1, 10 );
            image->SetScalarTypeToShort();
            image->AllocateScalars();

            short* dest = (short*) image->GetScalarPointer();
           if (dest) {
             memcpy(dest, RealtimeImageData.image_ptr, RealtimeImageData.size());
             image->Update();
           }
        }
    }
    else
    {
    }
}
