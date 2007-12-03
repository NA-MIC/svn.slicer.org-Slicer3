
#include "vtkIGTOpenTrackerStream.h"
#include "vtkObjectFactory.h"

#include "vtkKWTkUtilities.h"
#include "vtkKWApplication.h"
#include "vtkCommand.h"



#include <vtksys/SystemTools.hxx>
#include "vtkCallbackCommand.h"


vtkStandardNewMacro(vtkIGTOpenTrackerStream);
vtkCxxRevisionMacro(vtkIGTOpenTrackerStream, "$Revision: 1.0 $");

vtkIGTOpenTrackerStream::vtkIGTOpenTrackerStream()
{
    this->context = NULL;
}


vtkIGTOpenTrackerStream::~vtkIGTOpenTrackerStream()
{
    if (this->context)
    {
        delete this->context;
    }
}


void vtkIGTOpenTrackerStream::Init(const char *configFile)
{
    fprintf(stderr,"config file: %s\n",configFile);
    this->context = new Context(1); 
    // get callback module from the context
    CallbackModule * callbackMod = (CallbackModule *)context->getModule("CallbackConfig");

    // parse the configuration file
    context->parseConfiguration(configFile);  

    // sets the callback function
    // if we use NaviTrack (not opentracker), use this function
    // callbackMod->setCallback( "cb1", (OTCallbackFunction*)&callbackF ,this);    
    callbackMod->setCallback( "cb1", (OTCallbackFunction*)&callbackF ,this);

    // Add additional callbacks
    vtkIGTMessageAttributeSet::AttributeSetMap::iterator iter;
    for (iter = AttributeSetMap.begin(); iter != AttributeSetMap.end(); iter ++)
      {
      const char* cbname                 = iter->first.c_str();
      vtkIGTMessageAttributeSet* attrSet = iter->second;
      callbackMod->setCallback(cbname, (OTCallbackFunction*)&GenericCallback, (void*)attrSet);
      }

    context->start();
}



void vtkIGTOpenTrackerStream::callbackF(const Node&, const Event &event, void *data)
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

    VOT->QuaternionToXYZ(orientation, norm, transnorm);


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


void vtkIGTOpenTrackerStream::GenericCallback(const Node &node, const Event &event, void *data)
{

  if (!data)
    {
    return;
    }

  vtkIGTMessageAttributeSet* attrSet = (vtkIGTMessageAttributeSet*)data;
  
  // GET attributes
  //vtkIGTMessageAttributeSet::AttributeMapType& attrMap = attrSet->GetAttributeMap();
  //vtkIGTMessageAttributeSet::AttributeMapType::iterator iter;
  vtkIGTMessageAttributeSet::iterator iter;

  for (iter = attrSet->begin(); iter != attrSet->end(); iter ++)
    {
    std::string key = iter->first;
    vtkIGTMessageAttributeBase* attr = iter->second;

    if (event.hasAttribute(key))
      {

      //========== Macro for switch(attr->GetTypeID()) {} ==========
      #define CASE_GETATTRIB_TYPE(TYPE_ID, TYPE)       \
        case TYPE_ID:                                         \
          {                                                   \
          TYPE data = (TYPE)event.getAttribute<TYPE>((TYPE*)NULL, key); \
          attr->SetAttribute(&data);                          \
          }                                                   \
        break;
      //============================================================
      std::vector<float> fvec;
      switch(attr->GetTypeID())
        {
        CASE_GETATTRIB_TYPE(vtkIGTMessageAttributeSet::TYPE_BOOL,           bool               );//(bool*)          false);
        CASE_GETATTRIB_TYPE(vtkIGTMessageAttributeSet::TYPE_CHAR,           char               );//(char*)          0);
        CASE_GETATTRIB_TYPE(vtkIGTMessageAttributeSet::TYPE_SIGNED_CHAR,    signed char        );//(signed char*)   0);
        CASE_GETATTRIB_TYPE(vtkIGTMessageAttributeSet::TYPE_UNSIGNED_CHAR,  unsigned char      );//(unsigned char*) 0);
        CASE_GETATTRIB_TYPE(vtkIGTMessageAttributeSet::TYPE_INT,            int                );//(int*)           0);
        CASE_GETATTRIB_TYPE(vtkIGTMessageAttributeSet::TYPE_LONG,           long               );//(long*)          0);
        CASE_GETATTRIB_TYPE(vtkIGTMessageAttributeSet::TYPE_SHORT,          short              );//(short*)         0);
        CASE_GETATTRIB_TYPE(vtkIGTMessageAttributeSet::TYPE_UNSIGNED_INT,   unsigned int       );//(unsigned int*)  0);
        CASE_GETATTRIB_TYPE(vtkIGTMessageAttributeSet::TYPE_UNSIGNED_LONG,  unsigned long      );//(unsigned long*) 0);
        CASE_GETATTRIB_TYPE(vtkIGTMessageAttributeSet::TYPE_UNSIGNED_SHORT, unsigned short     );//(unsigned short*)0);
        CASE_GETATTRIB_TYPE(vtkIGTMessageAttributeSet::TYPE_DOUBLE,         double             );//(double*)        0.0);
        CASE_GETATTRIB_TYPE(vtkIGTMessageAttributeSet::TYPE_LONG_DOUBLE,    long double        );//(long double*)   0.0);
        CASE_GETATTRIB_TYPE(vtkIGTMessageAttributeSet::TYPE_FLOAT,          float              );//(float*)         0.0);
        CASE_GETATTRIB_TYPE(vtkIGTMessageAttributeSet::TYPE_STRING,         std::string        );//(char*)NULL);
        CASE_GETATTRIB_TYPE(vtkIGTMessageAttributeSet::TYPE_VECTOR_FLOAT,   std::vector<float> );//fvec);
        case vtkIGTMessageAttributeSet::TYPE_VTK_IMAGE_DATA:// vtkImageData
          {
          }
          break;

        /*
        case vtkIGTMessageAttributeSet::TYPE_BOOL:          // bool
          {
          bool data = (bool)event.getAttribute(key, false);
          attr->SetAttribute(&data);
          }
          break;
        case vtkIGTMessageAttributeSet::TYPE_CHAR:          // char
          {
          char data = (char)event.getAttribute(key, 0);
          attr->SetAttribute(&data);
          }
          break;
        case vtkIGTMessageAttributeSet::TYPE_SIGNED_CHAR:   // signed char
          {
          signed char data = (signed char)event.getAttribute(key, 0);
          attr->SetAttribute(&data);
          }
          break;
        case vtkIGTMessageAttributeSet::TYPE_UNSIGNED_CHAR: // unsigned char
          {
          unsigned char data = (unsigned char)event.getAttribute(key, 0);
          attr->SetAttribute(&data);
          }
          break;
        case vtkIGTMessageAttributeSet::TYPE_INT:           // int
          {
          int data = (int)event.getAttribute(key, 0);
          attr->SetAttribute(&data);
          }
          break;
        case vtkIGTMessageAttributeSet::TYPE_LONG:          // long
          {
          long data = (long)event.getAttribute(key, 0);
          attr->SetAttribute(&data);
          }
          break;
        case vtkIGTMessageAttributeSet::TYPE_SHORT:         // short
          {
          short data = (short)event.getAttribute(key, 0);
          attr->SetAttribute(&data);
          }
          break;
        case vtkIGTMessageAttributeSet::TYPE_UNSIGNED_INT:  // unsigned int
          {
          unsigned int data = (unsigned int)event.getAttribute(key, 0);
          attr->SetAttribute(&data);
          }
          break;
        case vtkIGTMessageAttributeSet::TYPE_UNSIGNED_LONG: // unsigned long
          {
          unsigned long data = (unsigned long)event.getAttribute(key, 0);
          attr->SetAttribute(&data);
          }
          break;
        case vtkIGTMessageAttributeSet::TYPE_UNSIGNED_SHORT:// unsigned short
          {
          unsigned short data = (unsigned short)event.getAttribute(key, 0);
          attr->SetAttribute(&data);
          }
          break;
        case vtkIGTMessageAttributeSet::TYPE_DOUBLE:        // double
          {
          double data = (double)event.getAttribute(key, 0);
          attr->SetAttribute(&data);
          }
          break;
        case vtkIGTMessageAttributeSet::TYPE_LONG_DOUBLE:   // long double
          {
          long double data = (long double)event.getAttribute(key, 0);
          attr->SetAttribute(&data);
          }
          break;
        case vtkIGTMessageAttributeSet::TYPE_FLOAT:         // float
          {
          float data = (float)event.getAttribute(key, 0);
          attr->SetAttribute(&data);
          }
          break;
        case vtkIGTMessageAttributeSet::TYPE_STRING:        // std::string
          {
          std::string data = (std::string)event.getAttribute(key, 0);
          attr->SetAttribute(&data);
          }
          break;
        case vtkIGTMessageAttributeSet::TYPE_VECTOR_FLOAT:  // std::vector<float>
          {
          std::vector<float> data = (std::vector<float>)event.getAttribute(key, 0);
          attr->SetAttribute(&data);
          }
          break;
        case vtkIGTMessageAttributeSet::TYPE_VTK_IMAGE_DATA:// vtkImageData
          {
          }
          break;
        */
        default:
          break;
        }

      }
    }
  
  vtkIGTMessageAttributeSet::MessageHandlingFunction* handler = attrSet->GetHandlerFunction();
  if (handler)
    {
    handler(attrSet);
    }

}


void vtkIGTOpenTrackerStream::AddCallback(const char* cbname,
                                          vtkIGTMessageAttributeSet::MessageHandlingFunction* func,
                                          vtkIGTMessageAttributeSet* attrSet)
{
    attrSet->SetOpenTrackerStream(this);
    attrSet->SetHandlerFunction(func);
    this->AttributeSetMap[cbname] = attrSet;
}


void vtkIGTOpenTrackerStream::StopPulling()
{
    context->close();
}



void vtkIGTOpenTrackerStream::PullRealTime()
{
    context->pushEvents();       // push event and
    context->pullEvents();       // pull event 
    context->stop();
}



void vtkIGTOpenTrackerStream::PrintSelf(ostream& os, vtkIndent indent)
{


}

