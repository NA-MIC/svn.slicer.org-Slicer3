/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkTrackFileData.cxx,v $
  Language:  C++
  Date:      $Date: 2008/12/29 18:39:45 $
  Version:   $Revision: 1.15 $

  =========================================================================*/
#include "vtkTrackFileData.h"
#include "vtkObjectFactory.h"
#include "vtkStructuredPoints.h"

#include "itkArchetypeSeriesFileNames.h"

#include <iostream>
#include <sstream>
#include <algorithm>

using namespace std;

vtkCxxRevisionMacro(vtkTrackFileData, "$Revision: 1.15 $");
vtkStandardNewMacro(vtkTrackFileData);

vtkTrackFileData::vtkTrackFileData() {

  this->Reader = vtkStructuredPointsReader::New();
  this->Writer = vtkStructuredPointsWriter::New();
  this->FileNames = vtkGlobFileNames::New();

  this->Reslice = vtkImageReslice::New();

  this->TrackingSystem = vtkTrackingSystem::New();
  this->TrackingSystem->SetEulerAnglesInDegrees();

  this->ImageCast  = vtkImageCast::New();

  this->ImageLuminance = vtkImageLuminance::New();

  this->NumberSourceFiles = 0;
  this->SourceImageData = NULL;
}

vtkTrackFileData::~vtkTrackFileData() { 
  if (this->Reader != NULL) {
    this->Reader->Delete();
  }
  if (this->Writer != NULL) {
    this->Writer->Delete();
  }
  if (this->Reslice != NULL) {
    this->Reslice->Delete();
  }
  if (this->FileNames != NULL) {
    this->FileNames->Delete();
  }
  if (this->TrackingSystem != NULL) {
    this->TrackingSystem->Delete();
  }
  if (this->ImageCast != NULL) {
    this->ImageCast->Delete();
  }
  if (this->ImageLuminance != NULL) {
    this->ImageLuminance->Delete();
  }
}


void vtkTrackFileData::PrintSelf(ostream &os, vtkIndent indent) 
{
  //SuperClass::PrintSelf(os,indent);
}

int 
vtkTrackFileData::InitilizeRead(const char *filename, const char *pattern)
{
  this->FileNames->SetDirectory(filename);
  this->FileNames->AddFileNames(pattern);
  this->NumberSourceFiles = this->FileNames->GetNumberOfFileNames();
  this->SortedFileNames.clear();

  this->ImageCast->SetOutputScalarTypeToUnsignedChar();

  if (this->NumberSourceFiles == 1)
    {
    int dims[3];
    this->SortedFileNames.push_back(this->FileNames->GetNthFileName(0));
    this->Reader->SetFileName(this->FileNames->GetNthFileName(0));
    this->Reader->Update();
    this->Reader->GetOutput()->UpdateData();

    this->SourceImageData = (vtkImageData *)this->Reader->GetOutput();
    //this->Reader->ReleaseDataFlagOn();
        
    this->SourceImageData->UpdateData();
    this->ImageLuminance->SetInput(this->SourceImageData);
    this->ImageLuminance->Update();
    this->SourceImageData = this->ImageLuminance->GetOutput();
    this->SourceImageData->UpdateData();

    this->SourceImageData->GetDimensions(dims);
    return dims[2];
    }
  else 
    {
    itk::ArchetypeSeriesFileNames::Pointer fit = itk::ArchetypeSeriesFileNames::New();
    fit->SetArchetype (this->FileNames->GetNthFileName(0));
    this->SortedFileNames = fit->GetFileNames();

    return this->NumberSourceFiles;
    }
}

  
void 
vtkTrackFileData::InitilizeWrite(const char *dir, const char *prefix)
{
  OutDirectory = dir;
  OutPrefix = prefix;
  this->Writer->SetFileTypeToBinary();
}

  
int 
vtkTrackFileData::ReadStep(int index, 
                           double &timeStamp, 
                           vtkMatrix4x4 **calibMatrix, 
                           vtkMatrix4x4 **regMatrix, 
                           vtkMatrix4x4 **sensorMatrix, 
                           vtkImageData **image)
{
  vtkFieldData *field;
  vtkDoubleArray *reg = NULL;
  vtkDoubleArray *calib = NULL;
  vtkDoubleArray *sensorData = NULL;
  vtkDoubleArray *times = NULL;
  int dims[3], i, j;
  double sp[3];
  int sensorIndex = index;

  timeStamp = 0;
  *calibMatrix = NULL;
  *regMatrix = NULL;
  *sensorMatrix = NULL;
  *image = NULL;

  if (this->NumberSourceFiles > 1)
    {
    sensorIndex = 0;
    if (index >= this->NumberSourceFiles)
      {
      std::cerr << "tkTrackFileData::ReadStep() index " << index << " exceeds number of files " << this->NumberSourceFiles << "\n";
      return 0;
      }
    this->Reader->SetFileName(this->SortedFileNames[index].c_str());
    this->Reader->Update();
    this->SourceImageData = (vtkImageData *)this->Reader->GetOutput();
    this->SourceImageData->UpdateData();
    (*image) = this->SourceImageData;
    }
  else if (this->SourceImageData &&  this->NumberSourceFiles == 1)
    {
    this->SourceImageData->GetDimensions(dims);
    this->SourceImageData->GetSpacing(sp);
    this->Reslice->SetResliceAxesDirectionCosines(1,0,0,0,1,0,0,0,1);
    this->Reslice->SetOutputSpacing(sp[0],sp[1],1.0);
    this->Reslice->SetOutputDimensionality(2);

    this->Reslice->SetInput(this->SourceImageData);
    this->Reslice->SetOutputExtent(0,dims[0]-1,0,dims[1]-1,0,1);
    this->Reslice->SetResliceAxesOrigin(0,0,index);
    //this->Reslice->Update();
    this->ImageCast->SetInput(this->Reslice->GetOutput());
    this->ImageCast->Update();

    //(*image) = this->Reslice->GetOutput();

    (*image) = this->ImageCast->GetOutput();
    (*image)->UpdateData();
    }
  else
    {
    return 0;
    }

  field = this->SourceImageData->GetFieldData();
  reg = (vtkDoubleArray *)field->GetArray("Registration");
  calib = (vtkDoubleArray *)field->GetArray("Calibration");
  sensorData = (vtkDoubleArray *)field->GetArray("SensorPositions");
  times = (vtkDoubleArray *)field->GetArray("TimeStamps");

  if (reg)
    {
    *regMatrix = vtkMatrix4x4::New();
    if (reg->GetNumberOfComponents() == 6)
      {
      this->TrackingSystem->BuildMatrixFromEulerAngles(reg, *regMatrix);
      }
    else if (reg->GetNumberOfComponents() == 16)
      {
      for (i=0;i<4;i++) 
        {
        for (j=0;j<4;j++) 
          {
          (*regMatrix)->SetElement(i,j,reg->GetComponent(0,j+4*i));
          }
        }
      }
    }

  if (calib)
    {
    *calibMatrix = vtkMatrix4x4::New();
    if (calib->GetNumberOfComponents() == 6)
      {
      this->TrackingSystem->BuildMatrixFromEulerAngles(calib, *calibMatrix);
      }
    else if (calib->GetNumberOfComponents() == 16)
      {
      for (i=0;i<4;i++) 
        {
        for (j=0;j<4;j++) 
          {
          (*calibMatrix)->SetElement(i,j,calib->GetComponent(0,j+4*i));
          }
        }
      }
    }


  if (times && sensorIndex < times->GetNumberOfTuples())
    {
    timeStamp = times->GetComponent(sensorIndex,0);
    }


  if (sensorData && sensorIndex < sensorData->GetNumberOfTuples()) 
    {
    *sensorMatrix = vtkMatrix4x4::New();

    if (sensorData->GetNumberOfComponents() == 6)
      {
      double r[6];
      double mat[4][4];
      for (i=0; i<6; i++)
        {
        r[i] = sensorData->GetComponent(index, i);
        }
      this->TrackingSystem->BuildMatrixFromEulerAngles(r[0],r[1],r[2],r[3],r[4],r[5], mat);

      std::cerr << "Data: " << r[0]<<" "<<r[1]<<" "<<r[2]<<" "<<r[3]<<" "<<r[4]<<" "<<r[5] << "\n";

      for (i=0;i<4;i++) 
        {
        for (j=0;j<4;j++) 
          {
          (*sensorMatrix)->SetElement(i,j,mat[i][j]);
          }
          std::cerr << "Matrix: ["<< i << "] = "  << mat[i][0]<<" "<<mat[i][1]<<" "<<mat[i][2]<<" "<<mat[i][3]<< "\n";
        }

      }
    else if (sensorData->GetNumberOfComponents() == 16)
      {
      for (i=0;i<4;i++) 
        {
        for (j=0;j<4;j++) 
          {
          (*sensorMatrix)->SetElement(i,j,sensorData->GetComponent(0,j+4*i));
          }
        }
      }
    }



  return 1;
}

int 
vtkTrackFileData::WriteStep(int index, 
                           double timeStamp, 
                           vtkMatrix4x4 *calibMatrix, 
                           vtkMatrix4x4 *regMatrix, 
                           vtkMatrix4x4 *sensorMatrix, 
                           vtkImageData *image)
{
  if (image == NULL)
    {
    return 0;
    }

  vtkFieldData *field;

  vtkStructuredPointsWriter *writer = vtkStructuredPointsWriter::New();
  vtkDoubleArray *reg = vtkDoubleArray::New();
  vtkDoubleArray *calib = vtkDoubleArray::New();
  vtkDoubleArray *sensorData = vtkDoubleArray::New();
  vtkDoubleArray *times = vtkDoubleArray::New();

  std::stringstream ss;
  ss << index;

  std::string filename = this->OutDirectory + "/" 
        + this->OutPrefix + ss.str() + ".vtk";
  this->Writer->SetFileName(filename.c_str());

  int i,j;

  reg->SetNumberOfComponents(16);
  calib->SetNumberOfComponents(16);
  sensorData->SetNumberOfComponents(16);
  times->SetNumberOfComponents(1);

  for (i=0;i<4;i++) {
    for (j=0;j<4;j++) {
      if (regMatrix) 
        {
        reg->InsertComponent(0,j+4*i, regMatrix->GetElement(i,j));
        }
      if (calibMatrix)
        {
        calib->InsertComponent(0,j+4*i, calibMatrix->GetElement(i,j));
        }
      if (sensorMatrix)
        {
        sensorData->InsertComponent(0,j+4*i, sensorMatrix->GetElement(i,j));
        }
    }
  }

  times->InsertComponent(0, 0, timeStamp);

  reg->SetName("Registration");
  calib->SetName("Calibration");
  sensorData->SetName("SensorPositions");
  times->SetName("TimeStamps");

  field = image->GetFieldData();
  field->AddArray(reg);
  field->AddArray(calib);
  field->AddArray(sensorData);
  field->AddArray(times);
  image->UpdateData();

  this->Writer->SetInput(image);
  this->Writer->Write();

  reg->Delete();
  calib->Delete();
  sensorData->Delete();

  return 1;
}


