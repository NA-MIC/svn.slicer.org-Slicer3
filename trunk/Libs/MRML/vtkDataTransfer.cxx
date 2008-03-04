#include "vtkDataTransfer.h"


vtkStandardNewMacro ( vtkDataTransfer );
vtkCxxRevisionMacro ( vtkDataTransfer, "$Revision: 1.0 $" );

//----------------------------------------------------------------------------
vtkDataTransfer::vtkDataTransfer()
{
  this->SourceURI = NULL;
  this->DestinationURI = NULL;
  this->Handler = NULL;
  this->TransferStatus = vtkDataTransfer::Unspecified;
  this->TransferID = -1;
  this->TransferType = vtkDataTransfer::Initialized;
  this->TransferNodeID = NULL;
  this->Progress = 0;
}


//----------------------------------------------------------------------------
vtkDataTransfer::~vtkDataTransfer()
{
   
  this->SourceURI = NULL;
  this->DestinationURI = NULL;
  this->Handler = NULL;
  this->TransferStatus = vtkDataTransfer::Unspecified;
  this->TransferID = -1;
  this->TransferType = vtkDataTransfer::Initialized;
  this->TransferNodeID = NULL;
  this->Progress = 0;
}


//----------------------------------------------------------------------------
void vtkDataTransfer::PrintSelf(ostream& os, vtkIndent indent)
{
  Superclass::PrintSelf ( os, indent );
  os << indent << "SourceURI: " <<
    ( this->SourceURI ? this->SourceURI : "(none)") << "\n";
  os << indent << "DestinationURI: " <<
    ( this->DestinationURI ? this->DestinationURI : "(none)") << "\n";
  os << indent << "Handler: " << this->GetHandler() << "\n";
  os << indent << "TransferStatus: " << this->GetTransferStatus() << "\n";
  os << indent << "TransferID: " << this->GetTransferID() << "\n";
  os << indent << "TransferType: " << this->GetTransferType() << "\n";
  os << indent << "TransferNodeID: " << this->GetTransferNodeID() << "\n";
  os << indent << "Progress: " << this->GetProgress() << "\n";
}
