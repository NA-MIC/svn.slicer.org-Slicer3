/*=auto=========================================================================

  Portions (c) Copyright 2008 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See Doc/copyright/copyright.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: vtkEventBroker.cxx,v $
  Date:      $Date: 2006/01/06 17:56:48 $
  Version:   $Revision: 1.58 $

=========================================================================auto=*/

#include "vtkObjectFactory.h"
#include "vtkCollection.h"
#include "vtkCallbackCommand.h"
#include "vtkEventBroker.h"
#include "vtkObservation.h"

vtkEventBroker *vtkEventBroker::Instance = NULL;

vtkCxxRevisionMacro(vtkEventBroker, "$Revision: 1.9.12.1 $");

//----------------------------------------------------------------------------
// Up the reference count so it behaves like New
vtkEventBroker* vtkEventBroker::New()
{
  vtkEventBroker* ret = vtkEventBroker::GetInstance();
  ret->Register(NULL);
  return ret;
}

//----------------------------------------------------------------------------
// Return the single instance of the vtkEventBroker
vtkEventBroker* vtkEventBroker::GetInstance()
{
  if(!vtkEventBroker::Instance)
    {
    // Try the factory first
    vtkEventBroker::Instance = (vtkEventBroker*)
      vtkObjectFactory::CreateInstance("vtkEventBroker");
    // if the factory did not provide one, then create it here
    if(!vtkEventBroker::Instance)
      {
      vtkEventBroker::Instance = new vtkEventBroker;
      }
    }
  // return the instance
  return vtkEventBroker::Instance;
}

//----------------------------------------------------------------------------
vtkEventBroker::vtkEventBroker()
{
  this->CallbackCommand = vtkCallbackCommand::New();
  this->CallbackCommand
  this->EventMode = vtkEventBroker::Synchronous;
  this->EventLogging = 0;
  this->LogFileName = NULL;
  this->ScriptHandler = NULL;
}

//----------------------------------------------------------------------------
vtkEventBroker::~vtkEventBroker()
{
  // clear out all the observation records
  std::vector< vtkObservation *>::iterator iter; 
  for(iter=this->Observations.begin(); iter != this->Observations.end(); iter++)  
    { 
    this->DetachObservation (*iter);
    (*iter)->Delete();
    }

  // close the event log if needed
  if ( this->LogFile.is_open() )
    {
    this->EventLoggingOff();
    this->LogEvent( NULL );
    }
}


//----------------------------------------------------------------------------
vtkObservation *vtkEventBroker::AddObservation (
  vtkObject *subject, unsigned long event, vtkObject *observer, vtkCallbackCommand *notify)
{
  vtkObservation *observation = vtkObservation::New();
  observation->SetEventBroker( this );
  this->Observations.push_back( observation );
  observation->AssignSubject( subject );
  observation->SetEvent( event );
  observation->AssignObserver( observer );
  observation->SetCallbackCommand( notify );

  this->AttachObservation( observation );
  return (observation);
}

//----------------------------------------------------------------------------
vtkObservation *vtkEventBroker::AddObservation (
  vtkObject *subject, const char *event, const char *script)
{
  vtkObservation *observation = vtkObservation::New();
  observation->SetEventBroker( this );
  this->Observations.push_back( observation );
  observation->AssignSubject( subject );

  // figure out event either as a predefined string, or
  // as an ascii number
  unsigned long eventID = vtkCommand::GetEventIdFromString( event );
  if ( eventID == vtkCommand::NoEvent )
    {
    eventID = static_cast<unsigned long> (atoi( event ));
    }
  observation->SetEvent( eventID );
  observation->SetScript( script );

  this->AttachObservation( observation );
  return (observation);
}

//----------------------------------------------------------------------------
void vtkEventBroker::AttachObservation ( vtkObservation *observation )
{
  //
  // set up to be notified if either of the objects is delete
  // and when the specified event takes place
  // (remove any old notifications first)
  //

  this->DetachObservation( observation );

  unsigned long tag;

  tag = observation->GetSubject()->AddObserver( vtkCommand::DeleteEvent, observation->GetObservationCallbackCommand() );
  observation->SetSubjectDeleteEventTag( tag );

  tag = observation->GetObserver()->AddObserver( vtkCommand::DeleteEvent, observation->GetObservationCallbackCommand() );
  observation->SetObserverDeleteEventTag( tag );

  tag = observation->GetSubject()->AddObserver( observation->GetEvent(), observation->GetObservationCallbackCommand() );
  observation->SetEventTag( tag );

}

//----------------------------------------------------------------------------
void vtkEventBroker::DetachObservation ( vtkObservation *observation )
{
  //
  // if the observation is attached, detach it
  //

  if ( observation->GetSubjectDeleteEventTag() )
    {
    observation->GetSubject()->RemoveObserver( observation->GetSubjectDeleteEventTag() );
    observation->SetSubjectDeleteEventTag( 0 );
    }

  if ( observation->GetObserverDeleteEventTag() )
    {
    observation->GetObserver()->RemoveObserver( observation->GetObserverDeleteEventTag() );
    observation->SetObserverDeleteEventTag( 0 );
    }

  if ( observation->GetEventTag() )
    {
    observation->GetSubject()->RemoveObserver( observation->GetEventTag() );
    observation->SetEventTag( 0 );
    }
}

//----------------------------------------------------------------------------
void vtkEventBroker::RemoveObservation ( vtkObservation *observation )
{
  std::vector< vtkObservation *> removeList;
  removeList.push_back( observation );
  this->RemoveObservations( removeList );
}

//----------------------------------------------------------------------------
void vtkEventBroker::RemoveObservations (std::vector< vtkObservation *>observations)
{
  // remove passed observations from:
  // - brokers observation list
  // - current event queue
  // - detach from subject (and observer)
  // - delete the observation

  // make a new broker observation list that doesn't include the passed observations
  std::vector< vtkObservation *> newObservations;
  std::vector< vtkObservation *>::iterator obsIter; 
  for(obsIter=this->Observations.begin(); obsIter != this->Observations.end(); obsIter++)  
    { 
    // foreach of the broker's observations see if it is in the list of items to be removed
    std::vector< vtkObservation *>::iterator searchIter;
    bool inRemoveList = false;
    for(searchIter=observations.begin(); searchIter != observations.end(); searchIter++)  
      {
      if (*obsIter == *searchIter)
        {
        inRemoveList = true;
        break;
        }
      }
    if ( !inRemoveList )
      {
      newObservations.push_back( *obsIter );
      }
    }
  this->Observations = newObservations;

  // remove from event queue
  std::deque< vtkObservation *> newEventQueue;
  std::deque< vtkObservation *>::iterator queueIter; 
  for(queueIter=this->EventQueue.begin(); queueIter != this->EventQueue.end(); queueIter++)  
    { 
    // foreach of the broker's observations see if it is in the list of items to be removed
    std::vector< vtkObservation *>::iterator searchIter;
    bool inRemoveList = false;
    for(searchIter=observations.begin(); searchIter != observations.end(); searchIter++)  
      {
      if (*queueIter == *searchIter)
        {
        inRemoveList = true;
        break;
        }
      }
    if ( !inRemoveList )
      {
      newEventQueue.push_back( *queueIter );
      }
    }
  this->EventQueue = newEventQueue;

  // detach and delete each of the observations
  std::vector< vtkObservation *>::iterator removeIter;
  for(removeIter=observations.begin(); removeIter != observations.end(); removeIter++)  
    {
    this->DetachObservation( *removeIter );
    (*removeIter)->Delete();
    }
}

//----------------------------------------------------------------------------
void vtkEventBroker::RemoveObservations (vtkObject *observer)
{
  // find matching observations to remove
  std::vector< vtkObservation *> removeList;
  std::vector< vtkObservation *>::iterator obsIter; 
  for(obsIter=this->Observations.begin(); obsIter != this->Observations.end(); obsIter++)  
    {
    if ( (*obsIter)->GetObserver() == observer )
      {
      removeList.push_back( *obsIter );
      }
    }
  this->RemoveObservations( removeList );
}

//----------------------------------------------------------------------------
void vtkEventBroker::RemoveObservations (vtkObject *subject, vtkObject *observer)
{
  // find matching observations to remove
  std::vector< vtkObservation *> removeList;
  std::vector< vtkObservation *>::iterator obsIter; 
  for(obsIter=this->Observations.begin(); obsIter != this->Observations.end(); obsIter++)  
    {
    if ( (*obsIter)->GetObserver() == observer && 
         (*obsIter)->GetSubject() == subject )
      {
      removeList.push_back( *obsIter );
      }
    }
  this->RemoveObservations( removeList );
}

//----------------------------------------------------------------------------
void vtkEventBroker::RemoveObservations (vtkObject *subject, unsigned long event, vtkObject *observer)
{
  // find matching observations to remove
  std::vector< vtkObservation *> removeList;
  std::vector< vtkObservation *>::iterator obsIter; 
  for(obsIter=this->Observations.begin(); obsIter != this->Observations.end(); obsIter++)  
    {
    if ( (*obsIter)->GetObserver() == observer && 
         (*obsIter)->GetSubject() == subject &&
         (*obsIter)->GetEvent() == event )
      {
      removeList.push_back( *obsIter );
      }
    }
  this->RemoveObservations( removeList );
}

//----------------------------------------------------------------------------
void vtkEventBroker::RemoveObservations (vtkObject *subject, unsigned long event, vtkObject *observer, vtkCallbackCommand *notify)
{
  // find matching observations to remove
  std::vector< vtkObservation *> removeList;
  std::vector< vtkObservation *>::iterator obsIter; 
  for(obsIter=this->Observations.begin(); obsIter != this->Observations.end(); obsIter++)  
    {
    if ( (*obsIter)->GetObserver() == observer && 
         (*obsIter)->GetSubject() == subject &&
         (*obsIter)->GetEvent() == event && 
         (*obsIter)->GetCallbackCommand() == notify )
      {
      removeList.push_back( *obsIter );
      }
    }
  this->RemoveObservations( removeList );
}
//----------------------------------------------------------------------------
vtkCollection *vtkEventBroker::GetObservationsForSubject ( vtkObject *subject )
{
  vtkCollection *collection = vtkCollection::New();
  std::vector< vtkObservation *>::iterator iter; 
  for(iter=this->Observations.begin(); iter != this->Observations.end(); iter++)  
    { 
    if ( (*iter)->GetSubject() == subject )
      {
      collection->AddItem( *iter );
      }
    }
  return collection;
}

//----------------------------------------------------------------------------
vtkCollection *vtkEventBroker::GetObservationsForObserver ( vtkObject *observer )
{
  vtkCollection *collection = vtkCollection::New();
  std::vector< vtkObservation *>::iterator iter; 
  for(iter=this->Observations.begin(); iter != this->Observations.end(); iter++)  
    { 
    if ( (*iter)->GetObserver() == observer )
      {
      collection->AddItem( *iter );
      }
    }
  return collection;
}

//----------------------------------------------------------------------------
int vtkEventBroker::GetNumberOfObservations ( )
{
  return (this->Observations.size());
}

//----------------------------------------------------------------------------
vtkObservation *vtkEventBroker::GetNthObservation ( int n )
{
  if ( n < 0 || n >= this->GetNumberOfObservations() )
    {
    return NULL;
    }
  return (this->Observations[n]);
}

//----------------------------------------------------------------------------
int vtkEventBroker::GenerateGraphFile ( const char *graphFile )
{
  std::ofstream file;

  file.open( graphFile, std::ios::out );

  if ( file.fail() )
    {
    vtkErrorMacro( "could not write to " << graphFile );
    return 1;
    }

  file << "digraph G {\n";

  vtkObservation *observation;
  int size = this->GetNumberOfObservations();
  for (int count = 0; count < size; count++)
    {
    observation = this->GetNthObservation( count );
    file << " " \
        << observation->GetObserver()->GetClassName() 
        << " -> "
        << observation->GetSubject()->GetClassName() 
        << " [ label = \"" 
        << vtkCommand::GetStringFromEventId( observation->GetEvent() )
        << "\" ]\n;" ;
    }


  file << "}\n";
  file.close();
  return 0;
}

//----------------------------------------------------------------------------
void vtkEventBroker::LogEvent ( vtkObservation *observation )
{
  if ( this->LogFileName == NULL )
    {
    // if we don't have a log file, we can't do anything
    return;
    }

  // start an event log if needed
  if ( this->EventLogging && !this->LogFile.is_open() )
    {
    this->LogFile.open( this->LogFileName, std::ios::out );
    this->LogFile << "strict digraph G {\n";
    }

  // close the log if done
  if ( !this->EventLogging && this->LogFile.is_open() )
    {
    this->LogFile << "}\n";
    this->LogFile.close();
    return;
    }

  if ( this->EventLogging && observation != NULL )
    {
    // log the actual event
    this->LogFile << " " \
        << observation->GetSubject()->GetClassName() 
        << " -> "
        << observation->GetObserver()->GetClassName() 
        << " [ label = \"" 
        << vtkCommand::GetStringFromEventId( observation->GetEvent() )
        << "\" ]\n;" ;
    }
}

//----------------------------------------------------------------------------
void vtkEventBroker::ProcessEvent ( vtkObservation *observation, vtkObject *caller, unsigned long eid, void *callData )
{
  //
  // we've got an event we care about - either invoke it
  // right away, or put it on the queue for later handling
  // - check first if it's the event the observer asked for (it might
  // be a delete event that the event broker asked for)
  //
  if ( eid == observation->GetEvent() )
    {
    if ( this->EventMode == vtkEventBroker::Synchronous )
      {
      this->InvokeObservation( observation );
      }
    else if ( this->EventMode == vtkEventBroker::Asynchronous )
      {
      this->QueueObservation( observation );
      }
    else
      {
      vtkErrorMacro ( "Bad EventMode " << this->EventMode );
      }
    }

  // 
  // for delete events, just clean up and get out
  // - we want to forget we ever had this observation
  //
  if ( eid == vtkCommand::DeleteEvent )
    {
    this->RemoveObservation( observation );
    return;
    }
}


//----------------------------------------------------------------------------
void vtkEventBroker::QueueObservation ( vtkObservation *observation )
{
  if ( observation->GetInEventQueue())
    {
    return;
    }
  this->EventQueue.push_back( observation );
  observation->SetInEventQueue(1);
}

//----------------------------------------------------------------------------
int vtkEventBroker::GetNumberOfQueuedObservations ()
{
  return( this->EventQueue.size() );
}

//----------------------------------------------------------------------------
vtkObservation *vtkEventBroker::GetNthQueuedObservation ( int n )
{
  if ( n < 0 || n >= this->GetNumberOfQueuedObservations() )
    {
    return NULL;
    }
  return (this->EventQueue[n]);
}

//----------------------------------------------------------------------------
vtkObservation *vtkEventBroker::DequeueObservation ()
{
  vtkObservation *observation = this->EventQueue.front();
  this->EventQueue.pop_front();
  observation->SetInEventQueue(0);
  return( observation );
}

//----------------------------------------------------------------------------
void vtkEventBroker::InvokeObservation ( vtkObservation *observation )
{
  // TODO record the timing before and after invocation
 
  // Invoke the observation
  observation->GetCallbackCommand()->Execute(
                                    observation->GetSubject(),
                                    observation->GetEvent(),
                                    NULL);

  // Write the to the log file if enabled
  this->LogEvent (observation);

}

//----------------------------------------------------------------------------
void vtkEventBroker::ProcessEventQueue ()
{
  while ( this->GetNumberOfQueuedObservations() > 0 )
    {
    this->InvokeObservation( this->DequeueObservation() );
    }
}

//----------------------------------------------------------------------------
void vtkEventBroker::PrintSelf(ostream& os, vtkIndent indent)
{
  this->vtkObject::PrintSelf(os, indent);
  
  os << indent << "NumberOfObservations: " << this->GetNumberOfObservations() << "\n";
  os << indent << "NumberOfQueueObservations: " << this->GetNumberOfQueuedObservations() << "\n";
  os << indent << "EventMode: " << this->GetEventModeAsString() << "\n";
  os << indent << "EventLogging: " << this->EventLogging << "\n";
  os << indent << "LogFileName: " <<
    (this->LogFileName ? this->LogFileName : "(none)") << "\n";
}

//----------------------------------------------------------------------------
// Description:
// the Callback is a static function to relay events 
//
void 
vtkEventBroker::Callback(vtkObject *caller, 
            unsigned long eid, void *clientData, void *callData)
{
  vtkObservation *observation = reinterpret_cast<vtkObservation *>(clientData);
  vtkEventBroker *self = observation->GetEventBroker();

  vtkDebugWithObjectMacro(self, "In vtkEvenBroker Callback");

  self->ProcessEvent(observation, caller, eid, callData);
}

