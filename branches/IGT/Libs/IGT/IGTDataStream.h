
#ifndef __IGTDATASTREAM_H
#define __IGTDATASTREAM_H


#include "vtkMatrix4x4.h"


class vtkIGTDataStream : public vtkMatrix4x4
{
public:

  // Constructors/Destructors
  //  


  /**
   * Empty Constructor
   */
  vtkIGTDataStream ( );

  /**
   * Empty Destructor
   */
  virtual ~vtkIGTDataStream ( );

  // Static Public attributes
  //  

  // Public attributes
  //  


  // Public attribute accessor methods
  //  


  // Public attribute accessor methods
  //  


protected:

  /**
   */
  void create_MRML_node ( );


  /**
   */
  void update_mrml ( );



private:

  // Static Private attributes
  //  

  // Private attributes
  //  

  int m_buffer_size;
  int m_LastInputNum;
  int m_LastInputTime;

  // Private attribute accessor methods
  //  


  // Private attribute accessor methods
  //  


  /**
   * Set the value of m_buffer_size
   * @param new_var the new value of m_buffer_size
   */
  void setBuffer_size ( int new_var );

  /**
   * Get the value of m_buffer_size
   * @return the value of m_buffer_size
   */
  int getBuffer_size ( );


  /**
   * Set the value of m_LastInputNum
   * @param new_var the new value of m_LastInputNum
   */
  void setLastInputNum ( int new_var );

  /**
   * Get the value of m_LastInputNum
   * @return the value of m_LastInputNum
   */
  int getLastInputNum ( );


  /**
   * Set the value of m_LastInputTime
   * @param new_var the new value of m_LastInputTime
   */
  void setLastInputTime ( int new_var );

  /**
   * Get the value of m_LastInputTime
   * @return the value of m_LastInputTime
   */
  int getLastInputTime ( );


  void initAttributes ( ) ;

};

#endif // __IGTDATASTREAM_H
