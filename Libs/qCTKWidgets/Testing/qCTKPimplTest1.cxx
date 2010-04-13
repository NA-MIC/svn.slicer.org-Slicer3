/*=========================================================================

  Library:   qCTK

  Copyright (c) Kitware Inc. 
  All rights reserved.
  Distributed under a BSD License. See LICENSE.txt file.

  This software is distributed "AS IS" WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the above copyright notice for more information.

=========================================================================*/

// qCTK includes
#include "qCTKPimpl.h"

// QT includes
#include <QApplication>

// STD includes
#include <stdlib.h>
#include <iostream>

class qCTKPimplHelper 
{
public:

};

int qCTKPimplTest1(int argc, char * argv [] )
{
  QApplication app(argc, argv);


  return EXIT_SUCCESS;
}

