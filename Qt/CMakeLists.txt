
# --------------------------------------------------------------------------
set(MODULE_NAME "EMSegment")
set(MODULE_TITLE "${MODULE_NAME}")

string(TOUPPER ${MODULE_NAME} MODULE_NAME_UPPER)

#-----------------------------------------------------------------------------
# Add subdirectories

# Logic
ADD_SUBDIRECTORY(Logic)

# Qt Widgets
ADD_SUBDIRECTORY(Widgets)

# The Python-based EMSegment wizard
ADD_SUBDIRECTORY(Py)

#-----------------------------------------------------------------------------
set(MODULE_EXPORT_DIRECTIVE "Q_SLICER_QTMODULES_${MODULE_NAME_UPPER}_EXPORT")

set(MODULE_INCLUDE_DIRECTORIES
  ${CMAKE_CURRENT_SOURCE_DIR}/Logic
  ${CMAKE_CURRENT_BINARY_DIR}/Logic
  )

set(MODULE_SRCS
  qSlicer${MODULE_NAME}Module.cxx
  qSlicer${MODULE_NAME}Module.h
  )

set(MODULE_MOC_SRCS
  qSlicer${MODULE_NAME}Module.h
  )

set(MODULE_UI_SRCS
  )

if(${VTK_VERSION_MAJOR} GREATER 5)
  set(VTK_LIBRARIES
    vtkChartsCore
    )
else()
  set(VTK_LIBRARIES
    vtkCharts
    )
endif()

set(MODULE_TARGET_LIBRARIES
  vtkSlicer${MODULE_NAME}ModuleLogic
  qSlicer${MODULE_NAME}ModuleWidgets
  MRMLLogic
  ${VTK_LIBRARIES}
  )

set(MODULE_RESOURCES
  Resources/qSlicer${MODULE_NAME}Module.qrc
  )

#-----------------------------------------------------------------------------
slicerMacroBuildLoadableModule(
  NAME ${MODULE_NAME}
  TITLE ${MODULE_TITLE}
  EXPORT_DIRECTIVE ${MODULE_EXPORT_DIRECTIVE}
  INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
  SRCS ${MODULE_SRCS}
  MOC_SRCS ${MODULE_MOC_SRCS}
  UI_SRCS ${MODULE_UI_SRCS}
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  RESOURCES ${MODULE_RESOURCES}
  )

if(BUILD_TESTING)
  SET(EMSegment_TESTING_DATA_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Testing/Data)

  #ADD_SUBDIRECTORY(Testing)
endif()

# --------------------------------------------------------------------------
set(MODULE_NAME "EMSegmentQuick")

string(TOUPPER ${MODULE_NAME} MODULE_NAME_UPPER)

set(MODULE_EXPORT_DIRECTIVE "Q_SLICER_QTMODULES_${MODULE_NAME_UPPER}_EXPORT")

set(MODULE_SRCS
  qSlicer${MODULE_NAME}Module.cxx
  qSlicer${MODULE_NAME}Module.h
  )

set(MODULE_MOC_SRCS
  qSlicer${MODULE_NAME}Module.h
  )

set(MODULE_UI_SRCS
  )

set(MODULE_RESOURCES
  )

#-----------------------------------------------------------------------------
slicerMacroBuildLoadableModule(
  NAME ${MODULE_NAME}
  EXPORT_DIRECTIVE ${MODULE_EXPORT_DIRECTIVE}
  INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
  SRCS ${MODULE_SRCS}
  MOC_SRCS ${MODULE_MOC_SRCS}
  UI_SRCS ${MODULE_UI_SRCS}
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  RESOURCES ${MODULE_RESOURCES}
  NO_TITLE # See note below
  )

# Note: Usually, one module per directory is defined. Trying to build two modules from the same
#       directory will lead to a multiple definition warning related to QTMODULE_TITLE macro.
#       By passing the option NO_TITLE, the macro won't be re-define and we will avoid the issue.
#       Doing so implies that we directly hardcode the title in qSlicerEMSegmentQuickModule.h


set(MODULE_NAME "EMSegment") # Reset module name so that logic/MRML/... build
                             # macros associate the expected folder with
                             # remaining module targets.
