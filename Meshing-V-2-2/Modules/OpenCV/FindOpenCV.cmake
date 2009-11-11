MESSAGE(STATUS "Looking for OpenCV")

MACRO(DBGMSG MSG)
    MESSAGE(STATUS ${MSG})
ENDMACRO(DBGMSG)

#
# Windows registry key for installation root.
# インストールルートパス取得用 Windows レジストリキー
#
GET_FILENAME_COMPONENT(OpenCV_ROOT_DIR_CANDIDATES
    "[HKEY_CURRENT_USER\\Software\\Intel\\OpenCV\\Settings;Path]"
    ABSOLUTE CACHE
)

#
# Installation root candidate list.
# インストールルートパス調査対象
#

SET(OpenCV_ROOT_DIR_CANDIDATES ${OpenCV_ROOT_DIR_CANDIDATES}
    ${OpenCV_DIR}
    $ENV{ProrgramFiles}/OpenCV
)

SET(OpenCV_INCLUDE_DIRS ${OpenCV_DIR}/include/opencv)
SET(OpenCV_BIN_DIR ${OpenCV_DIR}/bin)
DBGMSG("OpenCV_INC_DIR: ${OpenCV_INC_DIR}")
DBGMSG("OpenCV_BIN_DIR: ${OpenCV_BIN_DIR}")

IF(OpenCV_BIN_DIR)
  SET(OpenCV_FOUND ON)
ELSE(OpenCV_BIN_DIR)
  SET(OpenCV_FOUND 0)
  MESSAGE(FATAL_ERROR "Failed to find OpenCV lib.")
ENDIF(OpenCV_BIN_DIR)

#
# Sub directory candidate list for including.
# インクルードディレクトリ調査対象
#
SET(OpenCV_INCDIR_SUFFIXES
    include
    include/cv
    include/opencv
    cv/include
    cxcore/include
    cvaux/include
    otherlibs/cvcam/include
    otherlibs/highgui
    otherlibs/highgui/include
    otherlibs/_graphics/include
)

#
# Find sub directories for including.
# インクルードディレクトリ探索
#
FIND_PATH(OpenCV_CV_INCLUDE_DIR
    NAMES cv.h
    PATHS ${OpenCV_DIR}
    PATH_SUFFIXES ${OpenCV_INCDIR_SUFFIXES}
)
DBGMSG("OpenCV_CV_INCLUDE_DIR: ${OpenCV_CV_INCLUDE_DIR}")

FIND_PATH(OpenCV_CXCORE_INCLUDE_DIR
    NAMES cxcore.h
    PATHS ${OpenCV_DIR}
    PATH_SUFFIXES ${OpenCV_INCDIR_SUFFIXES}
)
DBGMSG("OpenCV_CXCORE_INCLUDE_DIR: ${OpenCV_CXCORE_INCLUDE_DIR}")

FIND_PATH(OpenCV_CVAUX_INCLUDE_DIR
    NAMES cvaux.h
    PATHS ${OpenCV_DIR}
    PATH_SUFFIXES ${OpenCV_INCDIR_SUFFIXES}
)
DBGMSG("OpenCV_CVAUX_INCLUDE_DIR: ${OpenCV_CVAUX_INCLUDE_DIR}")

FIND_PATH(OpenCV_HIGHGUI_INCLUDE_DIR
    NAMES highgui.h
    PATHS ${OpenCV_DIR}
    PATH_SUFFIXES ${OpenCV_INCDIR_SUFFIXES}
)
DBGMSG("OpenCV_HIGHGUI_INCLUDE_DIR: ${OpenCV_HIGHGUI_INCLUDE_DIR}")

FIND_PATH(OpenCV_CVCAM_INCLUDE_DIR
    NAMES cvcam.h
    PATHS ${OpenCV_DIR}
    PATH_SUFFIXES ${OpenCV_INCDIR_SUFFIXES}
)
DBGMSG("OpenCV_CVCAM_INCLUDE_DIR: ${OpenCV_CVCAM_INCLUDE_DIR}")



#
# Sub directory candidate list for libraries
# ライブラリディレクトリ調査対象
#
SET(OpenCV_LIBDIR_SUFFIXES
    lib
    OpenCV/lib
    otherlibs/_graphics/lib
)

#
# Find sub directories for libraries.
# ライブラリディレクトリ探索
#
FIND_LIBRARY(OpenCV_CV_LIBRARY
    NAMES cv opencv
    PATHS ${OpenCV_DIR}
    PATH_SUFFIXES ${OpenCV_LIBDIR_SUFFIXES}
)
DBGMSG("OpenCV_CV_LIBRARY: ${OpenCV_CV_LIBRARY}")

FIND_LIBRARY(OpenCV_CVAUX_LIBRARY
    NAMES cvaux
    PATHS ${OpenCV_DIR}
    PATH_SUFFIXES ${OpenCV_LIBDIR_SUFFIXES}
)
DBGMSG("OpenCV_CVAUX_LIBRARY: ${OpenCV_CVAUX_LIBRARY}")

FIND_LIBRARY(OpenCV_CVCAM_LIBRARY
    NAMES cvcam
    PATHS ${OpenCV_DIR}
    PATH_SUFFIXES ${OpenCV_LIBDIR_SUFFIXES}
)
DBGMSG("OpenCV_CVCAM_LIBRARY: ${OpenCV_CVCAM_LIBRARY}")

FIND_LIBRARY(OpenCV_CVHAARTRAINING_LIBRARY
    NAMES cvhaartraining
    PATHS ${OpenCV_DIR}
    PATH_SUFFIXES ${OpenCV_LIBDIR_SUFFIXES}
)
DBGMSG("OpenCV_CVHAARTRAINING_LIBRARY: ${OpenCV_CVHAARTRAINING_LIBRARY}")

FIND_LIBRARY(OpenCV_CXCORE_LIBRARY
    NAMES cxcore
    PATHS ${OpenCV_DIR}
    PATH_SUFFIXES ${OpenCV_LIBDIR_SUFFIXES}
)
DBGMSG("OpenCV_CXCORE_LIBRARY: ${OpenCV_CXCORE_LIBRARY}")

FIND_LIBRARY(OpenCV_CXTS_LIBRARY
    NAMES cxts
    PATHS ${OpenCV_DIR}
    PATH_SUFFIXES ${OpenCV_LIBDIR_SUFFIXES}
)
DBGMSG("OpenCV_CXTS_LIBRARY: ${OpenCV_CXTS_LIBRARY}")

FIND_LIBRARY(OpenCV_HIGHGUI_LIBRARY
    NAMES highgui
    PATHS ${OpenCV_DIR}
    PATH_SUFFIXES ${OpenCV_LIBDIR_SUFFIXES}
)
DBGMSG("OpenCV_HIGHGUI_LIBRARY: ${OpenCV_HIGHGUI_LIBRARY}")

FIND_LIBRARY(OpenCV_ML_LIBRARY
    NAMES ml
    PATHS ${OpenCV_DIR}
    PATH_SUFFIXES ${OpenCV_LIBDIR_SUFFIXES}
)
DBGMSG("OpenCV_ML_LIBRARY: ${OpenCV_ML_LIBRARY}")

FIND_LIBRARY(OpenCV_TRS_LIBRARY
    NAMES trs
    PATHS ${OpenCV_DIR}
    PATH_SUFFIXES ${OpenCV_LIBDIR_SUFFIXES}
)
DBGMSG("OpenCV_TRS_LIBRARY: ${OpenCV_TRS_LIBRARY}")


#
# SET CONFIGURATION VARIABLES
# 変数設定
#
SET(OpenCV_REQUIRED_COMPONENTS CV CXCORE CVAUX HIGHGUI)
FOREACH (NAME ${OpenCV_REQUIRED_COMPONENTS})
    IF (OpenCV_${NAME}_INCLUDE_DIR AND OpenCV_${NAME}_LIBRARY)
        LIST(APPEND OpenCV_INCLUDE_DIRS ${OpenCV_${NAME}_INCLUDE_DIR})
        LIST(APPEND OpenCV_LIBRARIES ${OpenCV_${NAME}_LIBRARY})
    ELSE (OpenCV_${NAME}_INCLUDE_DIR AND OpenCV_${NAME}_LIBRARY)
        SET(OpenCV_FOUND OFF)
    ENDIF (OpenCV_${NAME}_INCLUDE_DIR AND OpenCV_${NAME}_LIBRARY)
ENDFOREACH (NAME)

IF(OpenCV_FOUND)
    MESSAGE(STATUS "Found OpenCV")
ENDIF(OpenCV_FOUND)

IF (OpenCV_CV_LIBRARY)
    GET_FILENAME_COMPONENT(OpenCV_LINK_DIRECTORIES ${OpenCV_CV_LIBRARY} PATH)
ENDIF (OpenCV_CV_LIBRARY)

MARK_AS_ADVANCED(
    OpenCV_ROOT_DIR
    OpenCV_INCLUDE_DIRS
    OpenCV_CV_INCLUDE_DIR
    OpenCV_CXCORE_INCLUDE_DIR
    OpenCV_CVAUX_INCLUDE_DIR
    OpenCV_CVCAM_INCLUDE_DIR
    OpenCV_HIGHGUI_INCLUDE_DIR
    OpenCV_LIBRARIES
    OpenCV_CV_LIBRARY
    OpenCV_CXCORE_LIBRARY
    OpenCV_CVAUX_LIBRARY
    OpenCV_CVCAM_LIBRARY
    OpenCV_CVHAARTRAINING_LIBRARY
    OpenCV_CXTS_LIBRARY
    OpenCV_HIGHGUI_LIBRARY
    OpenCV_ML_LIBRARY
    OpenCV_TRS_LIBRARY
)

