# Install script for directory: /home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libImath.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libImath.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libImath.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/jacobacfr/src/beifong/build/ext_build/openexr/IlmBase/Imath/libImath.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libImath.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libImath.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libImath.so"
         OLD_RPATH "/home/jacobacfr/src/beifong/build/ext_build/openexr/IlmBase/Iex:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libImath.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/OpenEXR" TYPE FILE FILES
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathBoxAlgo.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathBox.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathColorAlgo.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathColor.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathEuler.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathExc.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathExport.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathForward.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathFrame.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathFrustum.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathFrustumTest.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathFun.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathGL.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathGLU.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathHalfLimits.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathInt64.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathInterval.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathLimits.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathLineAlgo.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathLine.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathMath.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathMatrixAlgo.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathMatrix.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathNamespace.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathPlane.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathPlatform.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathQuat.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathRandom.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathRoots.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathShear.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathSphere.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathVecAlgo.h"
    "/home/jacobacfr/src/beifong/ext/openexr/IlmBase/Imath/ImathVec.h"
    )
endif()

