cmake_minimum_required(VERSION 3.10)
project(openvdb_catkin)

find_package(catkin_simple REQUIRED)
catkin_simple(ALL_DEPS_REQUIRED)

include(ExternalProject)

# BUILD c-blosh
# Recommended dependency for improved disk compression
SET(C_BLOSC_SRC_PATH "c-blosc_src")
SET(C_BLOSC_BUILD_PATH "c-blosc_build")
SET(C_BLOSC_DOWNLOAD_PATH "c-blosc_download")
ExternalProject_Add(c-blosc_src
  GIT_REPOSITORY git@github.com:Blosc/c-blosc.git
  GIT_TAG v1.5.0
  UPDATE_COMMAND ""
  DOWNLOAD_DIR ${C_BLOSC_DOWNLOAD_PATH}
  SOURCE_DIR ${C_BLOSC_SRC_PATH}
  BINARY_DIR ${C_BLOSC_BUILD_PATH}
  CMAKE_ARGS -DCMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}
             -DCMAKE_INSTALL_PREFIX=${CATKIN_DEVEL_PREFIX}
)

# BUILD OpenVDB
SET(OPENVDB_SRC_PATH "openvdb_src")
SET(OPENVDB_BUILD_PATH "openvdb_build")
SET(OPENVDB_DOWNLOAD_PATH "openvdb_download")
ExternalProject_Add(openvdb_src
  URL https://github.com/AcademySoftwareFoundation/openvdb/archive/v7.0.0.zip
  URL_MD5 048f1f9e77aa517bd778686c0638ffa4
  UPDATE_COMMAND ""
  DOWNLOAD_DIR ${OPENVDB_DOWNLOAD_PATH}
  SOURCE_DIR ${OPENVDB_SRC_PATH}
  BINARY_DIR ${OPENVDB_BUILD_PATH}
  CMAKE_ARGS -DCMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}
             -DCMAKE_INSTALL_PREFIX=${CATKIN_DEVEL_PREFIX}
             -DLOG4CPLUS_INCLUDES=""
             -DLOG4CPLUS_LIBRARIES=""
             -DBlosc_ROOT=${CATKIN_DEVEL_PREFIX}
             # Needed due for newer GCC versions,
             # source: https://github.com/AcademySoftwareFoundation/openvdb/issues/70
             -DCMAKE_NO_SYSTEM_FROM_IMPORTED:BOOL=TRUE
)
add_dependencies(openvdb_src c-blosc_src)

install(DIRECTORY ${CATKIN_DEVEL_PREFIX}/include/openvdb
        DESTINATION ${CATKIN_GLOBAL_INCLUDE_DESTINATION}
        FILES_MATCHING PATTERN "*.h"
                       PATTERN "*.hpp")
install(DIRECTORY ${CATKIN_DEVEL_PREFIX}/lib/
        DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
      FILES_MATCHING PATTERN "libopenvdb*")

cs_export(INCLUDE_DIRS ${CATKIN_DEVEL_PREFIX}/include
          LIBRARIES openvdb
                    Half
                    IlmImf
                    tbb
)
