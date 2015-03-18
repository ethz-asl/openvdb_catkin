cmake_minimum_required(VERSION 2.8.3)
project(openvdb_catkin)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_PREFIX}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})

find_package(catkin_simple REQUIRED)
find_package(PythonLibs REQUIRED)
find_package(Boost REQUIRED)
find_package(OpenEXR REQUIRED)
find_package(TBB REQUIRED)
find_package(CppUnit REQUIRED)
find_package(PythonInterp REQUIRED)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${BOOST_INCLUDE_DIR}
  ${PYTHON_INCLUDE_DIR}
  ${OPENEXR_INCLUDE_DIR}
  ${TBB_INCLUDE_DIR}
  ${CPPUNIT_INCLUDE_DIR}
)

catkin_simple()

include(ExternalProject)

file(MAKE_DIRECTORY ${CATKIN_DEVEL_PREFIX}/include/openvdb)

ExternalProject_Add(openvdb_src
  URL https://github.com/dreamworksanimation/openvdb/archive/v3.0.0.zip
  UPDATE_COMMAND ""
  PATCH_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND cd ../openvdb_src/openvdb && make -j8
  INSTALL_COMMAND cd ../openvdb_src/openvdb && make install -j8
)

cs_add_library(${PROJECT_NAME} src/dependency_tracker.cc)
add_dependencies(${PROJECT_NAME} openvdb_src)
target_link_libraries(${PROJECT_NAME} 
  ${CATKIN_DEVEL_PREFIX}/lib/libopenvdb${CMAKE_SHARED_LIBRARY_SUFFIX})

cs_install()
cs_export(
  INCLUDE_DIRS ${CATKIN_DEVEL_PREFIX}/include/openvdb/
  CFG_EXTRAS openvdb-extras.cmake)
