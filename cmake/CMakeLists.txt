cmake_minimum_required(VERSION 2.8.3)
project(openvdb)

find_package(Boost REQUIRED)
find_package(CppUnit REQUIRED)
find_package(Half REQUIRED)
find_package(PythonLibs REQUIRED)
find_package(PythonInterp REQUIRED)
find_package(TBB REQUIRED)
find_package(OpenEXR REQUIRED)
find_package(Log4cplus REQUIRED)
find_package(Zlib REQUIRED)

include_directories(include)
include_directories(${OPENEXR_INCLUDE_DIRS})
include_directories(${HALF_INCLUDE_DIRS})
include_directories(${CPPUNIT_INCLUDE_DIRS})
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${BOOST_INCLUDE_DIRS})
include_directories(${TBB_INCLUDE_DIRS})
include_directories(${LOG4CPLUS_INCLUDES})
include_directories(${ZLIB_INCLUDE_DIRS})

include_directories(${PROJECT_SOURCE_DIR} )

set(BOOST_LIBRARIES boost_iostreams boost_system)

set(SRC_NAMES 
    openvdb/Grid.cc 
    openvdb/io/Archive.cc 
    openvdb/io/Compression.cc 
    openvdb/io/File.cc 
    openvdb/io/GridDescriptor.cc 
    openvdb/io/Queue.cc 
    openvdb/io/Stream.cc 
    openvdb/io/TempFile.cc 
    openvdb/math/Hermite.cc 
    openvdb/math/Maps.cc 
    openvdb/math/Proximity.cc 
    openvdb/math/QuantizedUnitVec.cc 
    openvdb/math/Transform.cc 
    openvdb/metadata/Metadata.cc 
    openvdb/metadata/MetaMap.cc 
    openvdb/openvdb.cc 
    openvdb/Platform.cc 
    openvdb/util/Formats.cc 
    openvdb/util/Util.cc)

add_library(${PROJECT_NAME} SHARED ${SRC_NAMES})
target_link_libraries(${PROJECT_NAME} ${TBB_LIBRARIES} ${BOOST_LIBRARIES} ${HALF_LIBRARIES} ${ZLIB_LIBRARIES})

install(TARGETS ${PROJECT_NAME}
	    RUNTIME DESTINATION bin
	    LIBRARY DESTINATION lib)
install(DIRECTORY openvdb/io openvdb/math openvdb/metadata openvdb/tools openvdb/tree openvdb/util DESTINATION include FILES_MATCHING PATTERN "*.h")
install(DIRECTORY openvdb/ DESTINATION include/openvdb
        FILES_MATCHING PATTERN "*.h")