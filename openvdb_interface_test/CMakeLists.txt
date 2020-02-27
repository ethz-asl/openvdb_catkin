cmake_minimum_required (VERSION 3.10)
project(openvdb_interface_test)

find_package(catkin_simple REQUIRED)
catkin_simple(ALL_DEPS_REQUIRED)

message(WARNING "OPENVDB_CATKIN_LIBRARIES: ${OPENVDB_CATKIN_INCLUDE}")

cs_add_library(${PROJECT_NAME}
                src/dependency-tracker.cc)

catkin_add_gtest(test_openvdb_interface test/test-openvdb-interface.cc)
target_link_libraries(test_openvdb_interface ${PROJECT_NAME})

cs_install()
cs_export()
