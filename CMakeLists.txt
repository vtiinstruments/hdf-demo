cmake_minimum_required(VERSION 3.15)
project(hdf_demo VERSION 1.1 LANGUAGES C CXX)

find_package(HDF5 REQUIRED COMPONENTS C CXX)

include_directories(${HDF5_INCLUDE_DIRS})

add_executable(hdf_demo hdf_demo.cpp hdf_inspection.cpp hdf_post_process.cpp hdf_demo.h hdf_inspection.h hdf_post_process.h hdf_utils.h)
target_link_libraries(hdf_demo PRIVATE ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
