cmake_minimum_required(VERSION 3.0.0)
project(Citrine5)

add_definitions(-Ofast)

# set(Eigen3_DIR "C:/Program Files/Eigen3/share/eigen3/cmake")
find_package(eigen3)

include_directories(Head)
include_directories(${EIGEN3_INCLUDE_DIR})
message(${EIGEN3_INCLUDE_DIR})

aux_source_directory(Cpp srcs)
add_executable(${PROJECT_NAME} ${srcs})
target_compile_options(Citrine5 PRIVATE -fno-var-tracking-assignments)