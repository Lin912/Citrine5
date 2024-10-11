cmake_minimum_required(VERSION 3.10)

# 设置项目名称
project(Citrine5)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/exe/)
#message(FATAL_ERROR  ${PROJECT_SOURCE_DIR})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-var-tracking-assignments")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-var-tracking-assignments")

# 设置编译器为g++
set(CMAKE_CXX_COMPILER g++)
# 设置编译器标准为C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 设置优化等级为ofast
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast")

# 添加头文件路径
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Head)

# 查找并包含Eigen库
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# 添加源文件
#file(GLOB SOURCE_FILES "Cpp/*.cpp")

set(SOURCE_FILES
    Cpp/main.cpp
    Cpp/Add.cpp
    Cpp/BC.cpp
    Cpp/Fiber.cpp
    Cpp/Fx.cpp
    Cpp/Iterator.cpp
    Cpp/Jacobian.cpp
    Cpp/Load.cpp
    Cpp/MNQ.cpp
    Cpp/ReadOut.cpp
)

# 添加可执行文件
add_executable(Citrine5 ${SOURCE_FILES})

# 链接Eigen库
target_link_libraries(Citrine5 Eigen3::Eigen)