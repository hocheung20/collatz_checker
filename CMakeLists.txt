cmake_minimum_required(VERSION 3.31)
project(collatz_checker CXX)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_C_STANDARD 11)

find_package(GMP REQUIRED)
find_package(Boost CONFIG REQUIRED system iostreams)
find_package(ParallelExecutorWithMaxPendingTasks REQUIRED)
find_package(Threads REQUIRED)

ADD_EXECUTABLE(collatz_checker main.cpp)

TARGET_INCLUDE_DIRECTORIES(collatz_checker PRIVATE ${GMP_INCLUDES} ${Boost_INCLUDE_DIR} ${ParallelExecutorWithMaxPendingTasks_INCLUDES})
TARGET_LINK_LIBRARIES(collatz_checker ${GMP_LIBRARIES} ${Boost_LIBRARIES} Threads::Threads)
