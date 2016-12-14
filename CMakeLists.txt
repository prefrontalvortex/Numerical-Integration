cmake_minimum_required(VERSION 3.6)
project(integrator)

set(CMAKE_CXX_STANDARD 11)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -lm -static-libgcc -lquadmath -lgfortran -std=c++11 ")


set(SOURCE_FILES main.c argparse.c argparse.h aux.c aux.h rng.c rng.h main.h functions.c functions.h)
add_executable(integrator ${SOURCE_FILES})

set_target_properties( integrator PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin/ )

target_link_libraries(integrator m Threads::Threads)