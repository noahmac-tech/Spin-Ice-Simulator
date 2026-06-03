cmake_minimum_required(VERSION 3.10)

project(SpinIceSimulator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

include_directories(src/include)

add_executable(spin_ice src/spin_ice.cpp)
target_link_libraries(spin_ice PRIVATE Threads::Threads)

add_executable(spin_ice_parallel src/spin_ice_parallel.cpp)
target_link_libraries(spin_ice_parallel PRIVATE Threads::Threads)

add_executable(low_temp_entropy src/low_temp_entropy.cpp)
target_link_libraries(low_temp_entropy PRIVATE Threads::Threads)

add_executable(low_temp_entropy_parallel src/low_temp_entropy_parallel.cpp)
target_link_libraries(low_temp_entropy_parallel PRIVATE Threads::Threads)