cmake_minimum_required (VERSION 3.15)
project (StabilizedRK)

file(GLOB Library_Sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.cpp" ) 

add_library (StabilizedRK STATIC ${Library_Sources})

target_compile_features    (StabilizedRK PUBLIC cxx_std_11)
target_include_directories (StabilizedRK PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


