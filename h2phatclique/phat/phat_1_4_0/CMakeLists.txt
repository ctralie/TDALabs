cmake_minimum_required(VERSION 2.8)

project("Phat")

set(CMAKE_BUILD_TYPE "Release")

INCLUDE(CheckCXXSourceCompiles)

include_directories (include)

FIND_PACKAGE(OpenMP)

set (CMAKE_REQUIRED_FLAGS ${OpenMP_CXX_FLAGS})
CHECK_CXX_SOURCE_COMPILES("
#include <omp.h>
int main() {
#if (_OPENMP >= 200805 || _MSC_VER >= 1500)
  return 0;
#else
  breaks_on_purpose
#endif
}
" OPENMP_VERSION)

if(OPENMP_VERSION)



set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
message("

OpenMP 3.0 not supported by the compiler (${CMAKE_CXX_COMPILER})! 
To get optimal performance for the \"chunk\" algorithm, use a C++ compiler with OpenMP 3.0 support (e.g., GCC >=4.4).
To use a different compiler, pass it to cmake in the variable CMAKE_CXX_COMPILER: 
 cmake . -DCMAKE_CXX_COMPILER=g++-4.7

")
endif()
 
FILE(GLOB_RECURSE all_headers "include/phat/*.h")

FILE(GLOB general_includes "include/phat/*.h")
FILE(GLOB algorithms_includes "include/phat/algorithms/*.h")
FILE(GLOB helpers_includes "include/phat/helpers/*.h")
FILE(GLOB representations_includes "include/phat/representations/*.h")

add_executable (simple_example src/simple_example.cpp ${all_headers})
add_executable (self_test src/self_test.cpp ${all_headers})
add_executable (phat src/phat.cpp ${all_headers})
add_executable (info src/info.cpp ${all_headers})
add_executable (benchmark src/benchmark.cpp ${all_headers})
add_executable (convert src/convert.cpp ${all_headers})

source_group(Header\ Files FILES ${general_includes})
source_group(Header\ Files\\helpers FILES ${helpers_includes})
source_group(Header\ Files\\representations FILES ${representations_includes})
source_group(Header\ Files\\algorithms FILES ${algorithms_includes})
