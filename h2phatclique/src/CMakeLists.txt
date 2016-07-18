#set(CMAKE_BUILD_TYPE "debug")
set(CMAKE_BUILD_TYPE "release")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fomit-frame-pointer -fexpensive-optimizations -O3 -std=c++11")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -std=c++11")

# maybe needed -- test on non-Linux
#set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    # Linux specific code
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # OSX specific code
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -headerpad_max_install_names")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

FIND_PACKAGE(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

message(STATUS "phatclique_SOURCE_DIR: ${phatclique_SOURCE_DIR}")

include_directories(${phatclique_SOURCE_DIR}/include ${phatclique_SOURCE_DIR}/phat/phat_1_4_0/include)

# maybe needed -- test on non-Linux
#SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
# maybe needed -- test on non-Linux
#SET(BUILD_SHARED_LIBRARIES OFF)

add_executable (phatclique phatclique.cpp bmatrix.cpp ck.cpp)
set_target_properties( phatclique PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin )
