cmake_minimum_required(VERSION 3.19)
project(CPP_lib)

set(CMAKE_CXX_STANDARD 17)

add_library(CPP_lib SHARED
        library.cpp
        )


target_include_directories(CPP_lib PUBLIC ~/miniconda3/include/python3.8)

FIND_PACKAGE( Boost COMPONENTS python numpy REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

TARGET_LINK_LIBRARIES( CPP_lib LINK_PUBLIC ${Boost_LIBRARIES} )

add_custom_command(TARGET CPP_lib POST_BUILD
        COMMAND "${CMAKE_COMMAND}" -E copy
        "libCPP_lib.so"
        "../ContactMapper.so"
        COMMENT "Copying to output directory")