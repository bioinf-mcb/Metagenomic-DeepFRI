cmake_minimum_required(VERSION 3.19)
project(ContactMapper)

set(CMAKE_CXX_STANDARD 17)

add_library(ContactMapper SHARED
        library.cpp
        )


target_include_directories(ContactMapper PUBLIC ~/miniconda3/include/python3.8)

FIND_PACKAGE( Boost COMPONENTS python numpy REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

TARGET_LINK_LIBRARIES( ContactMapper LINK_PUBLIC ${Boost_LIBRARIES} )

add_custom_command(TARGET ContactMapper POST_BUILD
        COMMAND "${CMAKE_COMMAND}" -E copy
        "libContactMapper.so"
        "../../libContactMapper.so"
        COMMENT "Copying to output directory")