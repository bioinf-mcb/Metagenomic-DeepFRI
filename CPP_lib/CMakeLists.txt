cmake_minimum_required(VERSION 3.19)
project(GenomicLibrary)

set(CMAKE_CXX_STANDARD 17)

add_library(AtomDistanceIO SHARED
        library.cpp
        AtomFileIO.cpp PythonUtils.cpp)


target_include_directories(AtomDistanceIO PUBLIC ~/miniconda3/include/python3.8)

FIND_PACKAGE( Boost COMPONENTS python numpy REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

TARGET_LINK_LIBRARIES( AtomDistanceIO LINK_PUBLIC ${Boost_LIBRARIES} )

add_custom_command(TARGET AtomDistanceIO POST_BUILD
        COMMAND "${CMAKE_COMMAND}" -E copy
        "libAtomDistanceIO.so"
        "../libAtomDistanceIO.so"
        COMMENT "Copying to output directory")