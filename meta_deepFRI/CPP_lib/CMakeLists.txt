cmake_minimum_required(VERSION 3.19)
project(GenomicLibrary)

set(CMAKE_CXX_STANDARD 17)

add_library(AtomDistanceIO SHARED
        library_definition.cpp
        atoms_file_io.h
        python_utils.h
        load_contact_maps.h)


target_include_directories(AtomDistanceIO PUBLIC $ENV{CONDA_PREFIX}/include/python3.10)

FIND_PACKAGE( Boost COMPONENTS python numpy REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

TARGET_LINK_LIBRARIES( AtomDistanceIO LINK_PUBLIC ${Boost_LIBRARIES} )

add_custom_command(TARGET AtomDistanceIO POST_BUILD
        COMMAND "${CMAKE_COMMAND}" -E copy
        "libAtomDistanceIO.so"
        "../libAtomDistanceIO.so"
        COMMENT "Copying to output directory")
