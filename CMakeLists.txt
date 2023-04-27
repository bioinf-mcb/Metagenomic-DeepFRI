cmake_minimum_required(VERSION 3.19)
project(GenomicLibrary)

set(CMAKE_CXX_STANDARD 17)

add_library(AtomDistanceIO SHARED
        meta_deepFRI/CPP_lib/library_definition.cpp
        meta_deepFRI/CPP_lib/atoms_file_io.h
        meta_deepFRI/CPP_lib/python_utils.h
        meta_deepFRI/CPP_lib/load_contact_maps.h)


target_include_directories(AtomDistanceIO PUBLIC ${PY_INCLUDE_PATH})

FIND_PACKAGE( Boost COMPONENTS python numpy REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

TARGET_LINK_LIBRARIES( AtomDistanceIO LINK_PUBLIC ${Boost_LIBRARIES} )
