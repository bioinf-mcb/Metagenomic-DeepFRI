//
// Created by soliareofastora on 27.05.2021.
//

#include <boost/python.hpp>

#include "atoms_file_io.h"
#include "load_contact_maps.h"
#include "python_utils.h"

namespace py = boost::python;
namespace np = py::numpy;

BOOST_PYTHON_MODULE (libAtomDistanceIO) {
  py::def("initialize", Initialize);
  py::def("save_atoms", SaveAtomsFile);
  py::def("load_contact_map", LoadContactMap);
  py::def("load_aligned_contact_map", LoadAlignedContactMap);
}
