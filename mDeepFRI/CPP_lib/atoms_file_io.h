#ifndef ATOMS_FILE_IO
#define ATOMS_FILE_IO

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <fstream>

namespace py = boost::python;
namespace np = py::numpy;

static void SaveAtomsFile(const np::ndarray &position_array,
                          const np::ndarray &groups_array,
                          const std::string &save_path);

static std::tuple<int, int *, float *>
LoadAtomsFile(const std::string &file_path);

#endif
