//
// Created by soliareofastora on 24.12.2021.
//

#ifndef PYTHON_UTILS
#define PYTHON_UTILS

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace py = boost::python;
namespace np = py::numpy;

// Call this function before using other functions!!
static void Initialize() {
  Py_Initialize();
  boost::python::numpy::initialize();
}

static void DestroyCapsule(PyObject* self) {
  auto* b = reinterpret_cast<bool*>(PyCapsule_GetPointer(self, nullptr));
  delete[] b;
}

static np::ndarray CreateNumpyArray(bool* array, int size){
  // handling over memory management to python
  // https://stackoverflow.com/questions/57068443/setting-owner-in-boostpythonndarray-so-that-data-is-owned-and-managed-by-pyt
  PyObject* capsule = PyCapsule_New((void*) array, nullptr, (PyCapsule_Destructor) &DestroyCapsule);
  py::handle<> capsule_handle{capsule};
  py::object capsule_owner{capsule_handle};

  return np::from_data(array,
                       np::dtype::get_builtin<bool>(),
                       py::make_tuple(size, size),
                       py::make_tuple(sizeof(bool) * size, sizeof(bool)),
                       capsule_owner);
}

#endif