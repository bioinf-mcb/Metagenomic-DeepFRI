//
// Created by soliareofastora on 27.05.2021.
//

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <cmath>
#include <fstream>
#include <iostream>

namespace py = boost::python;
namespace np = py::numpy;

static inline float _Distance(float* array, int i, int j) {
  return sqrtf(powf(array[i * 3] - array[j * 3], 2) + powf(array[i * 3 + 1] - array[j * 3 + 1], 2) + powf(array[i * 3 + 2] - array[j * 3 + 2], 2));
}

static inline void _DestroyCapsule(PyObject* self) {
  auto* b = reinterpret_cast<bool*>(PyCapsule_GetPointer(self, nullptr));
  delete[] b;
}

// Always run this function before using other functions!!
static void Initialize() {
  Py_Initialize();
  boost::python::numpy::initialize();
}

static void SaveAtoms(const np::ndarray &position_array, const np::ndarray &groups_array, const std::string &save_path) {
  int chain_length = (int) groups_array.shape(0);
  int* group_indexes = reinterpret_cast<int*>(groups_array.get_data());
  int atom_count = group_indexes[chain_length - 1];
  float* atom_positions = reinterpret_cast<float*>(position_array.get_data());

  std::ofstream writer(save_path, std::ios::out | std::ios::binary);
  writer.write(reinterpret_cast<const char*>(&chain_length), 4);
  writer.write(reinterpret_cast<const char*>(group_indexes), 4 * chain_length);
  writer.write(reinterpret_cast<const char*>(atom_positions), 4 * atom_count * 3);
  writer.close();
}

static np::ndarray LoadContactMap(const std::string& file_path, float angstrom_contact_threshold=6) {
  // parse file
  std::ifstream reader(file_path, std::ios::in | std::ios::binary);

  int chain_length;
  reader.read(reinterpret_cast<char*>(&chain_length), 4);

  int* group_indexes = new int[chain_length];
  reader.read(reinterpret_cast<char*>(group_indexes), 4 * chain_length);

  int atom_count = group_indexes[chain_length - 1];

  float* atoms_positions = new float[atom_count * 3];
  reader.read(reinterpret_cast<char*>(atoms_positions), 4 * atom_count * 3);
  reader.close();
  chain_length--;

  // allocate output array
  bool* const output_data = new bool[(int) pow(chain_length, 2)];
  std::memset(output_data, 0, (int) pow(chain_length, 2));

  // fill up output array with atom contacts
  for (int group_a = 0; group_a < chain_length; ++group_a) {
    output_data[group_a * chain_length + group_a] = true;

    for (int group_b = group_a + 1; group_b < chain_length; ++group_b) {
      bool group_connected = false;

      for (int atom_a = group_indexes[group_a]; atom_a < group_indexes[group_a + 1]; ++atom_a) {
        for (int atom_b = group_indexes[group_b]; atom_b < group_indexes[group_b + 1]; ++atom_b) {
          if (_Distance(atoms_positions, atom_a, atom_b) <= angstrom_contact_threshold) {
            group_connected = true;
            output_data[group_a * chain_length + group_b] = true;
            output_data[group_a + group_b * chain_length] = true;
            break;
          }
        }

        if (group_connected)
          break;
      }
    }
  }

  // workaround memory leak in numpy arrays inside python
  // https://stackoverflow.com/questions/57068443/setting-owner-in-boostpythonndarray-so-that-data-is-owned-and-managed-by-pyt
  PyObject* capsule = PyCapsule_New((void*) output_data, nullptr, (PyCapsule_Destructor) &_DestroyCapsule);
  py::handle<> capsule_handle{capsule};
  py::object capsule_owner{capsule_handle};

  return np::from_data(output_data,
                       np::dtype::get_builtin<bool>(),
                       py::make_tuple(chain_length, chain_length),
                       py::make_tuple(sizeof(bool) * chain_length, sizeof(bool)),
                       capsule_owner);
}

BOOST_PYTHON_MODULE (libAtomDistanceIO) {
  py::def("initialize", Initialize);
  py::def("save_atoms", SaveAtoms);
  py::def("load_contact_map", LoadContactMap);
}
