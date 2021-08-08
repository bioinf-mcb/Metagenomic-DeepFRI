#include "bit_set.hpp"
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <cmath>
#include <fstream>
#include <iostream>

namespace py = boost::python;
namespace np = py::numpy;

static inline float distance(float* array, int i, int j){
  return powf(array[i*3] - array[j*3],2) + powf(array[i*3+1] - array[j*3+1],2) + powf(array[i*3+2] - array[j*3+2],2);
}

inline void DestroyCapsule(PyObject* self) {
  auto * b = reinterpret_cast<bool*>( PyCapsule_GetPointer(self, nullptr) );
  delete [] b;
}

static void SaveAtoms(const np::ndarray &position_array, const np::ndarray &groups_array, std::string save_path){
  int groups_size = (int)groups_array.shape(0);
  int* groups_index = reinterpret_cast<int *>(groups_array.get_data());
  int atom_count = groups_index[groups_size-1];
  float* atom_positions = reinterpret_cast<float*>(position_array.get_data());

  std::ofstream writer(save_path, std::ios::out | std::ios::binary);
  writer.write(reinterpret_cast <const char *>(&groups_size), 4);
  writer.write(reinterpret_cast <const char *>(groups_index), 4 * groups_size);
  writer.write(reinterpret_cast <const char *>(atom_positions), 4 * atom_count * 3);
  writer.close();
}

static np::ndarray LoadContactMap(const std::string& load_path, float angstrom_contact_threshold=6) {
  std::ifstream reader(load_path, std::ios::in | std::ios::binary);
  int groups_size;
  reader.read(reinterpret_cast <char *>(&groups_size), 4);
  int* groups_index = new int[groups_size];
  reader.read(reinterpret_cast <char *>(groups_index), 4 * groups_size);
  int atom_count = groups_index[groups_size-1];
  float* atom_positions = new float[atom_count * 3];
  reader.read(reinterpret_cast <char *>(atom_positions), 4 * atom_count * 3);
  reader.close();
  int seq_size = groups_size - 1;

  bool * const output_data = new bool[(int)pow(seq_size,2)];
  std::memset(output_data, 0, (int)pow(seq_size,2));

  float distance_threshold = powf(angstrom_contact_threshold, 2);

  for (int group_a = 0; group_a < seq_size; ++group_a) {
    output_data[group_a * seq_size + group_a] = true;
    for (int group_b = group_a + 1; group_b < seq_size; ++group_b) {
      bool group_connected = false;

      for (int atom_a = groups_index[group_a]; atom_a < groups_index[group_a + 1]; ++atom_a) {
        for (int atom_b = groups_index[group_b]; atom_b < groups_index[group_b + 1]; ++atom_b) {
          if (distance(atom_positions, atom_a, atom_b) <= distance_threshold) {
            group_connected = true;
            output_data[group_a * seq_size + group_b] = true;
            output_data[group_a + group_b * seq_size] = true;
            break;
          }
        }
        if (group_connected)
          break;
      }

    }
  }

  // https://stackoverflow.com/questions/57068443/setting-owner-in-boostpythonndarray-so-that-data-is-owned-and-managed-by-pyt
  PyObject* capsule = ::PyCapsule_New((void*) output_data, nullptr, (PyCapsule_Destructor)&DestroyCapsule);
  py::handle<> capsule_handle{capsule};
  py::object capsule_owner{capsule_handle};

  return np::from_data(output_data,
                       np::dtype::get_builtin<bool>(),
                       py::make_tuple(seq_size, seq_size),
                       py::make_tuple(sizeof(bool) * seq_size, sizeof(bool)),
                       capsule_owner
  );
}

static void Initialize(){
  Py_Initialize();
  boost::python::numpy::initialize();
}

BOOST_PYTHON_MODULE (libAtomDistanceIO){
  py::def("initialize", Initialize);
  py::def("save_atoms", SaveAtoms);
  py::def("load_contact_map", LoadContactMap);
}