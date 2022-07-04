//
// Created by soliareofastora on 24.12.2021.
//

#ifndef ATOMS_FILE_IO
#define ATOMS_FILE_IO

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <fstream>

namespace py = boost::python;
namespace np = py::numpy;


static void SaveAtomsFile(const np::ndarray &position_array, const np::ndarray &groups_array, const std::string &save_path) {
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


static std::tuple<int, int*, float*> LoadAtomsFile(const std::string& file_path){
  std::ifstream reader(file_path, std::ios::in | std::ios::binary);

  int chain_length;
  reader.read(reinterpret_cast<char*>(&chain_length), 4);

  int* group_indexes = new int[chain_length];
  reader.read(reinterpret_cast<char*>(group_indexes), 4 * chain_length);

  chain_length--;
  int atom_count = group_indexes[chain_length];

  float* atoms_positions = new float[atom_count * 3];
  reader.read(reinterpret_cast<char*>(atoms_positions), 4 * atom_count * 3);
  reader.close();

  return std::make_tuple(chain_length, group_indexes, atoms_positions);
}

#endif
