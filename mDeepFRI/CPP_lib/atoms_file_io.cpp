//
// Created by soliareofastora on 24.12.2021.
//
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <fstream>

namespace py = boost::python;
namespace np = py::numpy;

static void SaveAtomsFile(const np::ndarray &position_array,
                          const np::ndarray &groups_array,
                          const std::string &save_path) {
  size_t chain_length = (size_t)groups_array.shape(0);
  uint32_t *group_indexes =
      reinterpret_cast<uint32_t *>(groups_array.get_data());
  size_t atom_count = group_indexes[chain_length - 1];
  float *atom_positions = reinterpret_cast<float *>(position_array.get_data());

  std::ofstream writer(save_path, std::ios::out | std::ios::binary);
  writer << chain_length << std::endl;
  writer << group_indexes << std::endl;
  writer << atom_positions << std::endl;

  writer.close();
};

static std::tuple<size_t, std::unique_ptr<size_t[]>, std::unique_ptr<float[]>>
LoadAtomsFile(const std::string &filepath) {
  std::ifstream reader(filepath, std::ios::in | std::ios::binary);

  size_t chain_length;
  reader.read(reinterpret_cast<char *>(&chain_length), 4);

  std::unique_ptr<size_t[]> group_indexes =
      std::make_unique<size_t[]>(chain_length);
  reader.read(reinterpret_cast<char *>(group_indexes.get()), 4 * chain_length);

  chain_length--;
  size_t atom_count = group_indexes[chain_length];

  std::unique_ptr<float[]> atoms_positions =
      std::make_unique<float[]>(atom_count * 3);
  reader.read(reinterpret_cast<char *>(atoms_positions.get()),
              4 * atom_count * 3);
  reader.close();

  return std::make_tuple(chain_length, std::move(group_indexes),
                         std::move(atoms_positions));
};
