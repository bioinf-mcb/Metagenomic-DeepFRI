//
// Created by soliareofastora on 24.12.2021.
//
#include <fstream>
#include <memory>
#include <string>

void SaveAtomsFile(float *position_array, size_t atom_count, int *groups_array,
                   size_t chain_length, const std::string &save_path) {

  std::ofstream writer(save_path, std::ios::out | std::ios::binary);
  writer.write(reinterpret_cast<const char *>(&chain_length), 4);
  writer.write(reinterpret_cast<const char *>(groups_array), 4 * chain_length);
  writer.write(reinterpret_cast<const char *>(position_array),
               4 * atom_count * 3);
  writer.close();
}

std::tuple<size_t, std::unique_ptr<size_t[]>, std::unique_ptr<float[]>>
LoadAtomsFile(const std::string &filepath) {

  std::ifstream reader(filepath, std::ios::in | std::ios::binary);
  if (!reader)
    throw std::runtime_error("Failed to open file: " + filepath);

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
