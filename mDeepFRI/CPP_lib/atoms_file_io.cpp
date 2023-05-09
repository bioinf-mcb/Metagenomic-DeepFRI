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
