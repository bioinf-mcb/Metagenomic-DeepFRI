//
// Created by soliareofastora on 04.01.2022.
//

#ifndef LOAD_CONTACT_MAPS
#define LOAD_CONTACT_MAPS

#include <cmath>
#include <iostream>
#include <queue>

#include "atoms_file_io.h"
#include "python_utils.h"

static float Distance(float* array, int i, int j) {
  return sqrtf(powf(array[i * 3] - array[j * 3], 2) + powf(array[i * 3 + 1] - array[j * 3 + 1], 2) + powf(array[i * 3 + 2] - array[j * 3 + 2], 2));
}


static std::pair<bool*, int> LoadDenseContactMap(const std::string& file_path, const float angstrom_contact_threshold){
  int chain_length;
  int* group_indexes;
  float* atoms_positions;
  std::tie(chain_length, group_indexes, atoms_positions) = LoadAtomsFile(file_path);

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
          if (Distance(atoms_positions, atom_a, atom_b) <= angstrom_contact_threshold) {
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

  delete[] group_indexes;
  delete[] atoms_positions;
  return std::make_pair(output_data, chain_length);
}


static std::vector<std::pair<int, int>>* LoadSparseContactMap(const std::string& file_path, const float angstrom_contact_threshold){
  int chain_length;
  int* group_indexes;
  float* atoms_positions;
  std::tie(chain_length, group_indexes, atoms_positions) = LoadAtomsFile(file_path);

  // fill up vector with sparse atom contacts
  auto* sparse_contacts = new std::vector<std::pair<int, int>>();
  sparse_contacts->reserve(chain_length * 10);

  for (int group_a = 0; group_a < chain_length; ++group_a) {
    for (int group_b = group_a + 1; group_b < chain_length; ++group_b) {
      bool group_connected = false;

      for (int atom_a = group_indexes[group_a]; atom_a < group_indexes[group_a + 1]; ++atom_a) {
        for (int atom_b = group_indexes[group_b]; atom_b < group_indexes[group_b + 1]; ++atom_b) {
          if (Distance(atoms_positions, atom_a, atom_b) <= angstrom_contact_threshold) {
            group_connected = true;
            sparse_contacts->emplace_back(group_a, group_b);
            break;
          }
        }

        if (group_connected)
          break;
      }
    }
  }

  delete[] group_indexes;
  delete[] atoms_positions;
  return sparse_contacts;
}


static np::ndarray LoadContactMap(const std::string& file_path, const float angstrom_contact_threshold) {
  bool* contact_map;
  int chain_length;
  std::tie(contact_map, chain_length) = LoadDenseContactMap(file_path, angstrom_contact_threshold);
  return CreateNumpyArray(contact_map, chain_length);
}


static np::ndarray LoadAlignedContactMap(const std::string& file_path, float angstrom_contact_threshold, const std::string& query_alignment, const std::string& target_alignment, const int generated_contacts) {
  std::vector<std::pair<int, int>>* sparse_target_contacts = LoadSparseContactMap(file_path, angstrom_contact_threshold);
  std::vector<std::pair<int, int>> sparse_query_contacts;
  sparse_query_contacts.reserve(sparse_target_contacts->size());

  int target_index = 0;
  int query_index = 0;
  std::vector<int> target_to_query_indexes;
  target_to_query_indexes.reserve(query_alignment.size());

  for (int i = 0; i < query_alignment.size(); ++i) {
    if (query_alignment[i] == '-') {
      target_to_query_indexes[target_index] = -1;
      ++target_index;
    } else if (target_alignment[i] == '-') {
      for (int j = 1; j < generated_contacts + 1; ++j) {
        sparse_query_contacts.emplace_back(query_index - j, query_index);
        sparse_query_contacts.emplace_back(query_index + j, query_index);
      }
      ++query_index;
    } else {
      target_to_query_indexes[target_index] = query_index;
      ++target_index;
      ++query_index;
    }
  }

  for (int i = 0; i < sparse_target_contacts->size(); ++i) {
    int contact_x = target_to_query_indexes[sparse_target_contacts->operator[](i).first];
    if (contact_x < 0) {
      continue;
    }
    int contact_y = target_to_query_indexes[sparse_target_contacts->operator[](i).second];
    if (contact_y < 0) {
      continue;
    }
    sparse_query_contacts.emplace_back(contact_x, contact_y);
  }

  bool* const output_data = new bool[(int) pow(query_index, 2)];
  std::memset(output_data, 0, (int) pow(query_index, 2));

  for (int i = 0; i < query_index; ++i) {
    output_data[i * query_index + i] = true;
  }
  for (std::pair<int, int> pair : sparse_query_contacts) {
    if (pair.first < 0) {
      continue;
    }
    if (pair.first >= query_index) {
      continue;
    }
    output_data[pair.first * query_index + pair.second] = true;
    output_data[pair.second * query_index + pair.first] = true;
  }

  delete sparse_target_contacts;
  return CreateNumpyArray(output_data, query_index);
}

#endif
