// #include <cmath>
// #include <iostream>
// #include <queue>
// #include <string>
// #include <sys/stat.h>

// #include "atoms_file_io.h"
// #include "load_contact_maps.h"
// #include "python_utils.h"

// // static np::ndarray LoadAlignedContactMap(const std::string &file_path,
// //                                          float angstrom_contact_threshold,
// //                                          const std::string
// //                                          &query_alignment, const
// //                                          std::string &target_alignment,
// //                                          const int generated_contacts) {
// //   std::vector<std::pair<int, int>> *sparse_target_contacts =
// //       LoadSparseContactMap(file_path, angstrom_contact_threshold);
// //   std::vector<std::pair<int, int>> sparse_query_contacts;
// //   sparse_query_contacts.reserve(sparse_target_contacts->size());

// //   int target_index = 0;
// //   int query_index = 0;
// //   std::vector<int> target_to_query_indexes;
// //   target_to_query_indexes.reserve(query_alignment.size());

// //   for (int i = 0; i < query_alignment.size(); ++i) {
// //     if (query_alignment[i] == '-') {
// //       target_to_query_indexes[target_index] = -1;
// //       ++target_index;
// //     } else if (target_alignment[i] == '-') {
// //       for (int j = 1; j < generated_contacts + 1; ++j) {
// //         sparse_query_contacts.emplace_back(query_index - j, query_index);
// //         sparse_query_contacts.emplace_back(query_index + j, query_index);
// //       }
// //       ++query_index;
// //     } else {
// //       target_to_query_indexes[target_index] = query_index;
// //       ++target_index;
// //       ++query_index;
// //     }
// //   }

// //   for (int i = 0; i < sparse_target_contacts->size(); ++i) {
// //     int contact_x =
// // target_to_query_indexes[sparse_target_contacts->operator[](i).first];
// //     if (contact_x < 0) {
// //       continue;
// //     }
// //     int contact_y =
// // target_to_query_indexes[sparse_target_contacts->operator[](i).second];
// //     if (contact_y < 0) {
// //       continue;
// //     }
// //     sparse_query_contacts.emplace_back(contact_x, contact_y);
// //   }

// //   bool *const output_data = new bool[(int)pow(query_index, 2)];
// //   std::memset(output_data, 0, (int)pow(query_index, 2));

// //   for (int i = 0; i < query_index; ++i) {
// //     output_data[i * query_index + i] = true;
// //   }
// //   for (std::pair<int, int> pair : sparse_query_contacts) {
// //     if (pair.first < 0) {
// //       continue;
// //     }
// //     if (pair.first >= query_index) {
// //       continue;
// //     }
// //     output_data[pair.first * query_index + pair.second] = true;
// //     output_data[pair.second * query_index + pair.first] = true;
// //   }

// //   delete sparse_target_contacts;
// //   return CreateNumpyArray(output_data, query_index);
// // }
