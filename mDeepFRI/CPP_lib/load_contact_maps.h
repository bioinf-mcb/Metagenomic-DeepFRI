//
// Created by soliareofastora on 04.01.2022.
//

#ifndef LOAD_CONTACT_MAPS
#define LOAD_CONTACT_MAPS

#include <cmath>
#include <iostream>
#include <queue>
#include <string>
#include <sys/stat.h>

#include "atoms_file_io.h"

std::pair<bool *, int> LoadDenseContactMap(std::string &file_path,
                                           float angstrom_contact_threshold);

std::vector<std::pair<int, int>> *
LoadSparseContactMap(std::string &file_path, float angstrom_contact_threshold);

// np::ndarray LoadAlignedContactMap(const std::string &file_path,
//                                          float angstrom_contact_threshold,
//                                          const std::string &query_alignment,
//                                          const std::string &target_alignment,
//                                          const int generated_contacts);

#endif
