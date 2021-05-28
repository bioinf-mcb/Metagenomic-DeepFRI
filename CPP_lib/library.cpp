#include "bit_set.hpp"
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <cmath>
#include <fstream>
#include <iostream>

static inline float distance(float* array, int i, int j){
  return powf(array[i*3] - array[j*3],2) + powf(array[i*3+1] - array[j*3+1],2) + powf(array[i*3+2] - array[j*3+2],2);
}

class ContactMapper{
 private:
  BitSet* bit_set;

 public:
  ContactMapper(){
    bit_set = new BitSet(100);
    Py_Initialize();
    boost::python::numpy::initialize();
  }

  ~ContactMapper(){
    delete bit_set;
  }

  void Destroy(){
    delete this;
  }

  void GenerateContactMap(const boost::python::numpy::ndarray &position_array, const boost::python::numpy::ndarray &groups_array, std::string name){
    int seq_size = (int)groups_array.shape(0) - 1;
    int bits_size = (int)(seq_size *(seq_size -1) / 2);
    int bytes_size = bits_size/8;
    if (bits_size % 8 > 0) {
      ++bytes_size;
    }

    if (bit_set->size < bytes_size){
      delete bit_set;
      bit_set = new BitSet(bits_size);
    }
    std::memset(bit_set->data, 0, bytes_size);

    auto atom_position = reinterpret_cast<float*>(position_array.get_data());
    auto start_group_index = reinterpret_cast<int *>(groups_array.get_data());

    for (int group_a = 0; group_a < seq_size; ++group_a) {
      for (int group_b = group_a + 1; group_b < seq_size; ++group_b) {

        bool group_connected = false;

        for (int atom_i = start_group_index[group_a]; atom_i < start_group_index[group_a + 1]; ++atom_i) {
          for (int atom_j = start_group_index[group_b]; atom_j < start_group_index[group_b + 1]; ++atom_j) {
            if (distance(atom_position, atom_i, atom_j) <= 36.0f) {
              group_connected = true;
              auto bit_index = int(bits_size - (seq_size - group_a) * ((seq_size - group_a) - 1) / 2 + group_b - group_a - 1);
              bit_set->set_bit(bit_index);
              break;
            }
          }

          if (group_connected) {
            break;
          }
        }

      }
    }

    std::ofstream writer(name, std::ios::out | std::ios::binary);
    writer.write((char*)bit_set->data, bytes_size);
    writer.close();
  }
};

boost::shared_ptr<ContactMapper> create_mapper(){
    return boost::shared_ptr<ContactMapper>(
        new ContactMapper(),
        boost::mem_fn(&ContactMapper::Destroy));
}

BOOST_PYTHON_MODULE (libContactMapper){
  boost::python::class_<ContactMapper, boost::shared_ptr<ContactMapper>,
      boost::noncopyable>("contact_mapper", boost::python::no_init)
      .def("__init__", boost::python::make_constructor(&create_mapper))
      .def("generate_contact_map", &ContactMapper::GenerateContactMap)
  ;
}