#include "bit_set.hpp"
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <cmath>
#include <fstream>
#include <iostream>

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

  float distance(float* array, int i, int j){
    return sqrtf(powf(array[i*3] - array[j*3],2)+powf(array[i*3+1] - array[j*3+1],2)+powf(array[i*3+2] - array[j*3+2],2));
  }

  void fun(const boost::python::numpy::ndarray &position_array, const boost::python::numpy::ndarray &groups_array, std::string name){

    int atom_size = (int)position_array.shape(0);
    int seq_size = (int)groups_array.shape(0);
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

    auto atom_positions = reinterpret_cast<float*>(position_array.get_data());
    auto group_counts = reinterpret_cast<int*>(groups_array.get_data());

    int group_i = 0;

    for (int i = 0; i < atom_size; ++i) {
      if (i== group_counts[group_i]){
        ++group_i;
      }
      if (group_i == seq_size)
        break;

      int group_j = group_i+1;

      for (int j = group_counts[group_i]; j < atom_size; ++j) {
        if (j == group_counts[group_j]){
          ++group_j;
        }
        if (distance(atom_positions,i,j) <= 6.0f){
          int bit_index = int(bits_size - (seq_size -group_i)*((seq_size -group_i)-1)/2 + group_j - group_i - 1);
          bit_set->set_bit(bit_index, true);
          j = group_counts[group_j];
          ++group_j;
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

BOOST_PYTHON_MODULE (ContactMapper){
  boost::python::class_<ContactMapper, boost::shared_ptr<ContactMapper>,
      boost::noncopyable>("contact_mapper", boost::python::no_init)
      .def("__init__", boost::python::make_constructor(&create_mapper))
      .def("fun", &ContactMapper::fun)
  ;
}