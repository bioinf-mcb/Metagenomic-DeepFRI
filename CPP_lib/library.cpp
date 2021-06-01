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

  void GenerateContactMap(const np::ndarray &position_array, const np::ndarray &groups_array, std::string save_path){
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

    std::ofstream writer(save_path, std::ios::out | std::ios::binary);
    writer.write((char*)bit_set->data, bytes_size);
    writer.close();
  }

  np::ndarray LoadCmap(const std::string& path) {
    std::ifstream reader(path, std::ios::in | std::ios::binary | std::ios::ate);
    auto bytes_size = reader.tellg();
    if (bit_set->size < bytes_size) {
      delete bit_set;
      bit_set = new BitSet(bytes_size*8);
    }
    std::memset(bit_set->data, 0, bytes_size);

    int n = int((1 + sqrt(64 * bytes_size + 1)) / 2);
    reader.seekg(0, std::ios::beg);
    reader.read((char*) bit_set->data, bytes_size);
    reader.close();

    bool * const data = new bool[n * n];
    std::memset(data, 0, n * n);

    for (int i = 0; i < n; ++i) {
      data[i*n + i] = true;
      for (int j = i+1; j < n; ++j) {
        auto bit_index = int(int((n * (n - 1)) / 2) - (n - i) * ((n - i) - 1) / 2 + j - i - 1);
        if (bit_set->get_bit(bit_index)) {
          data[i * n + j] = true;
          data[j * n + i] = true;
        }
      }
    }

    // https://stackoverflow.com/questions/57068443/setting-owner-in-boostpythonndarray-so-that-data-is-owned-and-managed-by-pyt
    PyObject* capsule = ::PyCapsule_New((void*)data, nullptr, (PyCapsule_Destructor)&DestroyCapsule);
    py::handle<> capsule_handle{capsule};
    py::object capsule_owner{capsule_handle};

    return np::from_data(data,
                         np::dtype::get_builtin<bool>(),
                         py::make_tuple(n, n),
                         py::make_tuple(sizeof(bool) * n, sizeof(bool)),
                         capsule_owner
                         );
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
      .def("load_cmap", &ContactMapper::LoadCmap)
  ;
}