# CPP_lib
Here you can find a few performance critical functions used in metagenomic_deepfri_pipeline

* `library_definition` contains functions definitions that are accessible in python using BOOST
* In `atoms_file_io` you can find how protein structures are stored in binary format
* `python_utils` implements quite interesting logic of [handling ownership of memory to python](https://stackoverflow.com/questions/57068443/setting-owner-in-boostpythonndarray-so-that-data-is-owned-and-managed-by-pyt)
* `load_contact_maps` contains the most interesting functions

### Build from source 
For more information please refer to the [cmake documentation](https://cmake.org/runningcmake/). 
You will probably have to edit `CMakeLists.txt :13` to set python include path. 
```
sudo snap install cmake
ccmake .
make
```