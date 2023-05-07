# CPP_lib
Here you can find a few performance critical functions used in metagenomic_deepfri_pipeline

* `library_definition` contains functions definitions that are accessible in python using BOOST
* In `atoms_file_io` you can find how protein structures are stored in binary format
* `python_utils` implements quite interesting logic of [handling ownership of memory to python](https://stackoverflow.com/questions/57068443/setting-owner-in-boostpythonndarray-so-that-data-is-owned-and-managed-by-pyt)
* `load_contact_maps` contains the most interesting functions

### Build from source
For more information please refer to the [cmake documentation](https://cmake.org/runningcmake/).

First we need to install boost
```bash
wget https://boostorg.jfrog.io/artifactory/main/release/1.81.0/source/boost_1_81_0.tar.gz
tar -xf boost_1_81_0.tar.gz
cd boost_1_81_0
./bootstrap.sh
./b2
```

After building boost you will have to copy and paste include path.
Additional libraries will be installed in `boost_1_81_0/stage/lib` directory.
You will need to find and paste path to `libboost_python310.so.1.81.0`
and `libboost_numpy310.so.1.81.0` files.

```bash
cd ../meta_deepFRI/CPP_lib
ccmake .
make
```
