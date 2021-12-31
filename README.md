# Metagenomic-DeepFRI

## About The Project

### Built With

* [DeepFRI](https://github.com/flatironinstitute/DeepFRI)
* [MMseqs2](https://github.com/soedinglab/MMseqs2)
* [Docker](https://www.docker.com/)

##Installation
### Docker
```
docker run -it -u $(id -u):$(id -g) -v /data:/data FUTURE_DOCKER_NAME
```

### Local setup
1. Setup python environment
    ```
    pip install .  
    ```
2. Install mmseqs2 
    ```
    sudo apt install mmseqs2
   ```
3. Install boost libraries
    ```
    sudo apt-get install libboost-numpy1.71 libboost-python1.71
   ```
4. Edit `CONFIG.py` to set up your folder structure
5. Run `post_setup.py` script to create folders and download DeepFRI weights
   ```
   python post_setup.py
   ```

## Usage

1. Upload structure files to specific folder
2. Run `create_mmaseqs_database.py` script
3. Upload `**/*.faa` files into query folder
4. Run `main_pipeline.py`
5. Collect results from finished folder

## Contributing

If you have a suggestion that would make this project better, fork the repo and create a pull request or email me.

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## Contact

Piotr Kucharski - soliareofastorauj@gmail.com
