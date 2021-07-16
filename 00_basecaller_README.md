# Basecaller

Instructions to setup a Linode GPU to run basecalling with guppy_basecaller.

# Linode GPU Setup

- Create a GPU VM with Ubuntu 18.04
- Check GPU is available
``` 
ssh root@69.164.211.148
lspci -vnn | grep NVIDIA
```
- Harden the server with the usual instructions.
- Install NVIDIA driver dependencies
```
sudo apt-get install build-essential
wget https://developer.download.nvidia.com/compute/cuda/11.3.1/local_installers/cuda_11.3.1_465.19.01_linux.run
sudo sh cuda_11.3.1_465.19.01_linux.run
```
- Select only the NVIDIA driver out of the 5 installation options
- Check driver installation successful
```
nvidia-smi
```

# Data, Github repository and Docker image

### Clone Bioliquid Nanopore repository
```
cd ~
git clone https://github.com/aretiandev/bioliquid-nanopore.git
```

### Get data
Data should be inside a fast5 folder in the data folder.
```
mkdir ~/data
cd ~/data
s3cmd get -r [data_filepath]
```

### Docker
- Install Docker
- Pull image
``` 
docker pull yufernando/bioaretian:guppy-gpu
```
- Run image
```
cd ~
docker run -d -v "$PWD":/home/jovyan/work -e GRANT_SUDO=yes --user root --name bioaretian yufernando/bioaretian:guppy-gpu
```

- Step into container if needed
```
docker exec -it bioaretian /bin/bash
```

# Samtols
```
sudo apt install samtools
```

# Guppy

Run guppy basecaller from the Makefile recipe.