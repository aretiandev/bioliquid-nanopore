# Setup Instructions

Instructions to 
- Setup a Linode GPU to run basecalling with guppy_basecaller.
- Setup a remote Jupyterlab instance.
- Install Python, R and Strspy.

# Setup a Linode GPU to run basecalling

## Setup GPU
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

## Get data
Data should be inside a fast5 folder in the data folder.
```
mkdir ~/data
cd ~/data
s3cmd get -r [data_filepath]
```

## Clone Bioliquid Nanopore repository
```
cd ~
git clone https://github.com/aretiandev/bioliquid-nanopore.git
```

## Setup Docker
- Install Docker
```
sudo apt-get update
sudo apt-get install \
    ca-certificates \
    curl \
    gnupg \
    lsb-release
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu \
  $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io
sudo groupadd docker
sudo usermod -aG docker fer
newgrp docker 
```

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

### Samtols
```
sudo apt install samtools
```

## Run Guppy Basecaller

Run guppy basecaller from the Makefile recipe.

## Run Qualimap
```
docker run --rm -v $PWD:/data pegi3s/qualimap qualimap bamqc -bam /data/bioliquid_run1.bam
```
If you get an “Out of memory error” try:
````
docker run --rm -v $PWD:/data pegi3s/qualimap qualimap bamqc -bam /data/bioliquid_run1.bam --java-mem-size=50G
```
Where 50G is the max memory. Try using 70-80% of your total memory.

# Setup remote JupyterLab

It is useful to run Jupyterlab from the server and connect to the port through ssh. From the remote server, run Jupyterlab and get the token:
```
docker run -d --rm -v $PWD:/home/jovyan/work -v /mnt:/mnt -p 8888:8888 --name bioaretian yufernando/bioaretian
docker logs bioaretian
```

From the local computer, open an ssh tunnel that connects port 8888 in the server with port 8888 in the local machine:
```
ssh -f -N -L 8888:localhost:8888 fer@aretian-genomics
```

Now open a browser in http://127.0.0.1:8888/ and use the token above or copy and paste the link from `docker logs`.

# Python
conda install pandas numpy
conda install -c conda-forge scikit-learn 

# R
sudo apt-get install libssl-dev

install.packages("tidyverse")
install.packages("dplyr")
install.packages("reshape2")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsamtools")
install.packages("FactoMineR")
install.packages("factoextra")

# Strspy
Clone repo
Follow installation instructions: https://github.com/unique379r/strspy
