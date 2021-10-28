# Basecalling Instructions

Basecalling steps:
- Setup a Linode GPU 
- Download data
- Clone bioliquid-nanopore repository
- Setup Docker
- Run guppy
- Run qualimap

## Setup GPU
- Create a Linode GPU VM with Ubuntu 18.04, 8 Cores, 32 GB Ram.
- Check GPU is available
``` 
ssh root@[IP-ADDRESS]
lspci -vnn | grep NVIDIA
```
- Harden the server with the usual instructions.
Copy the SSH key to the server.
```
scp ~/.ssh/id_rsa_linode.pub root@[IP-ADRESS]:~/.ssh/authorized_keys
```
Harden the server
```
ssh root@[IP-ADDRESS]
apt update && apt -y upgrade
apt -y install git make
git clone --single-branch --branch ubuntu https://github.com/yufernando/dotfiles ~/.dotfiles
cd ~/.dotfiles
make all host=aretian-genomics-gpu user=fer password=[PASSWORD]
```

- Install NVIDIA driver dependencies
```
ssh fer@[IP-ADDRESS]
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
mkdir ~/genomics/fast5
cd ~/genomics/fast5
s3cmd get -r [data_filepath]
```

## Clone Bioliquid Nanopore repository
```
cd ~/genomics
git clone https://github.com/aretiandev/bioliquid-nanopore.git
```

There should now be two folders in `~/genomics`: `bioliquid-nanopore/` and `fast5/`

## Install Guppy GPU

cd /tmp && \
    wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_5.0.7_linux64.tar.gz && \
    tar -xzvf ont-guppy_5.0.7_linux64.tar.gz && \
    rm /tmp/ont-guppy_5.0.7_linux64.tar.gz && \
    mv ont-guppy /opt/ont-guppy && \ 
    export PATH="${PATH}:/opt/ont-guppy/bin"

## Run Guppy Basecaller

You can run Guppy from the Makefile but it is more transparent to use the bash script. From within the docker container:
```
cd ~/genomics/bioliquid-nanopore
bash 00_basecaller.sh gpu /home/fer/genomics/fast5 /home/fer/genomics/basecall-latest
```
You can monitor GPU usage with `nvidia-smi`.

Concatenate all files
```
cat /home/fer/genomics/basecall-latest/pass/*fastq > /home/fer/genomics/basecall-latest/bioliquid_run4.fastq
```

## Align, index and sort
Install samtools
```
sudo apt install samtools
```
Align, index and sort
```
cd ~/genomics/bioliquid-nanopore
make align run=4
```

## Run Qualimap
```
cd ~/genomics
docker run --rm -v $PWD:/data pegi3s/qualimap qualimap bamqc -bam /data/bioliquid_run1.bam
```

If you get an “Out of memory error” try:
````
docker run --rm -v $PWD:/data pegi3s/qualimap qualimap bamqc -bam /data/bioliquid_run1.bam --java-mem-size=50G
```

Where 50G is the max memory. Try using 70-80% of your total memory.
