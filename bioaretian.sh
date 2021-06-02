#!/usr/bin/env bash

# bioaretian.sh
#
# This scripts runs the bioaretian image inside a docker container and exposes
# it on port 8888
# ARGS: MODE: gpu or blank: if using guppy_basecaller gpu version

MODE=$1

if ["$#" -eq 0 ];
then 
    docker run --rm -p 8888:8888 -v /home/fer/nanopore:/home/jovyan/work -e JUPYTER_ENABLE_LAB=yes -e GRANT_SUDO=yes --user root --name bioaretian yufernando/bioaretian
elif [ $MODE == "gpu" ];
then
    echo "Running the bioaretian image for a GPU architecture."
    echo ""
    docker run --rm -p 8888:8888 -v /home/fer/nanopore:/home/jovyan/work -e JUPYTER_ENABLE_LAB=yes -e GRANT_SUDO=yes --user root --name bioaretian yufernando/bioaretian:guppy-gpu
else
    docker run --rm -p 8888:8888 -v /home/fer/nanopore:/home/jovyan/work -e JUPYTER_ENABLE_LAB=yes -e GRANT_SUDO=yes --user root --name bioaretian yufernando/bioaretian
fi

# Run in the background
# docker run -d --rm -p 8888:8888 -v /home/fer/nanopore:/home/jovyan/work -e JUPYTER_ENABLE_LAB=yes -e GRANT_SUDO=yes --user root --name bioaretian yufernando/bioaretian
