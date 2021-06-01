#!/usr/bin/env bash

# bioaretian.sh
#
# This scripts runs the bioaretian image inside a docker container and exposes
# it on port 8888

docker run --rm -p 8888:8888 -v /home/fer/nanopore:/home/jovyan/work -e JUPYTER_ENABLE_LAB=yes -e GRANT_SUDO=yes --user root --name bioaretian yufernando/bioaretian

# Run in the background
# docker run -d --rm -p 8888:8888 -v /home/fer/nanopore:/home/jovyan/work -e JUPYTER_ENABLE_LAB=yes -e GRANT_SUDO=yes --user root --name bioaretian yufernando/bioaretian
