#! /bin/bash

source /mnt/c/Users/sturgeon/Desktop/sturgeon/venv2/bin/activate
sturgeon inputtobed -i $1 -o $1 -s guppy --probes-file $2
sturgeon predict -i $1 -o $1 -m $3
