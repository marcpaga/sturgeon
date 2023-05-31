#! /bin/bash

source ~/nanocns/software/sturgeon/venv/bin/activate;

sturgeon inputtobed --margin 50 -i $1 -o $1 -s guppy --probes-file /home/sturgeon/nanocns/data/probelocs_chm13.bed ;
sturgeon predict -p --i $1 -o $1 -m ~/nanocns/sturgeon_models/Sturgeon_general_V2/model.zip ;
