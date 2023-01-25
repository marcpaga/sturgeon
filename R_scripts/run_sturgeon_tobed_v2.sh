#! /bin/bash

source ~/nanocns/software/sturgeon/venv/bin/activate;

sturgeon inputtobed -i $1 -o $2 -s megalodon ;
sturgeon predict -p --i $2 -o $2 -m ~/nanocns/sturgeon_models/Sturgeon_general_V2/model.zip ;
