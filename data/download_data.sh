#!/bin/sh

mkdir current
cd current
# AlphaFold 2: if you wish to run TMvis for another proteome, please change link below
wget https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v3.tar

# TMbed: if you wish to run TMvis for another proteome, please go to https://github.com/Rostlab/TMbed and run TMbed
wget https://rostlab.org/public/tmbed/predictions/human_210422_tmbed.tar.gz
