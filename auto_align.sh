#!/bin/bash

name=$(basename $1 _unaligned.fna)
echo ${name}
mafft --auto --adjustdirection $1 > ${name}_aligned
mv ${name}_aligned ./results/aligned
