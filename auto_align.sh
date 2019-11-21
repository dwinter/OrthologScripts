#!/bin/bash

name=$(basename $1 _unaligned.fna)
echo "Aligning "${name}
mafft --auto --adjustdirection $1 > aligned/${name}_aligned.fna
