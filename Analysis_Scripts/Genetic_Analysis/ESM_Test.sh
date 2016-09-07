#!/bin/bash
#   Script to run the ESM Test on the PLINK2 permuted association test output.
#       PLINK2: https://www.cog-genomics.org/plink2
#       ESMTest: https://github.com/ThorntonLab/ESMtest
#       HDF5 Utils: https://www.hdfgroup.org/downloads/

export PERMOUT="${HOME}/DataDisk/tmp/Perms"
cd ${PERMOUT}

#   Define values for the ESM teset
export WIN_SIZE="4000"
export JUMP_SIZE="1000"
export NMARKERS="50"
export CMARKERS="50"
export CPERMS="1000"
export NPERMS="1500000"
#   Define a function that runs the ESM test on an input file
esm_test() {
    #   The chromosome name is the argument. We will build the input and output
    #   filenames from the chromosome.
    CHR="$1"
    INPUT="${PERMOUT}/${CHR}/${CHR}_Permutations.h5"
    OUTPUT="${PERMOUT}/${CHR}/${CHR}_ESM.txt"
    #   Run the ESM test
    esmk -o ${OUTPUT} -w ${WIN_SIZE} -j ${JUMP_SIZE} -k ${NMARKERS} -r 0.75 -n 1 --cmarkers ${CMARKERS} --cperms ${CPERMS} --nperms ${NPERMS} ${INPUT}
}

export -f esm_test
parallel esm_test ::: 1H 2H 3H 4H 5H 7H
