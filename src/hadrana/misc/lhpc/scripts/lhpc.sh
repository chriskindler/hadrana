#!/bin/bash

# h5ls --data A650r000n100_run13.h5/rqcd/seq3pt/meas0/NUCL_U_UNPOL0/seqsrc_id
# seqsrc_id                Dataset {SCALAR}
#     Data:
#         (0) "73b12658da373187bffaeec5dc7f5c266b167687"

# ensemble path specs
ensemble="A650"
base_path="/glurch/scratch/analysis/CLS/${ensemble}"
run_dir="run13_nog2"

# config file specs
run="run13"
rep="r000"
cfg="100"

echo ${base_path}

# FLAVOUR U UNPOLARISED
    echo "+++++++++++++++ POLARISED Y +++++++++++++++ "
for i in {0,1,2,3,8} 
do
    echo "--------------- FLAVOUR U ---------------"
    echo "                MEAS-ID $i               "
    h5ls -d ${base_path}/${ensemble}${rep}/${run_dir}/DATA/A650r000n100_run13.h5/rqcd/seq3pt/meas$i/NUCL_U_POLY0/seqsrc_id
done

# FLAVOUR D UNPOLARISED
for i in {4,5,6,7,9} 
do
    echo "+++++++++++++++ POLARISED Y +++++++++++++++ "
    echo "--------------- FLAVOUR D ---------------"
    echo "                MEAS-ID $i               "
    h5ls -d ${base_path}/${ensemble}${rep}/${run_dir}/DATA/A650r000n100_run13.h5/rqcd/seq3pt/meas$i/NUCL_D_POLY0/seqsrc_id
done
