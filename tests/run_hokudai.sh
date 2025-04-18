#!/bin/sh 
#------ pjsub option --------#
#PJM -L rscgrp=n22240a
#PJM -L node=1
#PJM -L elapse=24:00:00 
#PJM -g n22240
#PJM -j
#------- Program execution -------#
NUM_NODES=${PJM_VNODES}
NUM_CORES=40
NUM_PROCS=`expr $NUM_NODES "*" $NUM_CORES`

module load intel

PRG="${HOME}/vasp/vasp.6.4.3/bin/vasp_std"
export VASP_PP_PATH="${HOME}/vasp/potentials/"
export ASE_VASP_COMMAND="mpiexec.hydra -n ${NUM_PROCS} ${PRG}"

python test_phase_diagram.py

