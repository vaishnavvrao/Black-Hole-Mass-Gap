#!/bin/bash
#PBS -N mesa_PPISN
#PBS -q default
#PBS -l nodes=1:ppn=16
#PBS -o out.log
#PBH -e err.log
#PBS -S /bin/bash
#PBS -t 40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90

#Set variables
export MESA_DIR="/home/BHMG2/mesa-r12778"
export OMP_NUM_THREADS=16
export MESA_BASE="/home/BHMG2/ppisn_sim/many_jobs_old/base_axion"
export MESA_INLIST="$MESA_BASE/inlist"
export MESA_RUN="/home/BHMG2/ppisn_sim/many_jobs_old/runs_axion"

#CD to folder
cd $MESA_RUN/$PBS_ARRAYID

#Run MESA
$MESA_BASE/star
