#!/bin/bash
## ECGI
#SBATCH -J ECGI
## Resources: (nodes, procs, tasks, walltime, … etc)
#SBATCH -N 1
#SBATCH -n 18
#SBATCH --exclusive
#SBATCH -t15:00:00
# #  standard output message
#SBATCH -o batches/batch%j.out
# # output error message
#SBATCH -e batches/batch%j.err
module purge
module load compiler/gcc

echo “=====my job informations ====”
echo “Node List: ” $SLURM_NODELIST
echo “my jobID: ” $SLURM_JOB_ID
echo “Partition: ” $SLURM_JOB_PARTITION
echo “submit directory:” $SLURM_SUBMIT_DIR
echo “submit host:” $SLURM_SUBMIT_HOST
echo “In the directory: `pwd`”
echo “As the user: `whoami`”

cd /home/vpanneti/beegfs/niami/code_inverse 
make clean 
make 

./solveHeat 1 1 0.00 & ./solveHeat 1 2 0.00 & ./solveHeat 1 3 0.00 &
./solveHeat 2 1 0.00 & ./solveHeat 2 2 0.00 & ./solveHeat 2 3 0.00 &
./solveHeat 3 1 0.00 & ./solveHeat 3 2 0.00 & ./solveHeat 3 3 0.00 &
./solveHeat 1 1 0.02 & ./solveHeat 1 2 0.02 & ./solveHeat 1 3 0.02 &
./solveHeat 2 1 0.02 & ./solveHeat 2 2 0.02 & ./solveHeat 2 3 0.02 &
./solveHeat 3 1 0.02 & ./solveHeat 3 2 0.02 & ./solveHeat 3 3 0.02