#PBS -l nodes=1:ppn=48
#PBS -l walltime=1:00:00
#PBS -q cheme-tp
#PBS -m abe
#PBS -r n

# run as: make && qsub run.sh

cd $PBS_O_WORKDIR
source init.sh
mpirun -np 12 segregation2D
