#!/bin/bash

# pe request
#$ -pe mpich 64

# our Job name 
#$ -N antelope

#$ -S /bin/bash

#$ -q dque_ib # $ -V

#$ -cwd

# needs in 
#   $NSLOTS          
#       the number of tasks to be used
#   $TMPDIR/machines 
#       a valid machiche file to be passed to mpirun 
#   enables $TMPDIR/rsh to catch rsh calls if available

echo "Got $NSLOTS slots."
cat $TMPDIR/machines



#######################################################
### mpich 1.2.7p1 (w/ Intel-12.1 compiler)
#######################################################
#
# MPI_HOME=/opt/mpi/intel-12.1/mpich-1.2.7p1
# MPI_EXEC=$MPI_HOME/bin/mpirun
#
# cd $SGE_O_WORKDIR
#
# $MPI_EXEC -machinefile $TMPDIR/machines -np $NSLOTS $SGE_O_WORKDIR/exec_file



#######################################################
### mpich2 1.4.1p1 (w/ Intel-12.1 compiler)
#######################################################
#
# MPI_HOME=/opt/mpi/intel-12.1/mpich2-1.4.1p1
# MPI_EXEC=$MPI_HOME/bin/mpirun

# cd $SGE_O_WORKDIR

# $MPI_EXEC -machinefile $TMPDIR/machines -n $NSLOTS $SGE_O_WORKDIR/exec_file



#######################################################
### openmpi 1.4.4 (w/ Intel-12.1 compiler)
#######################################################
#
# MPI_HOME=/opt/mpi/intel-12.1/openmpi-1.4.4
# MPI_EXEC=$MPI_HOME/bin/mpirun
#
# cd $SGE_O_WORKDIR
#
# $MPI_EXEC -machinefile $TMPDIR/machines -n $NSLOTS /opt/vasp/vasp5.3.3/vasp


#######################################################
### openmpi 1.10.1 (w/ Intel-14.0.2 compiler)
#######################################################
#
 #MPI_HOME=/opt/mpi/intel-14.0/openmpi-1.10.1
 #MPI_HOME=/opt/mpi/intel-13.1/openmpi-1.10.1
 MPI_HOME=/opt/mpi/gcc-4.9.2/openmpi-1.10.1
 MPI_EXEC=$MPI_HOME/bin/mpirun
#
 cd $SGE_O_WORKDIR
#
 $MPI_EXEC -machinefile $TMPDIR/machines -n $NSLOTS ./exec_hhg_mpi > stdout.dat
