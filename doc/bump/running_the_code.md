# Running the code

To run the code, go the saber-bundle build directory and type:
 
    cd bin
    export OMP_NUM_THREADS=$NTHREAD
    mpirun -n $NTASK ./bump.x your_namelist

where $NTHREAD is the number of OpenMP threads and $NTASK is the number of MPI tasks that are desired.
