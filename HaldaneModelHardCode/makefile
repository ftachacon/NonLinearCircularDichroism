
#CFLAGS= -O3 -openmp
#INTEL= icpc

MPI= /usr/local/bin/mpicxx   # MPI C/C++ COMPIPLER PATH

#LIBS= -lpthread -lm -ldl -liomp5 -lmkl_rt    # OTHER PATHS TO LIBS


#My library a. chacon
MyPath=.
OWN_LIB=${MyPath}/include


ProgramMPI=hhg_rk5_hm_mpi_full_vel_nc.cpp
exnamempi=exec_hhg_mpi


comp: compileMPI 
#comp: compileSq compileMPI 


compileSq: ${Program0} 
	$(MPI)  ${Program0} -w -I${OWN_LIB} -o ${exname}


compileMPI: ${ProgramMPI}
	${MPI} ${ProgramMPI} -w -I${OWN_LIB} -o ${exnamempi}


clean: 
	rm a.out

#ROOT_PATH=/home/achacon
#LD_FFTW=$(ROOT_PATH)/libraries/fftw-3.3.4/LFFTW3


