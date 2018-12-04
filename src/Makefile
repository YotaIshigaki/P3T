#PS_PATH = ../FDPS-5.0/src/
PS_PATH = ../../FDPS-5.0/src/ 
INC = -I$(PS_PATH)

#CC = clang++ -std=c++11
#CC = time CC
CC = mpic++ -std=c++11
#CC = /usr/local/mpich2-1.4.1p1-gnu64-4.4.7/bin/mpic++ 
CFLAGS = -O3
CFLAGS += -Wall
CFLAGS += -ffast-math
CFLAGS += -funroll-loops
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CFLAGS += -march=skylake-avx512
CFLAGS += -mavx2
#CFLAGS += -mavx512f -mavx512dq

CFLAGS += -DUSE_QUAD
#CFLAGS += -DCALC_EP_64bit
#CFLAGS += -DPARALLEL_I8J2
#CFLAGS += -DCUTOFF_0
CFLAGS += -DCOLLISION
#CFLAGS += -DUSE_INDIVIDUAL_RADII
CFLAGS += -DOUTPUT_DETAIL
#CFLAGS += -DCHAMBERS
#CFLAGS += -DGAS_DRAG
CFLAGS += -DCALC_WTIME
CFLAGS += -DMONAR
#CFLAGS += -DFORDEBUG
CFLAGS += -g


SRC = main_p3t.cpp
PROGRAM = $(SRC:%.cpp=%.out)

$(PROGRAM):$(SRC)
	$(CC) $(INC) $(CFLAGS) -o $@ $<

clean:
	rm -f $(PROGRAM) *.o *~

