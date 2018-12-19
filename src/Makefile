PS_PATH = ../../FDPS-5.0a/src/ 
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
CFLAGS += -march=core-avx2
#CFLAGS += -mavx2
#CFLAGS += -mavx512f -mavx512dq
#CFLAGS += -DFAST_ALL_TO_ALL_FOR_K

CFLAGS += -DUSE_QUAD
#CFLAGS += -DCALC_EP_64bit
#CFLAGS += -DPARALLEL_I4J4
#CFLAGS += -DCUTOFF_0
CFLAGS += -DCOLLISION
#CFLAGS += -DUSE_INDIVIDUAL_RADII
CFLAGS += -DOUTPUT_DETAIL
#CFLAGS += -DCHAMBERS
#CFLAGS += -DGAS_DRAG
CFLAGS += -DCALC_WTIME
CFLAGS += -DMONAR
#CFLAGS += -DFORDEBUG
CFLAGS += -DCHECK_NEIGHBOR
CFLAGS += -g


SRC = main_p3t.cpp
PROGRAM = $(SRC:%.cpp=%.out)
HEADER = mathfunc.h kepler.h energy.h particle.h disk.h gravity.h gravity_kernel.hpp collisionA.h collisionB.h hermite.h hard.h read.h time.h func.h

$(PROGRAM): main_p3t.cpp gravity_kernel.o $(HEADER)
	$(CC) $(INC) $(CFLAGS) -o $@ $< gravity_kernel.o


gravity_kernel.o: gravity_kernel.cpp avx.h particle.h
	$(CC) $(INC) $(CFLAGS) -c $<

clean:
	rm -f *.out *.o *~

