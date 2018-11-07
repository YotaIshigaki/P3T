PS_PATH = ../../FDPS-4.0a/src/ 
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
#CFLAGS += -mavx2
CFLAGS += -mavx
CFLAGS += -DCOLLISION
CFLAGS += -DUSE_QUAD
CFLAGS += -DUSE_INDIVIDUAL_RADII
CFLAGS += -DOUTPUT_DETAIL
#CFLAGS += -DCHAMBERS
#CFLAGS += -DGAS_DRAG
CFLAGS += -DMONAR
#CFLAGS += -DFORDEBUG
CFLAGS += -g


SRC = main_p3t.cpp
PROGRAM = $(SRC:%.cpp=%.out)

$(PROGRAM):$(SRC)
	$(CC) $(INC) $(CFLAGS) -o $@ $<

clean:
	rm -f $(PROGRAM) *.o *~

