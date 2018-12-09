default: all

all:
	mpicc -O3 -std=c11 -march=native -o sudoku sudoku.c minheap.c dynarray.c -DMAX_WORKLOAD=4000 -DMIN_WORKLOAD=3500 -DMAX_NODES_PER_MESSAGE=20 -DCOORD_SEND_NODES_AMOUNT=12 -DN_THREADS=3 -lpthread -Wall

clean:
	rm sudoku
