default: all

all:
	mpicc -O3 -std=c11 -march=native -o sudoku sudoku.c minheap.c dynarray.c -Wall

clean:
	rm sudoku
