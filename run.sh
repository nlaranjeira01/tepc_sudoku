#!/bin/bash
i=0
csna=12
for maxw in {1750..4750..750}
do
	minw=`expr $maxw - 500`

	for mnpm in {18..50..8}
	do
		mpicc -O3 -std=c11 -march=native -o sudoku sudoku.c minheap.c dynarray.c -DMAX_WORKLOAD=$maxw -DMIN_WORKLOAD=$minw -DMAX_NODES_PER_MESSAGE=$mnpm -DCOORD_SEND_NODES_AMOUNT=$csna -DN_THREADS=3 -lpthread -Wall
		mpirun -np 6 ./sudoku < example_inputs/ex5.in
	done
done
