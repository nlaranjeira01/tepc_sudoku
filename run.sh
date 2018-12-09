#!/bin/bash
nt=3
csna=12

mpicc -O3 -std=c11 -march=native -o sudoku sudoku.c minheap.c dynarray.c -lpthread

for maxw in {1750..4750..750}
do
	minw=`expr $maxw - 500`

	for mnpm in {18..50..8}
	do
		mpirun -np 6 ./sudoku $maxw $minw $mnpm $csna $nt <example_inputs/ex5.in
	done
done
