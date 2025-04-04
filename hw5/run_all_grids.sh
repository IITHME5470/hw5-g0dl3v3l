#!/bin/bash

# Define grid configurations
grids=("2 2" "2 4" "4 4")

# Loop over each grid setting
for grid in "${grids[@]}"; do
    # Extract X and Y to calculate number of processes
    X=$(echo "$grid" | cut -d' ' -f1)
    Y=$(echo "$grid" | cut -d' ' -f2)
    NUM_PROCS=$((X * Y))

    echo "=========================================="
    echo "Running for grid: $grid with $NUM_PROCS processes"
    
    # Replace last line of input2d.in with new grid
    sed -i '$s/.*/'"$grid"'/' input2d.in

    # Compile and run
    mpicc parhc2d_skel.c -lm
    mpirun -n "$NUM_PROCS" ./a.out

    echo "Finished run for grid $grid"
    echo "=========================================="
done
