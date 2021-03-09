/**
 * A sequential Jacobi iteration program for solving partial differential equations.
 *
 * Usage: ./a.out <gridSize> <numIters>
 * Output:  Prints command-line arguments used, execution time and maximum error to stdout.
 *          Prints final grid values to a file 'seq-data.out'
 *
 * Authors: Lucas Villarroel, Erik Hanstad
 * Language: C
*/
#include <stdio.h>
#include <stdlib.h>

#define MAX_ITERATIONS 1000
#define MAX_SIZE 10000

int gridSize, numIters;

void printGrid(double** grid) {
  FILE* fp = fopen("seq-data.out", "w");
  for (int i = 0; i < gridSize; i++) {
    fprintf(fp, "[");
    for (int j = 0; j < gridSize; j++) {
      if (j == 0) {
        fprintf(fp, "%g", grid[i][j]);
      } else {
        fprintf(fp, ", %g", grid[i][j]);
      }
    }
    fprintf(fp, "]\n");
  }
}

int main(int argc, char* argv[]) {
  // take command line args
  gridSize = (argc > 1) ? ((atoi(argv[1]) < MAX_SIZE) ? atoi(argv[1]) : MAX_SIZE) : MAX_SIZE;
  numIters = (argc > 2) ? ((atoi(argv[2]) < MAX_ITERATIONS) ? atoi(argv[2]) : MAX_ITERATIONS) : MAX_ITERATIONS;

  // allocate memory for grids
  double** grid = malloc(gridSize*sizeof(double));
  double** new = malloc(gridSize*sizeof(double));
  for (int i = 0; i < gridSize; i++) {
    grid[i] = malloc(gridSize*sizeof(double));
    new[i] = malloc(gridSize*sizeof(double));
  }

  // initialize grids, boundary points are 1 the rest are 0
  for (int i = 0; i < gridSize; i++) {
    for (int j = 0; j < gridSize; j++) {
      if ((i == 0) || (j == 0) || (j == (gridSize-1)) || (i == (gridSize-1))) {
        grid[i][j] = 1;
        new[i][j] = 1;
      } else {
        grid[i][j] = 0;
        new[i][j] = 0;
      }
    }
  }
  printGrid(grid);
  return 0;
}