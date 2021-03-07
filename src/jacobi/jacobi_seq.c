/**
 * A sequential Jacobi iteration program for solving partial differential equations.
 *
 * Usage: ./a.out <size> <iterations>
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

int size, iterations;
double** grid;
double** new;

void printGrid(double** grid) {
  FILE* fp = fopen("seq-data.out", "w");
  for (int i = 0; i < size; i++) {
    fprintf(fp, "[");
    for (int j = 0; j < size; j++) {
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
  size = (argc > 1) ? ((atoi(argv[1]) < MAX_SIZE) ? atoi(argv[1]) : MAX_SIZE) : MAX_SIZE;
  iterations = (argc > 2) ? ((atoi(argv[2]) < MAX_ITERATIONS) ? atoi(argv[2]) : MAX_ITERATIONS) : MAX_ITERATIONS;

  // allocate memory for grids
  grid = malloc(size*sizeof(double));
  new = malloc(size*sizeof(double));
  for (int i = 0; i < size; i++) {
    grid[i] = malloc(size*sizeof(double));
    new[i] = malloc(size*sizeof(double));
  }

  // initialize grids, boundary points are 1 the rest are 0
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      if ((i == 0) || (j == 0) || (j == (size-1)) || (i == (size-1))) {
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