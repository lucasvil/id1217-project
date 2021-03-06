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
#include <math.h>
#include <time.h>

#define MAX_ITERATIONS 100000000
#define MAX_SIZE 10000

int gridSize, numIters;

double maxDiff(double** grid, double** new) {
  double maxDiff = 0.0;
  for (int i = 1; i < gridSize; i++) {
    for (int j = 1; j < gridSize; j++) {
      maxDiff = (fabs(grid[i][j] - new[i][j]) > maxDiff) ? fabs(grid[i][j] - new[i][j]) : maxDiff;
    }
  }
  maxDiff;
}

void jacobi(double** grid, double** new) {
  int boundary = gridSize - 1;
  for (int iter = 1; iter < numIters; iter++) {
    for (int i = 1; i < boundary; i++) {
      for (int j = 1; j < boundary; j++) {
        new[i][j] = (grid[i - 1][j] + grid[i + 1][j] + grid[i][j - 1] + grid[i][j + 1]) * 0.25;
      }
    }
    for (int i = 1; i < boundary; i++) {
      for (int j = 1; j < boundary; j++) {
        grid[i][j] = (new[i - 1][j] + new[i + 1][j] + new[i][j - 1] + new[i][j + 1]) * 0.25;
      }
    }
#ifdef debug
    printf("iter %d/%d\n", iter, numIters);
#endif // debug
  }
}

void printGrid(double** grid) {
  FILE* fp = fopen("seq-data.out", "w");
  for (int i = 0; i < gridSize; i++) {
    fprintf(fp, "[");
    for (int j = 0; j < gridSize; j++) {
      if (j == 0) {
        fprintf(fp, "%.10f", grid[i][j]);
      } else {
        fprintf(fp, ", %.10f", grid[i][j]);
      }
    }
    fprintf(fp, "]\n");
  }
}

int main(int argc, char* argv[]) {
  // take command line args
  gridSize = (argc > 1) ? ((atoi(argv[1]) < MAX_SIZE) ? atoi(argv[1]) : MAX_SIZE) : MAX_SIZE;
  numIters = (argc > 2) ? ((atoi(argv[2]) < MAX_ITERATIONS) ? atoi(argv[2]) : MAX_ITERATIONS) : MAX_ITERATIONS;

  gridSize += 2;
  // allocate memory for grids
  double** grid = malloc(gridSize * sizeof(double));
  double** new = malloc(gridSize * sizeof(double));
  for (int i = 0; i < gridSize; i++) {
    grid[i] = malloc(gridSize * sizeof(double));
    new[i] = malloc(gridSize * sizeof(double));
  }

  // initialize grids, boundary points are 1 the rest are 0
  for (int i = 0; i < gridSize; i++) {
    for (int j = 0; j < gridSize; j++) {
      if ((i == 0) || (j == 0) || (j == (gridSize - 1)) || (i == (gridSize - 1))) {
        grid[i][j] = 1;
        new[i][j] = 1;
      } else {
        grid[i][j] = 0;
        new[i][j] = 0;
      }
    }
  }

  clock_t start = clock();
  jacobi(grid, new);
  double diff = maxDiff(grid, new);
  clock_t end = clock();
  int msec = ((end - start) * 1000) / CLOCKS_PER_SEC;
  printf("Size: %d, Iterations: %d\n", gridSize-2, numIters);
  printf("Execution time: %d ms\n", msec);
  printf("Maximum error %g\n", diff);
  printGrid(grid);
  return 0;
}