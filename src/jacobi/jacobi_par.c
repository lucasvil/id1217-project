/**
 * A parallel Jacobi iteration program for solving partial differential equations.
 *
 * Usage: ./a.out <gridSize> <numIters> <numWorkers>
 * Output:  Prints command-line arguments used, execution time and maximum error to stdout.
 *          Prints final grid values to a file 'seq-data.out'
 *
 * Authors: Lucas Villarroel, Erik Hanstad
 * Language: C using the OpenMP API
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define MAX_ITERATIONS 100000000
#define MAX_SIZE 10000
#define MAX_WORKERS 4

int gridSize, numIters, numWorkers;

double maxDiff(double** grid, double** new) {
  int i, j;
  double maxDiff = 0.0;
  double myDiff = -1;
#pragma omp parallel for private(j, myDiff)
  for (i = 1; i < gridSize; i++) {
    for (j = 1; j < gridSize; j++) {
      myDiff = fabs(grid[i][j] - new[i][j]);
      if (myDiff > maxDiff) {
#pragma omp critical
        {
          if (myDiff > maxDiff)
            maxDiff = myDiff;
        }
      }
    }
  }
  return maxDiff;
}

void printGrid(double** grid) {
  FILE* fp = fopen("par-data.out", "w");
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

void jacobi(double** grid, double** new) {
  int i, j;
  int boundary = gridSize - 1;
  for (int iter = 1; iter < numIters; iter++) {
#pragma omp parallel for private(j)
    for (i = 1; i < boundary; i++) {
      for (j = 1; j < boundary; j++) {
        new[i][j] = (grid[i - 1][j] + grid[i + 1][j] + grid[i][j - 1] + grid[i][j + 1]) * 0.25;
      }
    }
#pragma omp parallel for private(j)
    for (i = 1; i < boundary; i++) {
      for (j = 1; j < boundary; j++) {
        grid[i][j] = (new[i - 1][j] + new[i + 1][j] + new[i][j - 1] + new[i][j + 1]) * 0.25;
      }
    }
#ifdef debug
    printf("iter %d/%d\n", iter, numIters);
#endif // debug
  }
  maxDiff(grid, new);
}

int main(int argc, char* argv[]) {
  // take command line args
  gridSize = (argc > 1) ? ((atoi(argv[1]) < MAX_SIZE) ? atoi(argv[1]) : MAX_SIZE) : MAX_SIZE;
  numIters = (argc > 2) ? ((atoi(argv[2]) < MAX_ITERATIONS) ? atoi(argv[2]) : MAX_ITERATIONS) : MAX_ITERATIONS;
  numWorkers = (argc > 3) ? ((atoi(argv[3]) < MAX_WORKERS) ? atoi(argv[3]) : MAX_WORKERS) : MAX_WORKERS;

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

  omp_set_num_threads(numWorkers);

  double start_time = omp_get_wtime();
  jacobi(grid, new);
  double diff = maxDiff(grid, new);
  double end_time = omp_get_wtime();

  int ms = (end_time - start_time) * 1000;
  printf("Size: %d, Iterations: %d\n", gridSize-2, numIters);
  printf("Execution time: %d ms\n", ms);
  printf("Maximum error %g\n", diff);
  printGrid(grid);
  return 0;
}