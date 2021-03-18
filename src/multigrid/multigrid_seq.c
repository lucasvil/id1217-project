/**
 * A sequential multigrid program for solving partial differential equations.
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
#define FEW_ITERATIONS 4

double maxDiff(double** grid, double** new, int size) {
  double maxDiff = 0.0;
  for (int i = 1; i <= size; i++) {
    for (int j = 1; j <= size; j++) {
      maxDiff = (fabs(grid[i][j] - new[i][j]) > maxDiff) ? fabs(grid[i][j] - new[i][j]) : maxDiff;
    }
  }
  return maxDiff;
}

void printGrid(double** grid, int size) {
  FILE* fp = fopen("seq-data.out", "w");
  for (int i = 0; i < size; i++) {
    fprintf(fp, "[");
    for (int j = 0; j < size; j++) {
      if (j == 0) {
        fprintf(fp, "%.10f", grid[i][j]);
      }
      else {
        fprintf(fp, ", %.10f", grid[i][j]);
      }
    }
    fprintf(fp, "]\n");
  }
}

void jacobi(double** grid, double** new, int size, int numIters) {
  for (int iter = 1; iter < numIters; iter++) {
    for (int i = 1; i <= size; i++) {
      for (int j = 1; j <= size; j++) {
        new[i][j] = (grid[i - 1][j] + grid[i + 1][j] + grid[i][j - 1] + grid[i][j + 1]) * 0.25;
      }
    }
    for (int i = 1; i <= size; i++) {
      for (int j = 1; j <= size; j++) {
        grid[i][j] = (new[i - 1][j] + new[i + 1][j] + new[i][j - 1] + new[i][j + 1]) * 0.25;
      }
    }
#ifdef debug
    printf("iter %d/%d\n", iter, numIters);
#endif // debug
  }
}

void restrictGrid(double** fine, double** coarse, int coarseSize) {
  int x, y, i, j;
  for (i = 1; i <= coarseSize; i++) {
    x = i * 2;
    for (j = 1; j <= coarseSize; j++) {
      y = j * 2;
      coarse[i][j] = fine[x][y] * 0.5 + (fine[x - 1][y] + fine[x][y - 1] + fine[x][y + 1] + fine[x + 1][y]) * 0.125;
    }
  }
}

void interpolate(double** fine, double** coarse, int fineSize, int coarseSize) {
  int x, y;
  // assign coarse grid points to corresponding fine grid points
  for (int i = 1; i <= coarseSize; i++) {
    x = i * 2;
    for (int j = 1; j <= coarseSize; j++) {
      y = j * 2;
      fine[x][y] = coarse[i][j];
    }
  }


  // update fine grid points in columns that were updated
  for (int i = 2; i <= fineSize; i += 2) {
    for (int j = 1; j <= fineSize; j += 2) {
      fine[i][j] = (fine[i - 1][j] + fine[i + 1][j]) * 0.5;
    }
  }

  // update other fine grid points
  for (int i = 1; i <= fineSize; i++) {
    for (int j = 2; j <= fineSize; j += 2) {
      fine[i][j] = (fine[i][j - 1] + fine[i][j + 1]) * 0.5;
    }
  }
}

double** initGrid(int size) {
  double** grid = malloc(size * sizeof(double));
  for (int i = 0; i < size; i++) {
    grid[i] = malloc(size * sizeof(double));
  }

  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      if ((i == 0) || (j == 0) || (j == (size - 1)) || (i == (size - 1))) {
        grid[i][j] = 1;
      }
      else {
        grid[i][j] = 0;
      }
    }
  }
  return grid;
}

int main(int argc, char* argv[]) {
  // take command line args
  int size1 = (argc > 1) ? ((atoi(argv[1]) < MAX_SIZE) ? atoi(argv[1]) : MAX_SIZE) : MAX_SIZE;
  int numIters = (argc > 2) ? ((atoi(argv[2]) < MAX_ITERATIONS) ? atoi(argv[2]) : MAX_ITERATIONS) : MAX_ITERATIONS;

  int size2 = (size1 * 2) + 1;
  int size3 = (size2 * 2) + 1;
  int size4 = (size3 * 2) + 1;

  // allocate memory for grids
  double** grid1 = initGrid(size1 + 2);
  double** new1 = initGrid(size1 + 2);
  double** grid2 = initGrid(size2 + 2);
  double** new2 = initGrid(size2 + 2);
  double** grid3 = initGrid(size3 + 2);
  double** new3 = initGrid(size3 + 2);
  double** grid4 = initGrid(size4 + 2);
  double** new4 = initGrid(size4 + 2);

  // start timer
  clock_t start = clock();

  // restrict from finest grid down
  jacobi(grid4, new4, size4, FEW_ITERATIONS);
  restrictGrid(grid4, grid3, size3);

  jacobi(grid3, new3, size3, FEW_ITERATIONS);
  restrictGrid(grid3, grid2, size2);

  jacobi(grid2, new2, size2, FEW_ITERATIONS);
  restrictGrid(grid2, grid1, size1);

  // coarsest grid reached, compute solution
  jacobi(grid1, new1, size1, numIters);


  // interpolate back to finest grid
  interpolate(grid2, grid1, size2, size1);
  jacobi(grid2, new2, size2, FEW_ITERATIONS);
  interpolate(grid3, grid2, size3, size2);
  jacobi(grid3, new3, size3, FEW_ITERATIONS);
  interpolate(grid4, grid3, size4, size3);
  jacobi(grid4, new4, size4, FEW_ITERATIONS);

  // calculate max difference
  double diff = maxDiff(grid4, new4, size4);

  // stop timer
  clock_t end = clock();
  int msec = ((end - start) * 1000) / CLOCKS_PER_SEC;

  printf("Size: %d, Iterations: %d\n", size1, numIters);
  printf("Execution time: %d ms\n", msec);
  printf("Maximum error %.5f\n", diff);

  printGrid(grid4, size4 + 2);

  return 0;
}