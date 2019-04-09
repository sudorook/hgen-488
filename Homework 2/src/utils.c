#include "utils.h"

#include <stdio.h>

int
max(int a, int b, int c)
{
  if (a > b) {
    if (a > c) {
      return (a);
    } else {
      return (c);
    }
  } else {
    if (b > c) {
      return (b);
    } else {
      return (c);
    }
  }
}

int
get_value(int* matrix, int row, int col, int size)
{
  return (matrix[row * size + col]);
}

void
set_value(int* matrix, int row, int col, int size, int val)
{
  matrix[row * size + col] = val;
}

/* Print out the contents of a matrix. */
void
print_matrix(int* matrix, int nrow, int ncol)
{
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      int val = get_value(matrix, i, j, ncol);
      printf("%3d ", val);
    }
    printf("\n");
  }
}

/* Zero out the entries of a matrix. */
void
zero_matrix(int* matrix, int nrow, int ncol)
{
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      set_value(matrix, i, j, ncol, 0);
    }
  }
}
