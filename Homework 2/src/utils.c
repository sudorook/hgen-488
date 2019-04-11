#include "utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Take a file name as input and create an output name replacing the file
 * extension with '.align.txt'
 */
char*
output_file(char* input)
{
  int len = strlen(input);

  int idx = 0;
  for (int i = len - 1; i != 0; i--) {
    char c = input[i];
    if (c == '.') {
      idx = i;
      break;
    }
  }

  char* output = malloc(sizeof(char) * (idx + 10)); // add extra characters
  strncpy(output, input, idx + 1);
  output[idx + 1] = 'a';
  output[idx + 2] = 'l';
  output[idx + 3] = 'i';
  output[idx + 4] = 'g';
  output[idx + 5] = 'n';
  output[idx + 6] = '.';
  output[idx + 7] = 't';
  output[idx + 8] = 'x';
  output[idx + 9] = 't';
  output[idx + 10] = '\0';

  return (output);
}

/* Implement a 3-way max function. */
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

/*
 * Wrapper function for matrix-like accessing array elements.
 */
int
get_value(int* matrix, int row, int col, int size)
{
  return (matrix[row * size + col]);
}

/*
 * Wrapper function for matrix-like assignment of values to array elements.
 */
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
