#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "globals.h"

/* Sequence */

/* static char *aa_sequence1 = "ARNDCQE\0";
 * static char *aa_sequence2 = "ARRDQE\0"; */

/* static char *aa_sequence1 = "ARNDCQEAASDGQ\0";
 * static char *aa_sequence2 = "ARRDQEAWERWD\0";  */

static char* aa_sequence1 =
  "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSL"
  "GTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHE"
  "DPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPRE"
  "PQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGN"
  "VFSCSVMHEALHNHYTQKSLSLSPGK";
static char* aa_sequence2 =
  "AKTTPPSVYPLAPGSAAQTNSMVTLGCLVKGYFPEPVTVTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVPSSTWP"
  "SQTVTCNVAHPASSTKVDKKIVPRDCGCKPCICTVPEVSSVFIFPPKPKDVLTITLTPKVTCVVVDISKDDPEVQF"
  "SWFVDDVEVHTAQTKPREEQFNSTFRSVSELPIMHQDWLNGKEFKCRVNSAAFPAPIEKTISKTKGRPKAPQVYTI"
  "PPPKEQMAKDKVSLTCMITDFFPEDITVEWQWNGQPAENYKNTQPIMDTDGSYFVYSKLNVQKSNWEAGNTFTCSV"
  "LHEGLHNHHTEKSLSHSPGK";

/*
 * Functions
 */

int
get_blosum_score(enum aa aa1, enum aa aa2)
{
  return (blosum_array[(int)aa1 * AA_TABLE_SIZE + (int)aa2]);
}

enum aa*
string_to_enum(const char* seq)
{
  int l = strlen(seq);
  enum aa* enum_seq = malloc(sizeof(enum aa) * l);
  for (int i = 0; i < l; i++) {
    char aa_char = seq[i];
    for (int j = 0; j < AA_TABLE_SIZE - 1; j++) {
      if (*aa_table[j] == aa_char) {
        enum_seq[i] = j;
      }
    }
  }
  return (enum_seq);
}

char*
enum_to_string(enum aa* seq, int l)
{
  char* char_seq = malloc(sizeof(char) * l);
  for (int i = 0; i < l; i++) {
    char_seq[i] = *aa_table[seq[i]];
  }
  return (char_seq);
}

void
print_enum_seq(const enum aa* seq, int len)
{
  for (int i = 0; i < len; i++) {
    printf("%.2d ", seq[i]);
  }
  printf("\n");
}

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

void
set_value(int* matrix, int row, int col, int size, int val)
{
  matrix[row * size + col] = val;
}

int
get_value(int* matrix, int row, int col, int size)
{
  return (matrix[row * size + col]);
}

enum aa*
reverse_sequence(enum aa* seq, int len)
{
  enum aa* newseq = malloc(sizeof(enum aa) * len);
  for (int i = 0; i < len; i++) {
    newseq[i] = seq[len - 1 - i];
  }
  return (newseq);
}

struct alignment
traceback(int* scores, const struct sequence seq1, const struct sequence seq2)
{
  int nrow = seq1.length + 1;
  int ncol = seq2.length + 1;
  int size = seq2.length + 1;

  int len = max(nrow, ncol, 0) * 2; /* reuse the max function defined earlier */
  enum aa align_seq1[len];
  enum aa align_seq2[len];

  int i = nrow - 1;
  int j = ncol - 1;
  int idx = 0;

  while ((i + j) != 0) {
    if ((j != 0) && (i != 0)) {
      int val1 = get_value(scores, i, j - 1, size);
      int val2 = get_value(scores, i - 1, j, size);
      int val3 = get_value(scores, i - 1, j - 1, size);

      if (val1 > val2) {
        if (val1 > val3) {
          align_seq1[idx] = _;
          align_seq2[idx] = seq2.enum_seq[j - 1];
          idx++;
          j--;
        } else {
          align_seq1[idx] = seq1.enum_seq[i - 1];
          align_seq2[idx] = seq2.enum_seq[j - 1];
          idx++;
          i--;
          j--;
        }
      } else {
        if (val2 > val3) {
          align_seq1[idx] = seq1.enum_seq[i - 1];
          align_seq2[idx] = _;
          idx++;
          i--;
        } else {
          align_seq1[idx] = seq1.enum_seq[i - 1];
          align_seq2[idx] = seq2.enum_seq[j - 1];
          idx++;
          i--;
          j--;
        }
      }
    } else if (j == 0) {
      align_seq1[idx] = seq1.enum_seq[i];
      align_seq2[idx] = _;
      idx++;
      i--;
    } else if (i == 0) {
      align_seq1[idx] = _;
      align_seq2[idx] = seq2.enum_seq[j];
      idx++;
      j--;
    }
  }
  align_seq1[idx] = '\0';
  align_seq2[idx] = '\0';

  struct alignment trace;
  trace.score = get_value(scores, nrow - 1, ncol - 1, size);

  enum aa* newseq;
  newseq = reverse_sequence(align_seq1, idx);
  struct sequence new_seq1 = { .enum_seq = newseq,
                               .char_seq = enum_to_string(newseq, idx),
                               .length = idx };
  trace.sequence1 = new_seq1;

  newseq = reverse_sequence(align_seq2, idx);
  struct sequence new_seq2 = { .enum_seq = newseq,
                               .char_seq = enum_to_string(newseq, idx),
                               .length = idx };
  trace.sequence2 = new_seq2;

  return (trace);
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

struct alignment
global_sequence_alignment(const struct sequence seq1,
                          const struct sequence seq2)
{

  int nrow = seq1.length + 1;
  int ncol = seq2.length + 1;
  int size = seq2.length + 1;
  int scores[nrow * ncol];

  /* Initialize the matrix */
  zero_matrix(scores, nrow, ncol);

  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      if ((i == 0) && (j == 0)) {
        set_value(scores, i, j, size, 0);
      } else if (i == 0) {
        int val = get_value(scores, i, j - 1, size) +
                  get_blosum_score(_, seq2.enum_seq[j - 1]);
        set_value(scores, i, j, size, val);
      } else if (j == 0) {
        int val = get_value(scores, i - 1, j, size) +
                  get_blosum_score(seq1.enum_seq[i - 1], _);
        set_value(scores, i, j, size, val);
      } else {
        int val1 = get_value(scores, i, j - 1, size) +
                   get_blosum_score(_, seq2.enum_seq[j - 1]);
        int val2 = get_value(scores, i - 1, j, size) +
                   get_blosum_score(seq1.enum_seq[i - 1], _);
        int val3 = get_value(scores, i - 1, j - 1, size) +
                   get_blosum_score(seq1.enum_seq[i - 1], seq2.enum_seq[j - 1]);
        set_value(scores, i, j, size, max(val1, val2, val3));
      }
    }
  }
  print_matrix(scores, nrow, ncol);
  printf("\n");

  return (traceback(scores, seq1, seq2));
}

/*
 * All happiness found here.
 */

int
main(int argv, char* argc[])
{

  const struct sequence seq1 = { .char_seq = aa_sequence1,
                                 .enum_seq = string_to_enum(aa_sequence1),
                                 .length = strlen(aa_sequence1) };

  const struct sequence seq2 = { .char_seq = aa_sequence2,
                                 .enum_seq = string_to_enum(aa_sequence2),
                                 .length = strlen(aa_sequence2) };

  printf("seq 1: %s\n", seq1.char_seq);
  printf("seq 2: %s\n", seq2.char_seq);
  print_enum_seq(seq1.enum_seq, seq1.length);
  print_enum_seq(seq2.enum_seq, seq2.length);

  struct alignment trace = global_sequence_alignment(seq1, seq2);

  printf("score: %d\n", trace.score);
  printf("seq 1: %s\n", trace.sequence1.char_seq);
  printf("seq 2: %s\n", trace.sequence2.char_seq);

  return (0);
}
