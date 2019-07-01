#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "globals.h"
#include "utils.h"

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
      if (aa_table[j] == aa_char) {
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
    char_seq[i] = aa_table[seq[i]];
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

/*
 * Utility function to reverse an array of AAs. This is needed because the
 * traceback function returns the optimal sequence in reverse order.
 */
enum aa*
reverse_sequence(enum aa* seq, int len)
{
  enum aa* newseq = malloc(sizeof(enum aa) * len);
  for (int i = 0; i < len; i++) {
    newseq[i] = seq[len - 1 - i];
  }
  return (newseq);
}

/*
 * Compute a traceback to return the optimal sequence alignment for a given
 * scoring matrix. Returns the score and the matches sequences (gaps included).
 */
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

/*
 * Populate a matrix of best alignments.
 */
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
  /* print_matrix(scores, nrow, ncol);
   * printf("\n"); */

  return (traceback(scores, seq1, seq2));
}

void
print_alignment(struct alignment a)
{
  printf("score: %d\n", a.score);
  printf("seq 1: %s\n", a.sequence1.char_seq);
  printf("seq 2: %s\n", a.sequence2.char_seq);
  printf("\n");
  return;
}

void
write_alignment(struct alignment a, FILE* f)
{
  fprintf(f, "score: %d\n", a.score);
  fprintf(f, "seq 1: %s\n", a.sequence1.char_seq);
  fprintf(f, "seq 2: %s\n", a.sequence2.char_seq);
  fprintf(f, "\n");
  return;
}

/*
 * All happiness found here.
 */

int
main(int argv, char* argc[])
{
  /* if (argv != 3) { */
  if (argv != 2) {
    fprintf(stderr,
            "Specify one input file containing protein sequences to align.\n");
    return (1);
  }
  printf("Loading sequences in the \'%s\' file.\n", argc[argv - 1]);

  char* input = argc[1];
  char* output = output_file(input);

  FILE* f = fopen(input, "r");
  if (f == NULL) {
    fprintf(stderr, "Error opening file \'%s.\'\n", input);
    exit(2);
  }

  /* Read the sequences provided in the input file. */
  int i = 0;
  int count = 0;
  char c;
  char* sequences[MAX_N_SEQUENCES];   // compare at most 1000 sequences.
  char sequence[MAX_SEQUENCE_LENGTH]; // million sequence limit.
  do {
    c = fgetc(f);
    if (feof(f)) {
      break;
    } else if (c == '\n') {
      /* Start reading a new sequences when a newline is reached. */
      sequence[i] = '\0';
      sequences[count] = malloc(strlen(sequence) + 1);
      strcpy(sequences[count], sequence);
      count++;
      i = 0;
    } else {
      sequence[i] = c;
      i++;
    }
  } while (1);

  fclose(f);

  /* Print the input sequences. */
  if (count > 9) {
    for (int i = 0; i < count; i++) {
      printf("seq %2d: %s\n", i + 1, sequences[i]);
    }
  } else {
    for (int i = 0; i < count; i++) {
      printf("seq %d: %s\n", i + 1, sequences[i]);
    }
  }
  printf("\n");

  /* Dynamically allocate memory for all the sequences and store the input
   * sequences as sequences tructs (members: length, char*, and enum aa*).
   */
  struct sequence* seqs = malloc(sizeof(sequence) * count);
  for (int i = 0; i < count; i++) {
    const struct sequence seq = { .char_seq = sequences[i],
                                  .enum_seq = string_to_enum(sequences[i]),
                                  .length = strlen(sequences[i]) };
    seqs[i] = seq;
  }

  /* Run pairwise sequence alignments for all sequences. */
  f = fopen(output, "w");
  if (f == NULL) {
    fprintf(stderr, "Error opening file \'%s.\'\n", output);
    exit(2);
  }

  for (int i = 0; i < count; i++) {
    for (int j = i + 1; j < count; j++) {
      struct alignment trace = global_sequence_alignment(seqs[i], seqs[j]);
      print_alignment(trace);
      write_alignment(trace, f);

      /* Free memory allocated on the heap. */
      free(trace.sequence1.char_seq);
      free(trace.sequence2.char_seq);
      free(trace.sequence1.enum_seq);
      free(trace.sequence2.enum_seq);
    }
  }

  fclose(f);

  /* Cleanup and free memory. */
  for (int i = 0; i < count; i++) {
    free(sequences[i]);
    free(seqs[i].enum_seq);
  }

  free(seqs);
  free(output);

  return (0);
}
