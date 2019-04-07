#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
/* #include <gsl/gsl_matrix.h>
 * #include <gsl/gsl_math.h> */

#define AA_TABLE_SIZE 24
#define NELEMS(x)  (sizeof(x) / sizeof((x)[0]))


/* 
 * Globals
 */

enum aa {
  A = 0,
  R,
  N,
  D,
  C,
  Q,
  E,
  G,
  H,
  I,
  L,
  K,
  M,
  F,
  P,
  S,
  T,
  W,
  Y,
  V,
  B,
  Z,
  X,
  _,
};

static char* aa_table[AA_TABLE_SIZE+1] = {
  "A",
  "R",
  "N",
  "D",
  "C",
  "Q",
  "E",
  "G",
  "H",
  "I",
  "L",
  "K",
  "M",
  "F",
  "P",
  "S",
  "T",
  "W",
  "Y",
  "V",
  "B",
  "Z",
  "X",
  "_",
  NULL
  /* '\0' */
};

struct alignment {
  char *sequence1;
  char *sequence2;
  int score;
};

/* int blosum_scores[] =  */

/*    A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
 * A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
 * R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
 * N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
 * D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
 * C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
 * Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
 * E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
 * G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
 * H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
 * I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
 * L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
 * K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
 * M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
 * F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
 * P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
 * S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
 * T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
 * W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
 * Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
 * V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
 * B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
 * Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
 * X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
 * * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 */
 
const int blosum_array[AA_TABLE_SIZE*AA_TABLE_SIZE] = { 
   4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -2, 
  -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -2, 
  -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -2, 
  -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -2, 
   0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -2, 
  -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -2, 
  -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -2, 
   0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -2, 
  -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -2, 
  -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -2, 
  -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -2, 
  -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -2, 
  -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -2, 
  -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -2, 
  -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -2, 
   1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -2, 
   0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -2, 
  -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -2, 
  -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -2, 
   0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -2, 
  -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -2, 
  -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -2, 
   0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -2, 
  -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  1
};

/* gsl_matrix *blosum_scores = gsl_matrix_alloc(AA_TABLE_SIZE, AA_TABLE_SIZE); */
/* const gsl_matrix_const_view blosum_scores = gsl_matrix_const_view_array( blosum_array, AA_TABLE_SIZE, AA_TABLE_SIZE ); */


/*
 * Functions
 */

int get_blosum_score(enum aa aa1, enum aa aa2) {
  return(blosum_array[(int)aa1 * AA_TABLE_SIZE + (int)aa2]);
}

enum aa* string_to_enum(char* seq) {
  int l = strlen(seq);
  enum aa *enum_seq = malloc(sizeof(enum aa)*l);
  for (int i=0; i < l; i++) {
    char aa_char = seq[i];
    for (int j=0; j < AA_TABLE_SIZE-1; j++) {
      if (*aa_table[j] == aa_char) {
        enum_seq[i] = j;
      }
    }
  }
  return(enum_seq);
}

char* enum_to_string(enum aa *seq, int l) {
  char* char_seq = malloc(sizeof(char)*l);
  for (int i=0; i < l; i++) {
    char_seq[i] = *aa_table[seq[i]];
  }
  return(char_seq);
}

void print_enum_seq(const enum aa* seq, int len) {
  for (int i=0; i < len; i++) {
    printf("%.2d ", seq[i]);
  }
  printf("\n"); 
}

void print_enum_seq2(const enum aa* seq, int len) {
  int i = 0;
  /* while(seq[i] != NULL) { */
  while(seq[i] != '\0') {
    printf("%.2d ", seq[i]);
    i++;
  }
  printf("\n"); 
}

int max(int a, int b, int c) {
  if (a > b) {
    if (a > c) {
      return(a);
    } else {
      return(c);
    }
  } else {
    if (b > c) {
      return(b);
    } else {
      return(c);
    }
  }
}

void set_value(int* matrix, int row, int col, int size, int val) {
  matrix[row*size + col] = val;
}

int get_value(int* matrix, int row, int col, int size) {
  return(matrix[row*size + col]);
}

enum aa* reverse_sequence(enum aa *seq, int len) {
  /* static enum aa newseq[len]; */
  enum aa *newseq = malloc(sizeof(enum aa)*len);
  for (int i=0; i < len; i++) {
    newseq[i] = seq[len-1-i];
  }
  return(newseq);
}

/* char* traceback(int *scores, const enum aa *seq1, const enum aa *seq2, int nrow, int ncol) { */
struct alignment traceback(int *scores, const enum aa *seq1, const enum aa *seq2, int nrow, int ncol) {
  int len = max(nrow, ncol, 0)*2; /* reuse the max function defined earlier */
  enum aa align_seq1[len];
  enum aa align_seq2[len];

  struct alignment trace;
  trace.score = get_value(scores, nrow-1, ncol-1, ncol);

  /* print_enum_seq(seq1, nrow-1);
   * print_enum_seq(seq2, ncol-1); */
  
  int i = nrow-1;
  int j = ncol-1;
  int idx = 0;

  while( (i + j) != 0 ) {
    /* printf("score: %d\n", get_value(scores, i, j, ncol)); */

    if ( (j != 0) && (i != 0) ) {
      int val1 = get_value(scores, i, j-1, ncol);
      int val2 = get_value(scores, i-1, j, ncol);
      int val3 = get_value(scores, i-1, j-1, ncol);

      if (val1 > val2) {
        if (val1 > val3) {
          /* printf("seq1 gap: (%d, %d)\n", _, seq2[j-1]); */
          align_seq1[idx] = _;
          align_seq2[idx] = seq2[j-1];
          idx++;
          j--;
        } else {
          /* printf("alignment: (%d, %d)\n", seq1[i-1], seq2[j-1]); */
          align_seq1[idx] = seq1[i-1];
          align_seq2[idx] = seq2[j-1];
          idx++;
          i--;
          j--;
        }
      } else {
        if (val2 > val3) {
          /* printf("seq2 gap: (%d, %d)\n", _, seq1[i-1]); */
          align_seq1[idx] = seq1[i-1];
          align_seq2[idx] = _;
          idx++;
          i--;
        } else {
          /* printf("alignment: (%d, %d)\n", seq1[i-1], seq2[j-1]); */
          align_seq1[idx] = seq1[i-1];
          align_seq2[idx] = seq2[j-1];
          idx++;
          i--;
          j--;
        }
      }
    } else if (j == 0) {
      align_seq1[idx] = seq1[i];
      align_seq2[idx] = _;
      idx++;
      i--;
    } else if (i == 0) {
      align_seq1[idx] = _;
      align_seq2[idx] = seq2[j];
      idx++;
      j--;
    }
  }
  align_seq1[idx] = '\0';
  align_seq2[idx] = '\0';

  /* print_enum_seq(align_seq1, idx);
   * print_enum_seq(align_seq2, idx); */
  
  /* printf("idx: %d\n", idx); */

  enum aa *newseq;
  newseq = reverse_sequence(align_seq1, idx);
  /* print_enum_seq(newseq, idx); */
  trace.sequence1 = enum_to_string(newseq, idx);
  /* printf("%s\n", trace.sequence1); */
  free(newseq);
  
  newseq = reverse_sequence(align_seq2, idx);
  /* print_enum_seq(newseq, idx); */
  trace.sequence2 = enum_to_string(newseq, idx);
  /* printf("%s\n", trace.sequence2); */
  free(newseq);

  return(trace);
}

/* void print_matrix(int* matrix, int nrow, int ncol, int size) { */
void print_matrix(int* matrix, int nrow, int ncol) {
  for (int i=0; i < nrow; i++) {
    for (int j=0; j < ncol; j++) {
      /* int val = get_value(matrix, i, j, size); */
      int val = get_value(matrix, i, j, ncol);
      printf("%3d ", val);
    }
    printf("\n");
  }
}

/* void zero_matrix(int* matrix, int nrow, int ncol, int size) { */
void zero_matrix(int* matrix, int nrow, int ncol) {
  for (int i=0; i < nrow; i++) {
    for (int j=0; j < ncol; j++) {
      /* set_value(matrix, i, j, size, 0); */
      set_value(matrix, i, j, ncol, 0);
    }
  }
}

struct alignment full_sequence_alignment(const char *seq1, const char *seq2) {
  const enum aa *enum_seq1;
  const int len_seq1 = strlen(seq1) + 1;
  enum_seq1 = string_to_enum(seq1);
  
  const enum aa *enum_seq2;
  const int len_seq2 = strlen(seq2) + 1;
  enum_seq2 = string_to_enum(seq2);

  int scores[(len_seq1)*(len_seq2)];

  /* Initialize the matrix */
  zero_matrix(scores, len_seq1, len_seq2);
  /* print_matrix(scores, len_seq1, len_seq2);
   * printf("\n"); */

  for (int i=0; i < len_seq1; i++) {
    for (int j=0; j < len_seq2; j++) {
      /* printf("%d\t%d\n", enum_seq1[i], enum_seq2[j]); */
      if (( i == 0) && (j == 0)) {
        set_value(scores, i, j, len_seq2, 0);
      } else if (i == 0) {
        /* printf("seq2: %d\t score: %d\n", enum_seq2[j-1], get_blosum_score(_, enum_seq2[j-1])); */
        int val = get_value(scores, i, j-1, len_seq2) + get_blosum_score(_, enum_seq2[j-1]);
        /* int val = get_value(scores, i, j-1, len_seq2) + -2; */
        set_value(scores, i, j, len_seq2, val);
      } else if (j == 0) {
        /* printf("seq1: %d\t score: %d\n", enum_seq1[i-1], get_blosum_score(enum_seq1[i-1], _)); */
        int val = get_value(scores, i-1, j, len_seq2) + get_blosum_score(enum_seq1[i-1], _);
        /* int val = get_value(scores, i-1, j, len_seq2) + -2; */
        set_value(scores, i, j, len_seq2, val);
      } else {
        int val1 = get_value(scores, i, j-1, len_seq2) + get_blosum_score(_, enum_seq2[j-1]);
        int val2 = get_value(scores, i-1, j, len_seq2) + get_blosum_score(enum_seq1[i-1], _);
        int val3 = get_value(scores, i-1, j-1, len_seq2) + get_blosum_score(enum_seq1[i-1], enum_seq2[j-1]);
        set_value(scores, i, j, len_seq2, max(val1, val2, val3));
      }
    }
  }
  print_matrix(scores, len_seq1, len_seq2);
  printf("\n");

  /* traceback(scores, enum_seq1, enum_seq1, len_seq1-1, len_seq2-1); */
  struct alignment trace = traceback(scores, enum_seq1, enum_seq2, len_seq1, len_seq2);
  return(trace);
}


/*
 * All happiness found here.
 */

int main(int argv, char* argc[]) {
  printf("Hi, how are you?\n");

  /* const char *sequence1 = "ARNDCQE"; */
  const char *sequence1 = "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK";
  /* const char *sequence2 = "ARRDQE"; */
  const char *sequence2 = "AKTTPPSVYPLAPGSAAQTNSMVTLGCLVKGYFPEPVTVTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVPSSTWPSQTVTCNVAHPASSTKVDKKIVPRDCGCKPCICTVPEVSSVFIFPPKPKDVLTITLTPKVTCVVVDISKDDPEVQFSWFVDDVEVHTAQTKPREEQFNSTFRSVSELPIMHQDWLNGKEFKCRVNSAAFPAPIEKTISKTKGRPKAPQVYTIPPPKEQMAKDKVSLTCMITDFFPEDITVEWQWNGQPAENYKNTQPIMDTDGSYFVYSKLNVQKSNWEAGNTFTCSVLHEGLHNHHTEKSLSHSPGK";
  
  const enum aa *enum_seq1;
  const int len_seq1 = strlen(sequence1);
  enum_seq1 = string_to_enum(sequence1);
  printf("%s\n", sequence1);
  print_enum_seq(enum_seq1, len_seq1);
  
  const enum aa *enum_seq2;
  const int len_seq2 = strlen(sequence2);
  enum_seq2 = string_to_enum(sequence2);
  printf("%s\n", sequence2);
  print_enum_seq(enum_seq2, len_seq2);

  struct alignment trace = full_sequence_alignment(sequence1, sequence2);

  printf("score: %d\n", trace.score);
  printf("%s\n", trace.sequence1);
  printf("%s\n", trace.sequence2);

  /* printf("%d\n", get_blosum_score(A, A));
   * printf("%d\n", get_blosum_score(A, R));
   * printf("%d\n", get_blosum_score(A, N));
   * printf("%d\n", get_blosum_score(A, D));
   * printf("%d\n", get_blosum_score(A, C));
   * printf("%d\n", get_blosum_score(A, Q));
   * printf("%d\n", get_blosum_score(A, E));
   * printf("%d\n", get_blosum_score(A, G));
   * printf("%d\n", get_blosum_score(A, H));
   * printf("%d\n", get_blosum_score(A, I));
   * printf("%d\n", get_blosum_score(A, L));
   * printf("%d\n", get_blosum_score(A, K));
   * printf("%d\n", get_blosum_score(A, M));
   * printf("%d\n", get_blosum_score(A, F));
   * printf("%d\n", get_blosum_score(A, P));
   * printf("%d\n", get_blosum_score(A, S));
   * printf("%d\n", get_blosum_score(A, T));
   * printf("%d\n", get_blosum_score(A, W));
   * printf("%d\n", get_blosum_score(A, Y));
   * printf("%d\n", get_blosum_score(A, V));
   * printf("%d\n", get_blosum_score(A, B));
   * printf("%d\n", get_blosum_score(A, Z));
   * printf("%d\n", get_blosum_score(A, X));
   * printf("%d\n", get_blosum_score(A, _)); */

  /* gsl_matrix *blosum_scores = gsl_matrix_alloc(AA_TABLE_SIZE, AA_TABLE_SIZE);
   * copy(blosum_array, AA_TABLE_SIZE*AA_TABLE_SIZE, blosum_scores->data); */
  /* const gsl_matrix_const_view blosum_scores = gsl_matrix_const_view_array( blosum_array, AA_TABLE_SIZE, AA_TABLE_SIZE ); */

  /* free(enum_seq1);
   * free(enum_seq2); */
  
  return(0);
}
