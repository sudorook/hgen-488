# Homework 2

## Requirements

Requires standard C libraries.


## Usage

**Input:** text file (ASCII) containing sequences and nothing more.

**Output:** a text file named \<input\>.align.txt, containing pairwise
alignments and scores for all sequences in the input file. Format:

```
score
seq 1 alignment
seq 2 alignment

score
seq 1 alignment
seq 2 alignment

.
.
.
```

## Compilation

Run:
```
make
./align_blosum62 <input>
```
