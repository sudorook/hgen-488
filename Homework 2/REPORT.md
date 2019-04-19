# Report

## Program description

This program will read in a sequence of amino acids from a text file and output
pairwise sequence alignments for all the provided sequences in a file called
\<input\>.align.txt.


### utils.h, utils.c

Implements utility functions for:

1. Reading/writing files

2. Finding the maximum value of three items

3. Implementing wrappers for matrix-like operations on arrays

   There is no matrix data structure in standard C, so matrices are simply left
   as arrays, and the get_value, set_value, use the size of the matrix to allow
   access to array locations as if they were indexed by rows and columns.

4. Printing out the score matrix. (For diagnostic purposes only.)


### globals.h

Contains constants and data types:

1. `enum aa`: Creates an enum for amino acids based on their 1-character codes.
   This makes the code more readable, because one can query the BLOSUM62 matrix
   with the amino acid character codes to get pairwise alignment scores.

2. `aa_table`: Contains amino acid codes as strings. Useful for converting
   alignments between enums and back to strings.

3. `struct sequence`: A struct for a sequence containing its string and enum
   representation, as well as its length.

4. `struct alignment`: A struct for storing the results of a pairwise
   alignment. Contains two aligned sequences and their score.

5. `blosum_array`: A copy-pasted static array for all the BLOSUM62 scores.


### main.c

Where the magic happens. Contains function used for actually computing the
alignment.

1. `get_blosum_score`: returns the BLOSUM62 score for a pair of amino acids.

2. `string_to_enum`: converts enum representation of amino acid sequence to a
   string.

3. `enum_to_string`: vice versa

4. `print_enum_seq`: print out the enum representation

5. `reverse_sequence`: used to reverse the traceback from the alignment score
   matrix. This is necessary because the traceback contains a reversed aligned
   sequence.

6. `traceback`: Takes a reference to a matrix of alignment scores and traverses
   the matrix to compute the aligned sequences and their score. It starts at
   the bottom right corner of the alignment matrix and moves up, left, or
   up-left based on the maximum score. An up- or left-movement signifies a gap
   in one sequence, and an up-left one signifies a sequence match. The
   traceback function continues to travel up and left until it reaches the
   upper left corner. Then the sequences and score are returned.

7. `global_sequence_alignment`: Computes a global sequence alignment for a pair
   of sequences of length m and n. It iterates down and left (row by column)
   from the upper right corner of a (m+1)x(n+1) matrix. The first row and
   column correspond to nothing but gaps, and the remaining cells are filled by
   taking the maximum of alignment scores for prior possible subsequences.
   Specifically, cell A[i,j] is the max(A[i-1,j-1] + M[i,j], A[i-1,j] + M[gap],
   A[i,j-1] + M[gap]), where M corresponds to the matrix of BLOSUM62 scores.


## Results

### Alignment

The sequence alignment for:

1: ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK
2: AKTTPPSVYPLAPGSAAQTNSMVTLGCLVKGYFPEPVTVTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVPSSTWPSQTVTCNVAHPASSTKVDKKIVPRDCGCKPCICTVPEVSSVFIFPPKPKDVLTITLTPKVTCVVVDISKDDPEVQFSWFVDDVEVHTAQTKPREEQFNSTFRSVSELPIMHQDWLNGKEFKCRVNSAAFPAPIEKTISKTKGRPKAPQVYTIPPPKEQMAKDKVSLTCMITDFFPEDITVEWQWNGQPAENYKNTQPIMDTDGSYFVYSKLNVQKSNWEAGNTFTCSVLHEGLHNHHTEKSLSHSPGK

yields:

1: ASTKGP_SVFPLAPS_S_K_STS_GGT_AALGCLVKDYFPEPVTVSWNSGALTS_GVHTFPAVLQSSGLYSLSS_VVTVPSSSLGTQTYIC_NV_NHKPS_NTKVDKKVEPKSCD_KTHTC__P__PCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHED_PEVKFNWYVDGVEVHNAKTKPREEQYNSTYR_VVS_VLTVLHQDWLNGKEYKCKVSNKAL_PAPIEKTISKAKGQPREPQVYTLPPSRDE_LTKNQVSLTCLVKGFYPSDIAVEWE_SNGQPE_NNYKT_TP_PVLDSDGSF_FLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK
2: AKTT_PPSVYPLAPGS_AAQT_NS_MVT__LGCLVKGYFPEPVTVTWNSGSLS_SGVHTFPAVLQSD_LYTLSSSV_TVPSSTWPSQTV_TCNVAH_P_ASSTKVDKKIVPRDCGCKP_C_ICTVP_____EV_SS__VFIFPPKPKDVLTITLTPKVTCVVVDISKD_DPEVQFSWFVDDVEVHTAQTKPREEQFNSTFRSV_SEL_PIMHQDWLNGKEFKCRVN_SAAFPAPIEKTISKTKGRPKAPQVYTIPPPKE_QMAKDKVSLTCMITDFFPEDITVEWQWN_GQPAEN_YKNT_QP_IMDTDGSYF_VYSKLNVQKSNWEAGNTFTCSVLHEGLHNHHTEKSLSHSPGK

with alignment score 1198.


### Comparison with BLASTp

The output score differs from BLASTp scores produced by the online tool because
of differences in how the scoring system was implemented. This alignment
program implements a flat -2 score for a gap, but BLASTp does not. It allows
users to choose among several penalty scores, which generally involve a high
penalty (-9, -12, etc.) for a gap initiation and a smaller penalty (-2, etc.)
for a gap extension.

For this program to recapitulate the online one, it would have to be modified
to keep track of gap lengths.

> As a side note, if it is possible for the global sequence alignment to yield
> the same score, this algorithm will not return all of them. It will simply
> return the first it finds.
