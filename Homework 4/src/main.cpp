#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#define CODON_TABLE_SIZE 64

static const std::string S1 = "GCTGCTGCTGCTAAACGTTTGGGGCAGTCGAT";
static const std::string S2 = "GGTGCTCCAAGCTTTTGAGTCTGCTAGTGTCAACCCT";
static const std::string S3 = "GTGGGCCCCCTAGCTAGCTAGCTGGGGCAC";
static const std::string S4 = "TGTCGCTGGCTGGACTGCTGATCGTAGTAG";

static const std::string S = "TAGCTAGCT";

struct global_hash_table
{
  std::vector<std::string> sequences;
  std::array<std::vector<std::tuple<int, int>>, CODON_TABLE_SIZE> table;
  int count = 0;
};

struct subsequence_match
{
  int id;
  std::string sequence;
  std::string match;
};

static std::map<char, int> dna = {
  { 'A', 0 },
  { 'C', 1 },
  { 'G', 2 },
  { 'T', 3 },
};

int
dna_hash(std::string triplet)
{
  if (triplet.size() != 3) {
    fprintf(stderr, "ERROR: non-triplet passed to dna_hash");
    exit(1);
  }

  int code =
    dna[triplet[0]] * (4 * 4) + dna[triplet[1]] * (4) + dna[triplet[2]];
  return (code);
}

void
print_global_hash_table(global_hash_table* hash_table)
{
  for (int i = 0; i != CODON_TABLE_SIZE; i++) {
    if (hash_table->table[i].size() > 0) {
      std::cout << i << ": ";
      for (auto j = hash_table->table[i].begin();
           j != hash_table->table[i].end();
           j++) {
        std::cout << "(" << std::get<0>(*j) << ',' << std::get<1>(*j) << ") ";
      }
      std::cout << std::endl;
    }
  }
  return;
}

void
append_global_hash_table(struct global_hash_table* hash_table, std::string seq)
{
  int len = seq.size();

  hash_table->count = hash_table->count + 1;
  hash_table->sequences.push_back(seq);

  for (int i = 0; i != len - 2; i++) {
    std::string substring = seq.substr(i, 3);
    hash_table->table[dna_hash(substring)].push_back(
      std::tuple<int, int>(hash_table->count, i));
  }
  return;
}

std::vector<struct subsequence_match>
substring_search(struct global_hash_table* hash_table, std::string substring)
{
  int len = substring.size();
  std::vector<std::vector<std::tuple<int, int, int>>> search_table;

  // Subset the global alignment by the motifs present in the sequence being
  // searched, and create a table of triplets (seq id, position of start of
  // sequence, position of start of motif).
  for (int i = 0; i != len - 2; i++) {
    std::string s = substring.substr(i, 3);
    std::vector<std::tuple<int, int, int>> v;
    for (auto j = hash_table->table[dna_hash(s)].begin();
         j != hash_table->table[dna_hash(s)].end();
         j++) {
      v.push_back(std::tuple<int, int, int>(
        std::get<0>(*j), (std::get<1>(*j)) - i, std::get<1>(*j)));
    }
    search_table.push_back(v);

    /* // Print the output
     * std::cout << i << ": ";
     * for (auto j = v.begin(); j != v.end(); j++) {
     *   std::cout << "(" << std::get<0>(*j) << ',' << std::get<1>(*j) << ','
     *             << std::get<2>(*j) << ") ";
     * }
     * std::cout << std::endl; */
  }

  std::vector<struct subsequence_match> matches;

  for (int i = 0; i < (int)search_table[0].size(); i++) {
    bool flag = true;

    int id = std::get<0>(search_table[0][i]);
    int start = std::get<1>(search_table[0][i]);

    for (int j = 1; j < len - 2; j++) {
      std::tuple<int, int, int> node;
      node = *std::find(search_table[j].begin(),
                        search_table[j].end(),
                        std::tuple<int, int, int>(id, start, start + j));
      if (node == (*search_table[j].end())) {
        flag = false;
        break;
      }
    }

    if (flag == true) {
      struct subsequence_match match;
      match.sequence = hash_table->sequences[id - 1];
      match.id = id;
      match.match =
        std::string(start, '.') + substring +
        std::string(match.sequence.size() - substring.size() - start, '.');
      matches.push_back(match);
    }
  }

  return (matches);
}

void
print_substring_matches(std::vector<struct subsequence_match> matches)
{
  for (auto match = matches.begin(); match != matches.end(); match++) {
    std::cout << "Sequence " << match->id << ": " << std::endl;
    std::cout << match->sequence << std::endl;
    std::cout << match->match << std::endl;
    std::cout << std::endl;
  }
}

int
main(void)
{
  struct global_hash_table table;

  append_global_hash_table(&table, S1);
  append_global_hash_table(&table, S2);
  append_global_hash_table(&table, S3);
  append_global_hash_table(&table, S4);
  // print_global_hash_table(&table);

  std::vector<struct subsequence_match> matches = substring_search(&table, S);
  print_substring_matches(matches);

  return (0);
}
