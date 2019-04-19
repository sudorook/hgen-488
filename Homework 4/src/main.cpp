#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <array>
#include <map>

#define CODON_TABLE_SIZE 64

static const std::string S1 = "GCTGCTGCTGCTAAACGTTTGGGGCAGTCGAT";
static const std::string S2 = "GGTGCTCCAAGCTTTTGAGTCTGCTAGTGTCAACCCT";
static const std::string S3 = "GTGGGCCCCCTAGCTAGCTAGCTGGGGCAC";
static const std::string S4 = "TGTCGCTGGCTGGACTGCTGATCGTAGTAG";

static const std::string S = "TAGCTAGCT";

struct global_hash_table {
  std::vector<std::string> sequences;
  std::array<std::vector<std::tuple<int,int>>, CODON_TABLE_SIZE> table;
  int count=0;
};

static std::map<char, int> dna = {
  { 'A', 0 },
  { 'C', 1 },
  { 'G', 2 },
  { 'T', 3 },
};

int dna_hash(std::string triplet) {
  if (triplet.size() != 3) {
    fprintf(stderr, "ERROR: non-triplet passed to dna_hash");
    exit(1);
  }

  int code = dna[triplet[0]]*(4*4) + dna[triplet[1]]*(4) + dna[triplet[2]];
  return(code);
}

void print_global_hash_table(global_hash_table *table) {
  for (int i=0; i!=CODON_TABLE_SIZE; i++) {
    if (table->table[i].size() > 0) {
      std::cout << i << ": ";
      for (auto j = table->table[i].begin(); j != table->table[i].end(); j++) {
        std::cout << "(" << std::get<0>(*j) << ',' << std::get<1>(*j) << ") ";
      }
      std::cout<<std::endl;
    }
  }
  return;
}

void print_global_hash_table(std::array<std::vector<std::tuple<int,int>>, CODON_TABLE_SIZE> *table) {
  for (int i=0; i!=CODON_TABLE_SIZE; i++) {
    if ((*table)[i].size() > 0) {
      std::cout << i << ": ";
      for (auto j = (*table)[i].begin(); j != (*table)[i].end(); j++) {
        std::cout << "(" << std::get<0>(*j) << ',' << std::get<1>(*j) << ") ";
      }
      std::cout<<std::endl;
    }
  }
  return;
}

void append_global_hash_table(struct global_hash_table *table, std::string seq) {
  int len = seq.size();

  table->count = table->count+1;
  table->sequences.push_back(seq);

  for (int i=0; i!=len-2; i++) {
    std::string substring = seq.substr(i, 3);
    table->table[dna_hash(substring)].push_back(std::tuple<int,int>(table->count,i));
  }
  return;
}

void substring_search(struct global_hash_table *table, std::string substring) {
  int len = substring.size();
  std::vector<std::vector<std::tuple<int,int,int>>> search_table;
  
  // Subset the global alignment by the motifs present in the sequence being
  // searched, and create a table of triplets (seq id, position of start of
  // sequence, position of start of motif).
  for (int i=0; i!=len-2; i++) {
    std::string s = substring.substr(i, 3);
    std::vector<std::tuple<int,int,int>> v;
    for (auto j = table->table[dna_hash(s)].begin(); j != table->table[dna_hash(s)].end(); j++) {
      v.push_back(std::tuple<int,int,int>(std::get<0>(*j), (std::get<1>(*j))-i, std::get<1>(*j)));
    }
    search_table.push_back(v);

    // print the output
    std::cout << i << ": ";
    for (auto j = v.begin(); j != v.end(); j++) {
      std::cout << "(" << std::get<0>(*j) << ','<< std::get<1>(*j) << ',' << std::get<2>(*j) << ") ";
    }
    std::cout<<std::endl;
  }

  std::vector<std::tuple<int,int>> indices;

  for (int i=0; i<(int)search_table[0].size(); i++) {
    bool flag = true;

    int id = std::get<0>(search_table[0][i]);
    int start = std::get<1>(search_table[0][i]);

    for (int j=1; j<len-2; j++) {
      std::tuple<int,int,int> node;
      node = *std::find(search_table[j].begin(),
                        search_table[j].end(),
                        std::tuple<int,int,int>(id,start,start+j));
      if ( node == (*search_table[j].end()) ) {
        flag = false;
        break;
      }
    }

    if (flag == true) {
      indices.push_back(std::tuple<int,int>(id,start));
    }
  }

  for (auto j = indices.begin(); j != indices.end(); j++) {
    int id = std::get<0>(*j);
    int start = std::get<1>(*j);
    std::string sequence = table->sequences[std::get<0>(*j)-1];

    std::cout<<"Sequence "<<id<<": "<<std::endl;
    std::cout<<sequence<<std::endl;
    std::cout<<std::string(start,'.');
    std::cout<<substring;
    std::cout<<std::string(sequence.size() - substring.size() - start,'.')<<std::endl;
  }

  return;
}

int main(void) {
  struct global_hash_table table;

  append_global_hash_table(&table, S1);
  append_global_hash_table(&table, S2);
  append_global_hash_table(&table, S3);
  append_global_hash_table(&table, S4);
  print_global_hash_table(&table);

  substring_search(&table, S);

  return(0);
}
