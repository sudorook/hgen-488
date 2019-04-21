#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

// static const std::string string = "ABAABA";
// static const std::string string = "BANANA";
static const std::string string = "DOGWOOD";
// static const std::string string = "DOGWORD";

struct LF_map {
  std::vector<std::tuple<char, int>> F;
  std::vector<std::tuple<char, int>> L;
};

struct LF_map make_LF_map(std::string bwt_string) {
  std::string sorted_string = bwt_string;
  std::sort(sorted_string.begin(), sorted_string.end());
  
  std::vector<std::tuple<char, int>> sorted;
  for (int i = 0; i < (int)sorted_string.size(); i++) {
    int counter=0;
    for (int j=0; j < i; j++) {
      if (sorted_string[i] == sorted_string[j]) {
        counter++;
      }
    }
    sorted.push_back(std::tuple<char, int>(sorted_string[i], counter));
  }

  std::vector<std::tuple<char, int>> bwt;
  for (int i = 0; i<(int)bwt_string.size(); i++) {
    int counter=0;
    for (int j=0; j < i; j++) {
      if (bwt_string[i] == bwt_string[j]) {
        counter++;
      }
    }
    bwt.push_back(std::tuple<char, int>(bwt_string[i], counter));
  }

  struct LF_map map;
  map.F = sorted;
  map.L = bwt;

  return(map);
}

void print_LF_map(struct LF_map map) {
  std::cout<<"  F   ...   L  "<<std::endl;
  for (int i = 0; i<(int)map.F.size(); i++) {
    std::cout<<"("<<std::get<0>(map.F[i])<<","<<std::get<1>(map.F[i])<<")"<<" ... ";
    std::cout<<"("<<std::get<0>(map.L[i])<<","<<std::get<1>(map.L[i])<<")"<<std::endl;
  }
  std::cout<<std::endl; 
}

void print_LF_map(std::string bwt_string) {
  struct LF_map map = make_LF_map(string);
  std::cout<<"  F   ...   L  "<<std::endl;
  for (int i = 0; i<(int)map.F.size(); i++) {
    std::cout<<"("<<std::get<0>(map.F[i])<<","<<std::get<1>(map.F[i])<<")"<<" ... ";
    std::cout<<"("<<std::get<0>(map.L[i])<<","<<std::get<1>(map.L[i])<<")"<<std::endl;
  }
  std::cout<<std::endl; 
}

std::string BWT(std::string string) {
  std::string pre = string + "$";
  int length = pre.size();
  std::vector<std::string> BWM;

  for (int i=0; i<(int)pre.size(); i++) {
    std::string s = pre.substr(i, length) + pre.substr(0, i);
    BWM.push_back(s);
  }

  std::sort(BWM.begin(), BWM.end());

  std::string bwt;
  for (int i=0; i<(int)BWM.size(); i++) {
    bwt = bwt + BWM[i][length-1];
  }

  return(bwt);
}

std::string inverse_BWT(struct LF_map map) {

  // Generate the original string.
  std::string original;
  std::vector<std::tuple<char, int>>::iterator i,j;
  
  char c;
  i = map.F.begin();
  j = map.L.begin();

  while(std::get<0>(*j) != '$') {
    int idx = std::distance(map.F.begin(), std::find(map.F.begin(), map.F.end(), *i));
    j = map.L.begin() + idx;

    c = std::get<0>(*j);
    original = c + original;

    i = std::find(map.F.begin(), map.F.end(), *j);
  }

  return(original.erase(0,1));
}

std::string inverse_BWT(std::string bwt_string) {
  struct LF_map map = make_LF_map(bwt_string);

  // Generate the original string.
  std::string original;
  std::vector<std::tuple<char, int>>::iterator i,j;
  
  char c;
  i = map.F.begin();
  j = map.L.begin();

  while(std::get<0>(*j) != '$') {
    int idx = std::distance(map.F.begin(), std::find(map.F.begin(), map.F.end(), *i));
    j = map.L.begin() + idx;

    c = std::get<0>(*j);
    original = c + original;

    i = std::find(map.F.begin(), map.F.end(), *j);
  }

  return(original.erase(0,1));
}

std::vector<std::string> substring_search(LF_map map, std::string substring) {
  std::vector<std::vector<int>> indices;
  // std::cout<<substring<<std::endl;

  // Get rows in the LF map that match the end of the query string.
  for (int j=0; j<(int)map.F.size(); j++) {
    int end = substring.size()-1;
    if (std::get<0>(map.F[j]) == substring[end]) {
      indices.push_back(std::vector<int>{j});
      // std::cout<<substring[end]<<": "<<j<<std::endl;
    }
  }

  int counter = (int)substring.size()-2;
  while ((indices.size() != 0) && counter >= 0) {
  /* while ((indices.size() != 0) && counter < (int)substring.size()) { */
    for (int i=0; i<(int)indices.size(); i++) {
      std::vector<int> *index = &indices[i];
      int idx = (*index)[index->size()-1];

      // std::cout<<"idx: "<<idx<<std::endl;
      // std::cout<<"new: "<<std::get<0>(map.L[idx])<<std::endl;

      if (std::get<0>(map.L[idx]) == substring[counter]) {
        // std::cout<<"yay: "<<substring[counter]<<std::endl;
        int new_idx = std::distance(map.F.begin(), std::find(map.F.begin(),
              map.F.end(), map.L[idx]));
        index->push_back(new_idx);
      } 
      else {
        index->empty();
        indices.erase(indices.begin() + i);
        i--;
      }
    }
    counter--;

    /* std::cout<<"counter:    "<<counter<<std::endl;
     * std::cout<<"index size: "<<indices.size()<<std::endl;
     * for (int i=0; i < (int)indices.size(); i++) {
     *   std::cout<<"index "<<i<<" size: "<<indices[i].size()<<std::endl;
     *   std::cout<<"   ";
     *   for (int j=0; j < (int)indices[i].size(); j++) {
     *     std::cout<<indices[i][j]<<" ";
     *   }
     *   std::cout<<std::endl;
     * } */
  }

  char c;
  // std::cout<<"n matches: "<<indices.size()<<std::endl;
  // std::string orig_sequence = inverse_BWT(map);
  std::vector<std::string> results;
  for (int i=0; i < (int)indices.size(); i++) {
    /* match.empty(); */
    std::string match;
    for (int j=0; j < (int)substring.size(); j++) {
      // std::cout<<std::get<0>(map.F[indices[i][j]]);
      c = std::get<0>(map.F[indices[i][j]]);
      match = c + match;
    }

    results.push_back(match);

    // // std::string prefix = std::string(indices[i][indices[i].size()-1]-1, '.');
    // std::string prefix = std::string(indices[i][indices[i].size()-1], '.');
    // std::cout<<indices[i][indices[i].size()-1]<<std::endl;
    // match = prefix + match + std::string(orig_sequence.size() - prefix.size() - match.size(), '.');
    // [> match = prefix + match; <]
    // 
    // std::cout<<orig_sequence<<std::endl;
    // std::cout<<match<<std::endl;
    // std::cout<<std::endl;
  }
  return results;
}

void print_substrings(std::vector<std::string> substrings) {
  for (auto it = substrings.begin(); it != substrings.end(); it++) {
    std::cout<<*it<<std::endl;
  }
}

int main(void) {

  // Obtain BWT string for "DOGWORD"
  std::cout<<"BW Transform:"<<std::endl;
  std::string bwt_string = BWT(string);
  std::cout<<"orig: "<<string<<std::endl;
  std::cout<<"bwt:  "<<bwt_string<<std::endl;
  std::cout<<std::endl;

  // Generate LF map
  struct LF_map map = make_LF_map(bwt_string);
  print_LF_map(map);

  // Obtain "DOGWORD" from the LF mapping of BWT
  std::cout<<"Inverse BWT by LF mapping:"<<std::endl;
  std::string lfmap_string = inverse_BWT(bwt_string);
  std::cout<<"orig:    "<<string<<std::endl;
  std::cout<<"bwt:     "<<bwt_string<<std::endl;
  std::cout<<"inv bwt: "<<lfmap_string<<std::endl; 
  std::cout<<std::endl;

  // Search substring via LF map
  std::cout<<"Search for substring 'GWO'"<<std::endl;
  // substring_search(map, "GWO");
  // std::vector<std::string> results = substring_search(map, "ANA");
  std::vector<std::string> results = substring_search(map, "GWO");
  print_substrings(results);

  return 0;
}
