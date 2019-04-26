#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

/* Hard-code the string for the BWT. */
static const std::string main_string = "DOGWORD";

/* Create a struct for the LF mapping. There are two vectors of tuples, F and
 * L. F contains tuples for the first column, where the tuples are: (letter,
 * rank). L is the same but for the last column. */
struct LF_map
{
  std::vector<std::tuple<char, int>> F;
  std::vector<std::tuple<char, int>> L;
};

/* Create a function to generate an LF map from a string produced by the BWT.
 * */
struct LF_map
make_LF_map(std::string bwt_string)
{
  std::string sorted_string = bwt_string;
  std::sort(sorted_string.begin(), sorted_string.end());

  // Loop over the the array and use a counter to set the indices for all the
  // letters.
  std::vector<std::tuple<char, int>> sorted;
  for (int i = 0; i < (int)sorted_string.size(); i++) {
    int counter = 0;
    for (int j = 0; j < i; j++) {
      if (sorted_string[i] == sorted_string[j]) {
        counter++;
      }
    }
    sorted.push_back(std::tuple<char, int>(sorted_string[i], counter));
  }

  // Repeat the above, but this time for the last column (L).
  std::vector<std::tuple<char, int>> bwt;
  for (int i = 0; i < (int)bwt_string.size(); i++) {
    int counter = 0;
    for (int j = 0; j < i; j++) {
      if (bwt_string[i] == bwt_string[j]) {
        counter++;
      }
    }
    bwt.push_back(std::tuple<char, int>(bwt_string[i], counter));
  }

  struct LF_map map;
  map.F = sorted;
  map.L = bwt;

  return (map);
}

/* Utility function for printing the LF mapping. */
void
print_LF_map(struct LF_map map)
{
  std::cout << "  F   ...   L  " << std::endl;
  for (int i = 0; i < (int)map.F.size(); i++) {
    std::cout << "(" << std::get<0>(map.F[i]) << "," << std::get<1>(map.F[i])
              << ")"
              << " ... ";
    std::cout << "(" << std::get<0>(map.L[i]) << "," << std::get<1>(map.L[i])
              << ")" << std::endl;
  }
  std::cout << std::endl;
}

/* Utility function for printing the LF mapping. Take as string as input. */
void
print_LF_map(std::string bwt_string)
{
  struct LF_map map = make_LF_map(main_string);
  std::cout << "  F   ...   L  " << std::endl;
  for (int i = 0; i < (int)map.F.size(); i++) {
    std::cout << "(" << std::get<0>(map.F[i]) << "," << std::get<1>(map.F[i])
              << ")"
              << " ... ";
    std::cout << "(" << std::get<0>(map.L[i]) << "," << std::get<1>(map.L[i])
              << ")" << std::endl;
  }
  std::cout << std::endl;
}

/* Compute the BWT for a string. */
std::string
BWT(std::string string)
{
  std::string pre = string + "$";
  int length = pre.size();
  std::vector<std::string> BWM;

  /* Create the BWT matrix. */
  for (int i = 0; i < (int)pre.size(); i++) {
    std::string s = pre.substr(i, length) + pre.substr(0, i);
    BWM.push_back(s);
  }

  std::sort(BWM.begin(), BWM.end());

  /* Take the last column of the sorted BWT matrix. */
  std::string bwt;
  for (int i = 0; i < (int)BWM.size(); i++) {
    bwt = bwt + BWM[i][length - 1];
  }

  return (bwt);
}

/* Compute the inverse BWT from the LF mapping. */
std::string
inverse_BWT(struct LF_map map)
{

  // Generate the original string.
  std::string original;
  std::vector<std::tuple<char, int>>::iterator i, j;

  char c;
  i = map.F.begin();
  j = map.L.begin();

  while (std::get<0>(*j) != '$') {
    int idx =
      std::distance(map.F.begin(), std::find(map.F.begin(), map.F.end(), *i));
    j = map.L.begin() + idx;

    c = std::get<0>(*j);
    original = c + original;

    i = std::find(map.F.begin(), map.F.end(), *j);
  }

  return (original.erase(0, 1));
}

/* Compute the inverse BWT from the BWT string via LF mapping. */
std::string
inverse_BWT(std::string bwt_string)
{
  struct LF_map map = make_LF_map(bwt_string);

  // Generate the original string.
  std::string original;
  std::vector<std::tuple<char, int>>::iterator i, j;

  char c;
  i = map.F.begin();
  j = map.L.begin();

  while (std::get<0>(*j) != '$') {
    int idx =
      std::distance(map.F.begin(), std::find(map.F.begin(), map.F.end(), *i));
    j = map.L.begin() + idx;

    c = std::get<0>(*j);
    original = c + original;

    i = std::find(map.F.begin(), map.F.end(), *j);
  }

  return (original.erase(0, 1));
}

/* Use the LF mapping to perform a substring search. */
std::vector<std::string>
substring_search(LF_map map, std::string substring)
{
  std::vector<std::vector<int>> indices;

  // Get rows in the LF map that match the end of the query string.
  for (int j = 0; j < (int)map.F.size(); j++) {
    int end = substring.size() - 1;
    if (std::get<0>(map.F[j]) == substring[end]) {
      indices.push_back(std::vector<int>{ j });
    }
  }

  // Iterate backward over the search string. Loop over the LF table and until
  // the search sequence length is reached or a viable subsequence is not
  // possible to find.
  int counter = (int)substring.size() - 2;
  while ((indices.size() != 0) && counter >= 0) {
    for (int i = 0; i < (int)indices.size(); i++) {
      std::vector<int>* index = &indices[i];
      int idx = (*index)[index->size() - 1];

      // If the L value matches the character in the substring, update the
      // index for the F column and continue iterating. If there is no match,
      // remove the path from the list of viable paths.
      if (std::get<0>(map.L[idx]) == substring[counter]) {
        int new_idx = std::distance(
          map.F.begin(), std::find(map.F.begin(), map.F.end(), map.L[idx]));
        index->push_back(new_idx);
      } else {
        index->empty();
        indices.erase(indices.begin() + i);
        i--;
      }
    }
    counter--;
  }

  // Generate a list of paths from the viable paths found in the above loop.
  // Because the search sequences are found back-to-front, the character
  // sequences need to be re-inverted to get the correct match.
  char c;
  std::vector<std::string> results;
  for (int i = 0; i < (int)indices.size(); i++) {
    std::string match;
    for (int j = 0; j < (int)substring.size(); j++) {
      c = std::get<0>(map.F[indices[i][j]]);
      match = c + match;
    }

    results.push_back(match);
  }
  return results;
}

/* Utility function for iterating over a vector and printing the contents. */
void
print_substrings(std::vector<std::string> substrings)
{
  for (auto it = substrings.begin(); it != substrings.end(); it++) {
    std::cout << *it << std::endl;
  }
}

/* Filled with fluffy polar bears. */
int
main(void)
{

  // Obtain BWT string for "DOGWORD"
  std::cout << "BW Transform:" << std::endl;
  std::string bwt_string = BWT(main_string);
  std::cout << "orig: " << main_string << std::endl;
  std::cout << "bwt:  " << bwt_string << std::endl;
  std::cout << std::endl;

  // Generate LF map
  struct LF_map map = make_LF_map(bwt_string);
  print_LF_map(map);

  // Obtain "DOGWORD" from the LF mapping of BWT
  std::cout << "Inverse BWT by LF mapping:" << std::endl;
  std::string lfmap_string = inverse_BWT(bwt_string);
  std::cout << "orig:    " << main_string << std::endl;
  std::cout << "bwt:     " << bwt_string << std::endl;
  std::cout << "inv bwt: " << lfmap_string << std::endl;
  std::cout << std::endl;

  // Search substring via LF map
  std::cout << "Search for substring 'GWO'" << std::endl;
  std::vector<std::string> results = substring_search(map, "GWO");
  print_substrings(results);

  return 0;
}
