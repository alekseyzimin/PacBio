#include <vector>
#include <iostream>

#include "lis.hpp"

int main(int argc, char *argv[])
{
  typedef float input_type;
  std::vector<input_type> input;
  while(true) {
    input_type x;
    std::cin >> x;
    if(!std::cin.good())
      break;
    input.push_back(x);
  }
  size_t lis_size = lis::length(input.cbegin(), input.cend());
  std::cout << lis_size << "\n";

  auto lis_seq = lis::sequence(input.cbegin(), input.cend());
  for(auto it = lis_seq.cbegin(); it != lis_seq.cend(); ++it)
    std::cout << *it << " ";
  std::cout << "\n";

  auto lis_indices = lis::indices(input.cbegin(), input.cend());
  for(auto it = lis_indices.cbegin(); it != lis_indices.cend(); ++it)
    std::cout << input[*it] << " ";
  std::cout << "\n";

  return 0;
}
