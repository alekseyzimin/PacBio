#include <string>

namespace misc {
char rev_comp(const char c) {
  switch(c) {
  case 'A': case 'a': return 'T';
  case 'C': case 'c': return 'G';
  case 'G': case 'g': return 'C';
  case 'T': case 't': return 'A';
  }
  return 'N';
}
std::string rev_comp(const std::string s) {
  std::string res;
  for(auto it = s.crbegin(); it != s.crend(); ++it)
    res += rev_comp(*it);
  return res;
}
} // namespace misc
