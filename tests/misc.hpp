#ifndef __MISC_H__
#define __MISC_H__

#include <stdlib.h>

namespace misc {
char rev_comp(const char c);
std::string rev_comp(const std::string s);

struct remove_file {
  const char* path;
  bool do_unlink;
  remove_file() : path(getenv("TEST_TMP")), do_unlink(false) { }
  remove_file(const char* p, bool unlink = true) : path(p), do_unlink(unlink) { }
  ~remove_file() { if(do_unlink) unlink(path); }
};

} // namespace misc

#endif /* __MISC_H__ */
