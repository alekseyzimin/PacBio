#ifndef __MISC_H__
#define __MISC_H__

#include <sys/types.h>
#include <unistd.h>
#include <libgen.h>
#include <stdlib.h>
#include <string>
#include <vector>

namespace misc {
char rev_comp(const char c);
std::string rev_comp(const std::string s);

struct remove_file {
  const char* path;
  std::string str_path;
  bool do_unlink;
  remove_file() : do_unlink(false) {
    std::vector<char> buf(4096, '\0');
    std::string exe("/proc/self/exe");
    readlink(exe.c_str(), buf.data(), buf.size() - 1);
    str_path = std::string("tests/") + basename(buf.data()) + ".tmp";
    path     = str_path.c_str();
  }
  remove_file(const char* p, bool unlink = true) : path(p), do_unlink(unlink) { }
  ~remove_file() { if(do_unlink) unlink(path); }
};

} // namespace misc

#endif /* __MISC_H__ */
