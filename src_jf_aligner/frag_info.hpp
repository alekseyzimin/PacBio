#ifndef __FRAG_INFO_H__
#define __FRAG_INFO_H__

#include <src_jf_aligner/super_read_name.hpp>


// Multiple vector to hold names of fragments. There is a vector for
// each thread and the string is copied when push_back is called.
class frag_lists {
public:
  struct name_unitigs {
    std::string     name;
    super_read_name unitigs;
  };
  struct frag_info {
    unsigned int    len;
    name_unitigs    fwd;
    name_unitigs    bwd;
    frag_info(unsigned int l, const char* s)
      : len(l)
    {
      fwd.name    = s;
      fwd.unitigs = fwd.name;
      bwd.unitigs = fwd.unitigs;
      bwd.unitigs.reverse();
      bwd.name    = bwd.unitigs.name();
    }
  };

  frag_lists(size_t threads) : names_(threads) { }
  ~frag_lists() {
    for(auto it = names_.cbegin(); it != names_.cend(); ++it)
      for(auto it2 = it->cbegin(); it2 != it->cend(); ++it2)
        delete *it2;
    //        ::operator delete((void*)*it2);
  }

  size_t size() const { return names_.size(); }

  void ensure(size_t threads) {
    if(names_.size() < threads)
      names_.resize(threads);
  }

  void clear(size_t thid) {
    for(auto it2 = names_[thid].cbegin(); it2 != names_[thid].cend(); ++it2)
      delete *it2;
  }

  const frag_info* push_back(int thid, unsigned int len, const char* s) {
    auto fi = new frag_info(len, s);
    names_[thid].push_back(fi);
    return fi;
  }

  const frag_info* push_back(int thid, unsigned int len, const std::string& s) {
    return push_back(thid, len, s.c_str());
  }

  const std::vector<const frag_info*> operator[](int i) const {
    return names_[i];
  }

private:
  std::vector<std::vector<const frag_info*> > names_;
};

#endif /* __FRAG_INFO_H__ */
