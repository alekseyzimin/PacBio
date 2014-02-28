#ifndef __FRAG_INFO_H__
#define __FRAG_INFO_H__

// Multiple vector to hold names of fragments. There is a vector for
// each thread and the string is copied when push_back is called.
class frag_lists {
public:
  struct frag_info {
    unsigned int len;
    char         name[];
  };

  frag_lists(size_t threads) : names_(threads) { }
  ~frag_lists() {
    for(auto it = names_.cbegin(); it != names_.cend(); ++it)
      for(auto it2 = it->cbegin(); it2 != it->cend(); ++it2)
        ::operator delete((void*)*it2);
  }

  size_t size() const { return names_.size(); }

  void ensure(size_t threads) {
    if(names_.size() < threads)
      names_.resize(threads);
  }

  const frag_info* push_back(int thid, unsigned int len, const char* s) {
    auto fi = (frag_info*) ::operator new(sizeof(frag_info) + strlen(s) + 1);
    fi->len = len;
    strcpy(fi->name, s);
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
