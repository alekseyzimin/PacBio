/******************************************
Copyright University of Maryland 2015
******************************************/
#ifndef __UNION_FIND_H__
#define __UNION_FIND_H__

namespace union_find {
struct set {
  int  rank;
  set* parent;
  set() : rank(0), parent(this) { }
  // set(const set& s) = delete;
  // set(set&& s) = delete;
  inline set* root(); // equivalent to find_root
  bool operator==(set& rhs) { return root() == rhs.root(); }
  inline set& operator|=(set& rhs); // equivalent to union_sets
  void reset() {
    rank = 0;
    parent = this;
  }
};

void union_sets(set* s1, set* s2);
set* find_root(set* s);
inline void union_sets(set& s1, set& s2) { union_sets(&s1, &s2); }
inline set* find_root(set& s) { return find_root(&s); }
inline set* set::root() { return find_root(this); }
inline set& set::operator|=(set& rhs) {
  union_sets(this, &rhs);
  return *this;
}
} // namespace union_find

#endif /* __UNION_FIND_H__ */
