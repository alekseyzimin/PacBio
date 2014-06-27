#include <src_jf_aligner/union_find.hpp>

void union_find::union_sets(union_find::set* s1, union_find::set* s2) {
  auto r1 = find_root(s1);
  auto r2 = find_root(s2);

  if(r1->rank > r2->rank) {
    r2->parent = s1;
  } else if(r1->rank < r2->rank) {
    r1->parent = s2;
  } else if(s1 != s2) {
    r2->parent = s1;
    ++r1->rank;
  }
}

union_find::set* union_find::find_root(union_find::set* s) {
  if(s != s->parent)
    s->parent = find_root(s->parent);
  return s->parent;
}
