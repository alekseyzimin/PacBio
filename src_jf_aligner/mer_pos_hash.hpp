#ifndef _MER_POS_HASH_HPP_
#define _MER_POS_HASH_HPP_

#include <stdexcept>

#include <jellyfish/jellyfish.hpp>
#include <src_jf_aligner/lf_forward_list.hpp>

template<typename mer_type = jellyfish::mer_dna>
class mer_pos_hash {
  // The offset is 1-based, >0 if forward, <0 if backward. This is for
  // the canonical representation of a k-mer with respect to the
  // sequence in the fragment.
  struct elt {
    const char* frag;
    int         offset;
    // Ordering of the element based on there offset. Only meaningful
    // if the frags are equal, which is assumed.
    bool operator<(const elt& rhs) const { return offset < rhs.offset; }
    elt() = default;
    //    elt(const std::string& s, int o) : frag(&s), offset(o) { }
    elt(const char* s, int o) : frag(s), offset(o) { }
  };

public:
  typedef elt                  position_type;
  typedef mer_type             key_type;
  typedef lf_forward_list<elt> mapped_type;

  mer_pos_hash(size_t size) :
    mers_(size, 2 * mer_type::k(), 0, 126),
    pos_(mers_.size())
  { }

  /** Returns a reference to the list of positions associated with mer
      m */
  mapped_type& operator[](const mer_type& m) {
    bool   is_new;
    size_t id;
    if(!mers_.set(m, &is_new, &id))
      throw std::runtime_error("Mer hash table is full");
    return pos_[id];
  }

  /** Push a position for mer m. */
  void push_front(const mer_type& m, const char* s, int o) {
    (*this)[m].push_front(position_type(s, o));
  }

  /** Find the mer m in the hash. Return NULL if not found, a pointer
      to the list of positions if found. */
  mapped_type* find_pos(const mer_type& m) {
    key_type tmp_m;
    return find_pos(m, tmp_m);
  }

  const mapped_type* find_pos(const mer_type& m) const {
    key_type tmp_m;
    return find_pos(m, tmp_m);
  }

  /** Find the mer m in the hash. Return NULL if not found, a pointer
      to the list of positions if found. Optimization. */
  mapped_type* find_pos(const mer_type& m, mer_type& tmp_m) {
    size_t id;
    if(mers_.get_key_id(m, &id, tmp_m))
      return &pos_[id];
    return 0;
  }

  const mapped_type* find_pos(const mer_type& m, mer_type& tmp_m) const {
    size_t id;
    if(mers_.get_key_id(m, &id, tmp_m))
      return &pos_[id];
    return 0;
  }

private:
  mer_array                mers_;
  std::vector<mapped_type> pos_;
};

#endif /* _MER_POS_HASH_HPP_ */
