#include <stdexcept>

#include <jellyfish/jellyfish.hpp>
#include <src_jf_aligner/lf_forward_list.hpp>

template<typename mer_type = jellyfish::mer_dna>
class mer_pos_hash {
  struct elt {
    const std::string* frag;
    unsigned int offset;
    // Ordering of the element based on there offset. Only meaningful
    // if the frags are equal, which is assumed.
    bool operator<(const elt& rhs) const { return offset < rhs.offset; }
    elt() = default;
    elt(const std::string& s, unsigned int o) : frag(&s), offset(o) { }
    elt(const std::string* s, unsigned int o) : frag(s), offset(o) { }
  };

  // List type containing the positions of the mers


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
  void push_front(const mer_type& m, const std::string& s, unsigned int o) {
    (*this)[m].push_front(position_type(s, o));
  }

private:
  mer_array                mers_;
  std::vector<mapped_type> pos_;
};
