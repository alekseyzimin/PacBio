#ifndef _MER_POS_HASH_HPP_
#define _MER_POS_HASH_HPP_

#include <cstring>
#include <stdexcept>

#include <jellyfish/hash_counter.hpp>
#include <jellyfish/allocators_mmap.hpp>
#include <src_jf_aligner/lf_forward_size_list.hpp>

template<typename T, typename mer_type = jellyfish::mer_dna, typename mer_array_type = jellyfish::large_hash::array<mer_type> >
class mer_pos_hash {
  // The offset is 1-based, >0 if forward, <0 if backward. This is for
  // the canonical representation of a k-mer with respect to the
  // sequence in the fragment.
  struct elt {
    const T* frag;
    int      offset;
    // Ordering of the element based on there offset. Only meaningful
    // if the frags are equal, which is assumed.
    bool operator<(const elt& rhs) const { return offset < rhs.offset; }
    elt() = default;
    //    elt(const std::string& s, int o) : frag(&s), offset(o) { }
    elt(const T* s, int o) : frag(s), offset(o) { }
  };

public:
  typedef elt                                  position_type;
  typedef mer_type                             key_type;
  typedef lf_forward_size_list_base<elt>       mapped_type;
  typedef typename mapped_type::iterator       pos_iterator;
  typedef typename mapped_type::const_iterator const_pos_iterator;

private:
  typedef typename mapped_type::head_node      head_node;
  typedef typename mapped_type::node           node;

  // For paged per thread memory allocation
  struct data_page {
    node*  data;
    size_t used;
  };
  static const size_t node_per_page = 1024 * 1024;

public:
  mer_pos_hash(size_t size, uint16_t max_count = std::numeric_limits<uint16_t>::max()) :
    mers_(size, 2 * mer_type::k(), 0, 126),
    data_(mers_.size() * sizeof(head_node)),
    pos_((head_node*)data_.get_ptr()),
    max_count_(max_count)
  {
    if(!pos_)
      throw std::runtime_error("Can't allocate mer pos hash");
  }

  ~mer_pos_hash() {
    //    delete [] pos_;
    for(auto it = data_pages_.begin(); it != data_pages_.end(); ++it)
      delete [] *it;
  }

  head_node* insert_mer(const mer_type& m) {
    bool   is_new;
    size_t id;
    if(!mers_.set(m, &is_new, &id))
      throw std::runtime_error("Mer hash table is full");
    return &pos_[id];
  }

  pos_iterator pos_end() { return pos_iterator(); }
  const_pos_iterator pos_end() const { return const_pos_iterator(); }

  /** Find the mer m in the hash. Return NULL if not found, a pointer
      to the list of positions if found. */
  pos_iterator find_pos(const mer_type& m) {
    key_type tmp_m;
    return find_pos(m, tmp_m);
  }

  const_pos_iterator find_pos(const mer_type& m) const {
    key_type tmp_m;
    return find_pos(m, tmp_m);
  }

  /** Find the mer m in the hash. Return NULL if not found, a pointer
      to the list of positions if found. Optimization. */
  pos_iterator find_pos(const mer_type& m, mer_type& tmp_m) {
    size_t id;
    if(mers_.get_key_id(m, &id, tmp_m))
      return pos_iterator(pos_[id].next_);
    return pos_iterator();
  }

  pos_iterator find_pos(const mer_type& m, mer_type& tmp_m) const {
    size_t id;
    if(mers_.get_key_id(m, &id, tmp_m))
      return pos_iterator(pos_[id].next_);
    return pos_iterator();
  }

  std::pair<pos_iterator, uint16_t> find_pos_size(const mer_type& m) const {
    mer_type tmp_m;
    size_t   id;
    if(mers_.get_key_id(m, &id, tmp_m)) {
      const auto& list = pos_[id];
      return std::pair<pos_iterator, uint16_t>(pos_iterator(list.next_), list.size_);
    }
    return std::pair<pos_iterator, uint16_t>(pos_iterator(), 0);
  }

  /** Same than find_pos but interface compatible with multimap.
   */
  std::pair<pos_iterator, pos_iterator> equal_range(const mer_type& m) {
    return std::make_pair(find_pos(m), pos_iterator());
  }

  std::pair<const_pos_iterator, const_pos_iterator> equal_range(const mer_type& m) const {
    return std::make_pair(find_pos(m), const_pos_iterator());
  }

  class thread {
    mer_pos_hash& hash_;
    data_page     data_;

  public:
    thread(mer_pos_hash& hash) : hash_(hash), data_({0, 0}) { }
    thread(mer_pos_hash* hash) : hash_(*hash), data_({0, 0}) { }

    /** Push a position for mer m. */
    void push_front(const mer_type& m, const T* s, int o) {
      if(!data_.data || data_.used >= mer_pos_hash::node_per_page) {
        data_.used = 0;
        data_.data = new node[node_per_page];
        hash_.data_pages_.push_front(data_.data);
      }
      node* cnode  = &data_.data[data_.used++];
      cnode->val_.frag   = s;
      cnode->val_.offset = o;
      auto list = hash_.insert_mer(m);
      mapped_type::insert_after_(list, cnode, hash_.max_count_);
    }
  };
  friend class thread;

private:
  mer_array_type                   mers_;
  allocators::mmap                 data_;
  typename mapped_type::head_node* pos_;
  lf_forward_list<node*>           data_pages_;
  const uint16_t                   max_count_;
};

#endif /* _MER_POS_HASH_HPP_ */

