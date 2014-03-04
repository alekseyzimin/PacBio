#ifndef _MER_POS_HASH_HPP_
#define _MER_POS_HASH_HPP_

#include <cstring>
#include <stdexcept>

#include <jellyfish/jellyfish.hpp>
#include <src_jf_aligner/lf_forward_list.hpp>

template<typename T, typename mer_type = jellyfish::mer_dna>
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
  typedef lf_forward_list_base<elt>            mapped_type;
  typedef typename mapped_type::iterator       pos_iterator;
  typedef typename mapped_type::const_iterator const_pos_iterator;

private:
  typedef typename mapped_type::head_node      head_node;
  typedef typename mapped_type::node           node;

public:
  mer_pos_hash(size_t size) :
    mers_(size, 2 * mer_type::k(), 0, 126),
    pos_(0)
    //    pos_(mers_.size())
  {
    pos_ = new head_node[mers_.size()];
    memset(pos_, '\0', sizeof(head_node) * mers_.size());
  }

  ~mer_pos_hash() {
    delete [] pos_;
    for(auto it = data_pages.begin(); it != data_pages.end(); ++it)
      delete [] *it;
  }

  // /** Returns a reference to the list of positions associated with mer
  //     m */
  // mapped_type operator[](const mer_type& m) {
  //   bool   is_new;
  //   size_t id;
  //   if(!mers_.set(m, &is_new, &id))
  //     throw std::runtime_error("Mer hash table is full");
  //   return mapped_type(pos_[id]);
  // }

  head_node* insert_mer(const mer_type& m) {
    bool   is_new;
    size_t id;
    if(!mers_.set(m, &is_new, &id))
      throw std::runtime_error("Mer hash table is full");
    return &pos_[id];
  }

  /** Push a position for mer m. */
  void push_front(const mer_type& m, const T* s, int o) {
    struct data_page {
      node*  data;
      size_t used;
    };
    static const size_t       node_per_page = 1024 * 1024;
    static __thread data_page thread_data   = { 0, 0 };

    if(!thread_data.data || thread_data.used >= node_per_page) {
      thread_data.used = 0;
      thread_data.data = new node[node_per_page];
      data_pages.push_front(thread_data.data);
    }
    node* cnode  = &thread_data.data[thread_data.used++];
    cnode->val_.frag   = s;
    cnode->val_.offset = o;
    mapped_type::insert_after_(insert_mer(m), cnode);
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

private:
  mer_array                         mers_;
  typename mapped_type::head_node*  pos_;
  lf_forward_list<node*>            data_pages;
};

#endif /* _MER_POS_HASH_HPP_ */
