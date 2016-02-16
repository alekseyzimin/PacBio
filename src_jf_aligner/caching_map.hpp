/******************************************
Copyright University of Maryland 2015
******************************************/
#ifndef __CACHING_MAP_H__
#define __CACHING_MAP_H__

#include <vector>
#include <unordered_map>
#include <iterator>
#include <utility>

/**
 * Un-ordered map with caching. The value objects in the map are not
 * deleted on a clear() or erase. When inserting a new key, the value
 * may already contain a recycled object. It is the responsibility of
 * the user to properly initialize the value.
 *
 * On the other hand, the value objects are deleted upon destruction
 * of the map.
 */

template<class Key, class T,
         class Hash  = std::hash<Key>,
         class Pred  = std::equal_to<Key>,
         class Alloc = std::allocator< std::pair<const Key, int> >
         >
class caching_map {
  //  typedef std::vector<T, typename Alloc::rebind<T>::other> vector_type;
  typedef std::vector<T> vector_type; // Missing allocator stuff
  typedef std::unordered_map<Key, int, Hash, Pred, Alloc>  map_type;

public:
  typedef Key                                      key_type;
  typedef const Key                                const_key_type;
  typedef T                                        mapped_type;
  typedef std::pair<const_key_type&, mapped_type&> value_type;
  typedef Hash                                     hasher;
  typedef Pred                                     key_equal;
  typedef Alloc                                    allocator_type;
  typedef value_type&                              reference;
  typedef const value_type&                        const_reference;
  typedef value_type*                              pointer;
  typedef const value_type*                        const_pointer;
  typedef typename map_type::size_type             size_type;
  typedef typename map_type::difference_type       difference_type;

private:
  vector_type values_;
  map_type    map_;
  size_type   index_;

public:
  class const_iterator;
  class iterator : public std::iterator<std::forward_iterator_tag, value_type> {
    vector_type*                values_;
    typename map_type::iterator it_;
  public:
    iterator() = default;
    iterator(vector_type* v, typename map_type::iterator it) : values_(v), it_(it) { }
    iterator(const iterator& rhs) : values_(rhs.values_), it_(rhs.it_) { }
    iterator& operator=(const iterator& rhs) {
      values_ = rhs.values_;
      it_     = it_.values_;
      return *this;
    }
    bool operator==(const iterator& rhs) const { return it_ == rhs.it_ && values_ == rhs.values_; }
    bool operator!=(const iterator& rhs) const { return it_ != rhs.it_ || values_ != rhs.values_; }
    bool operator==(const const_iterator& rhs) const { return it_ == rhs.it_ && values_ == rhs.values_; }
    bool operator!=(const const_iterator& rhs) const { return it_ != rhs.it_ || values_ != rhs.values_; }
    value_type operator*() { return value_type(it_->first, (*values_)[it_->second]); }
    //    value_type* operator->() { return &operator*(); }
    iterator& operator++() { ++it_; return *this; }
    iterator operator++(int) { iterator res(*this); operator++(); return res; }
  };

  class const_iterator : public std::iterator<std::forward_iterator_tag, const value_type> {
    const vector_type*                values_;
    typename map_type::const_iterator it_;
  public:
    const_iterator() = default;
    const_iterator(const vector_type* v, const typename map_type::const_iterator it) : values_(v), it_(it) { }
    const_iterator(const const_iterator& rhs) : values_(rhs.values_), it_(rhs.it_) { }
    const_iterator(const iterator& rhs) : values_(rhs.values_), it_(rhs.it_) { }
    const_iterator& operator=(const const_iterator& rhs) {
      values_ = rhs.values_;
      it_     = it_.values_;
      return *this;
    }
    bool operator==(const const_iterator& rhs) const { return it_ == rhs.it_ && values_ == rhs.values_; }
    bool operator!=(const const_iterator& rhs) const { return it_ != rhs.it_ || values_ != rhs.values_; }
    bool operator==(const iterator& rhs) const { return it_ == rhs.it_ && values_ == rhs.values_; }
    bool operator!=(const iterator& rhs) const { return it_ != rhs.it_ || values_ != rhs.values_; }
    const value_type operator*() const { return value_type(it_->first, (*values_)[*it_]); }
    //    const value_type* operator->() const { return &operator*(); }
    iterator& operator++() { ++it_; return *this; }
    iterator operator++(int) { iterator res(*this); operator++(); return res; }
  };

  caching_map() :
    values_(), map_(), index_(0)
  { }
  explicit caching_map(size_type n, const hasher& hf = hasher(), const key_equal& eql = key_equal(),
                       const allocator_type& alloc = allocator_type()) :
    values_(alloc),
    map_(n, hf, eql, alloc),
    index_(0)
  { }
  explicit caching_map(const allocator_type& alloc) :
    values_(alloc),
    map_(alloc),
    index_(0)
  { }
  caching_map(const caching_map& rhs) :
    values_(rhs.values_),
    map_(rhs.map_),
    index_(rhs.index_)
  { }
  caching_map(const caching_map& rhs, const allocator_type& alloc) :
    values_(rhs.values_, alloc),
    map_(rhs.map_, alloc),
    index_(rhs.index_)
  { }
  caching_map(caching_map&& rhs) :
    values_(std::move(rhs.values_)),
    map_(std::move(rhs.map_)),
    index_(rhs.index_)
  { }
  caching_map(caching_map&& rhs, const allocator_type& alloc) :
    values_(std::move(rhs.values_), alloc),
    map_(std::move(rhs.map_), alloc),
    index_(rhs.index_)
  { }

  iterator find(const key_type& k) { return iterator(&values_, map_.find(k)); }
  const_iterator find(const key_type& k) const { return const_iterator(&values_, map_.find(k)); }

  iterator begin() { return iterator(&values_, map_.begin()); }
  const_iterator begin() const { return const_iterator(&values_, map_.begin()); }
  const_iterator cbegin() const { return const_iterator(&values_, map_.cbegin()); }

  iterator end() { return iterator(&values_, map_.end()); }
  const_iterator end() const { return const_iterator(&values_, map_.end()); }
  const_iterator cend() const { return const_iterator(&values_, map_.cend()); }

  void clear() noexcept { map_.clear(); index_ = 0; }

  std::pair<iterator, bool> insert(const value_type& val) {
    auto it = map_.find(val.first);
    if(it != map_.end()) return std::make_pair(iterator(&values_, it), false);
    if(index_ >= values_.size()) {
      values_.push_back(val.second);
    } else {
      values_[index_] = val.second;
    }
    std::cout << "inserted at " << (void*)&values_[index_] << std::endl;
    ++index_;
    return std::make_pair(iterator(&values_, map_.insert(it, std::make_pair(&val.first, index_ - 1))), true);
  }
  template<typename P>
  std::pair<iterator, bool> insert(P&& val) {
    auto it = map_.find(val.first);
    if(it != map_.end()) return std::make_pair(iterator(&values_, it), false);
    if(index_ >= values_.size()) {
      values_.push_back(std::move(val.second));
    } else {
      values_[index_] = std::move(val.second);
    }
    std::cout << "inserted at " << (void*)&values_[index_] << std::endl;
    ++index_;
    return std::make_pair(iterator(&values_, map_.insert(it, std::make_pair(std::move(val.first), index_ - 1))), true);
  }
};

#endif /* __CACHING_MAP_H__ */
