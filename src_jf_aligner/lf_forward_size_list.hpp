/******************************************
Copyright University of Maryland 2015
******************************************/
#ifndef _LF_FORWARD_SIZE_LIST_HPP_
#define _LF_FORWARD_SIZE_LIST_HPP_

#include <memory>
#include <atomic>

#include <src_jf_aligner/lf_forward_list.hpp>
#include <src_jf_aligner/bounded_counter.hpp>


template<typename T, bool LF = true, typename U = uint16_t>
struct lf_forward_size_list_base {
  struct base_node {
    mutable base_node* next_;
  };
  struct head_node : base_node {
    bounded_counter<LF, U> size_;
  };
  struct node : base_node {
    T val_;
  };

  lf_forward_size_list_base(base_node* h)  { head_.next_ = h; }

  class iterator;
  // Const iterator
  class const_iterator : public std::iterator<std::forward_iterator_tag, T> {
  public:
    friend class iterator;
    friend struct lf_forward_size_list_base<T, LF>;
    explicit const_iterator(const base_node* head = 0) : head_(head) { }
    const_iterator(const const_iterator& rhs) : head_(rhs.head_) { }
    const_iterator(const iterator& rhs) : head_(rhs.head_) { }
    const_iterator& operator=(const const_iterator& rhs) {
      head_ = rhs.head_;
      return *this;
    }
    bool operator==(const const_iterator& rhs) const { return head_ == rhs.head_; }
    bool operator!=(const const_iterator& rhs) const { return head_ != rhs.head_; }
    const T& operator*() const { return static_cast<const node*>(head_)->val_; }
    const T* operator->() const { return &static_cast<const node*>(head_)->val_; }
    const_iterator& operator++() { head_ = (*const_cast<base_node* volatile*>(&head_))->next_; return *this; }
    const_iterator operator++(int) {
      const_iterator res(*this);
      ++*this;
      return res;
    }
  private:
    const base_node* head_;
  };

  // Iterator
  class iterator : public std::iterator<std::forward_iterator_tag, T> {
  public:
    friend class const_iterator;
    friend struct lf_forward_size_list_base<T, LF>;
    explicit iterator(base_node* head = 0) : head_(head) { }
    iterator(const iterator& rhs) : head_(rhs.head_) { }
    iterator& operator=(const iterator& rhs) {
      head_ = rhs.head_;
      return *this;
    }
    bool operator==(const iterator& rhs) const { return head_ == rhs.head_; }
    bool operator!=(const iterator& rhs) const { return head_ != rhs.head_; }
    bool operator==(const const_iterator& rhs) const { return head_ == rhs.head_; }
    bool operator!=(const const_iterator& rhs) const { return head_ != rhs.head_; }
    T& operator*() { return static_cast<node*>(head_)->val_; }
    const T& operator*() const { return static_cast<const node*>(head_)->val_; }
    T* operator->() { return &static_cast<node*>(head_)->val_; }
    const T* operator->() const { return &static_cast<const node*>(head_)->val_; }
    iterator& operator++() { head_ = (*const_cast<base_node* volatile*>(&head_))->next_; return *this; }
    iterator operator++(int) {
      iterator res(*this);
      ++*this;
      return res;
    }
  private:
    base_node* head_;
  };

  /** Get an iterator */
  iterator begin() { return iterator(head_.next_); }
  const_iterator begin() const { return const_iterator(head_.next_); }
  const_iterator cbegin() const { return const_iterator(head_.next_); }
  iterator before_begin() { return iterator(&head_); }
  const_iterator before_begin() const { return const_iterator(&head_); }
  const_iterator cbefore_begin() const { return const_iterator(&head_); }

  iterator end() { return iterator(); }
  const_iterator end() const { return const_iterator(); }
  const_iterator cend() const { return const_iterator(); }

  static void insert_after_(const base_node* head, base_node* n) {
    if(LF) {
      base_node* hn = *const_cast<base_node* volatile*>(&(head->next_));
      do {
        n->next_ = hn;
        hn = __sync_val_compare_and_swap(&head->next_, hn, n);
      } while(hn != n->next_);
    } else {
      n->next_    = head->next_;
      head->next_ = n;
    }
  }

  static void insert_after_(head_node* head, base_node* n, const U max_count) {
    if(head->size_.inc(max_count))
      insert_after_(static_cast<base_node*>(head), n);
  }

  static void insert_after_(const_iterator pos, base_node* n) {
    insert_after_(pos.head_, n);
  }

  iterator push_front_(base_node* n, const U max) {
    insert_after_(&head_, n);
    return iterator(n);
  }

  // Remove elements after position up to last. Somewhat thread
  // safe. It can be mixed with insert operation but not with erase
  // operations on ranges that overlap. The base class does not
  // allocate/deallocate the elements, hence it does not control the
  // life span of the element removed.
  // Returns a pointer to the first node removed from the list (position->next_).
  static base_node* erase_after_(base_node* position, base_node* last) {
    base_node* nn;
    if(LF) {
      base_node* onn;
      nn = *const_cast<base_node* volatile*>(&position->next_);
      do {
        onn = nn;
        nn = __sync_val_compare_and_swap(&position->next_, nn, last);
      } while(nn != onn);
      return nn;
    } else {
      nn              = position->next_;
      position->next_ = last;
      return nn;
    }
  }
  static base_node* erase_after_(iterator position, iterator last) {
    return erase_after_(position.head_, last.head_);
  }

  head_node head_;
};

#endif /* _LF_FORWARD_SIZE_LIST_HPP_ */
