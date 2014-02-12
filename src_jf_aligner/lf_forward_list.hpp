#ifndef _LF_FORWARD_LIST_HPP_
#define _LF_FORWARD_LIST_HPP_

#include <memory>

template<typename T, class Alloc = std::allocator<T>>
class lf_forward_list {
  struct head_node {
    head_node* next_;
  };
  struct node : head_node {
    T     val_;
  };

public:
  typedef Alloc allocator_type;
  typedef T     value_type;
  class iterator;
  class const_iterator;

  /** Construct an empty list */
  explicit lf_forward_list (const allocator_type& alloc = allocator_type()) :
    T_alloc_(alloc), node_alloc_(alloc), head_(0)
  { }

  ~lf_forward_list() {
    clear();
  }

  /** Remove all element from the list. Not safe if other threads are
      accessing elements in the list, but it can be mixed with other
      insert operations. */
  void clear() {
    head_node* hn = *const_cast<head_node* volatile*>(&head_);
    head_node* ohn;

    do {
      ohn = hn;
      hn  = __sync_val_compare_and_swap (&head_, hn, 0);
    } while(hn != ohn);

    while(hn) {
      ohn = hn;
      T_alloc_.destroy(&static_cast<node*>(hn)->val_);
      node_alloc_.deallocate(static_cast<node*>(hn), 1);
      hn = ohn->next_;
    }
  }

  /** Push an element to the front of the list (copy) */
  void push_front(const value_type& val) {
    node* n = node_alloc_.allocate(1);
    T_alloc_.construct(&n->val_, val);
    push_front_(n);
  }

  /** Push an element to the front of the list (move) */
  void push_front(value_type&& val) {
    node* n = node_alloc_.allocate(1);
    T_alloc_.construct(&n->val_, std::move(val));
    push_front_(n);
  }

  /** Get an iterator */
  iterator begin() { return iterator(head_); }
  const_iterator begin() const { return const_iterator(head_); }
  const_iterator cbegin() const { return const_iterator(head_); }
  iterator end() { return iterator(); }
  const_iterator end() const { return const_iterator(); }
  const_iterator cend() const { return const_iterator(); }

protected:
  void push_front_(head_node* n) {
    head_node* hn = *const_cast<head_node* volatile*>(&head_);
    do {
      n->next_ = hn;
      hn = __sync_val_compare_and_swap (&head_, hn, n);
    } while(hn != n->next_);
  }

public:
  class iterator : public std::iterator<std::forward_iterator_tag, T> {
  public:
    explicit iterator(head_node* head = 0) : head_(head) { }
    iterator(const iterator& rhs) : head_(rhs.head_) { }
    iterator& operator=(const iterator& rhs) {
      head_ = rhs.head_;
      return *this;
    }
    bool operator==(const iterator& rhs) { return head_ == rhs.head_; }
    bool operator!=(const iterator& rhs) { return head_ != rhs.head_; }
    T& operator*() { return static_cast<node*>(head_)->val_; }
    const T& operator*() const { return static_cast<const node*>(head_)->val_; }
    T* operator->() { return &static_cast<node*>(head_)->val_; }
    const T& operator->() const { return &static_cast<const node*>(head_)->val_; }
    iterator& operator++() { head_ = (*const_cast<head_node* volatile*>(&head_))->next_; return *this; }
    iterator operator++(int) {
      iterator res(*this);
      ++*this;
      return res;
    }
  private:
    head_node* head_;
  };

  class const_iterator : public std::iterator<std::forward_iterator_tag, T> {
  public:
    explicit const_iterator(head_node* head = 0) : head_(head) { }
    const_iterator(const const_iterator& rhs) : head_(rhs.head_) { }
    const_iterator(const iterator& rhs) : head_(rhs.head_) { }
    const_iterator& operator=(const const_iterator& rhs) {
      head_ = rhs.head_;
      return *this;
    }
    bool operator==(const const_iterator& rhs) const { return head_ == rhs.head_; }
    bool operator!=(const const_iterator& rhs) const { return head_ != rhs.head_; }
    const T& operator*() const { return static_cast<const node*>(head_)->val_; }
    const T& operator->() const { return &static_cast<const node*>(head_)->val_; }
    const_iterator& operator++() { head_ = (*const_cast<head_node* volatile*>(&head_))->next_; return *this; }
    const_iterator operator++(int) {
      const_iterator res(*this);
      ++*this;
      return res;
    }
  private:
    head_node* head_;
  };

private:
  typedef typename Alloc::template rebind<T>::other    T_alloc_type;
  typedef typename Alloc::template rebind<node>::other node_alloc_type;
  T_alloc_type    T_alloc_;
  node_alloc_type node_alloc_;
  head_node*      head_;
};

#endif /* _LF_FORWARD_LIST_HPP_ */
