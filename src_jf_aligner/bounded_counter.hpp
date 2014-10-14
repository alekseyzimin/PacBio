#ifndef __BOUNDED_COUNTER_H__
#define __BOUNDED_COUNTER_H__

#include <limits>
#include <atomic>


template<bool LF, typename T>
struct bounded_counter;

template<typename T>
struct bounded_counter<false, T> {
  T counter_;
  bounded_counter() : counter_(0) { }

  bool inc(const T max = std::numeric_limits<T>::max()) {
    if(counter_ < max) {
      ++counter_;
      return true;
    }
    return false;
  }

  operator T() const noexcept { return counter_; }
};

template<typename T>
struct bounded_counter<true, T> {
  std::atomic<T> counter_;
  bounded_counter() : counter_(0) { }

  bool inc(const T max = std::numeric_limits<T>::max()) {
    T val = counter_.load();
    while(!max || val < max) {
      if(counter_.compare_exchange_weak(val, val + 1))
        return true;
    }
    return false;
  }

  operator T() const noexcept { return counter_.load(); }
};

#endif /* __BOUNDED_COUNTER_H__ */
