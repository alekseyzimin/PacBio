/******************************************
Copyright University of Maryland 2015
******************************************/
#ifndef __BOUNDED_COUNTER_H__
#define __BOUNDED_COUNTER_H__

#include <limits>
#include <atomic>
#include <type_traits>


template<bool LF, typename T>
struct bounded_counter;

template<typename T>
struct bounded_counter<false, T> {
  T counter_;

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
  T counter_;

  bool inc(const T max = std::numeric_limits<T>::max()) {
    T oval = counter_;
    while(!max || oval < max) {
      T nval = __sync_val_compare_and_swap(&counter_, oval, oval + 1);
      if(nval == oval)
        return true;
      oval = nval;
    }
    return false;
  }

  operator T() const noexcept { return counter_; }
};

#endif /* __BOUNDED_COUNTER_H__ */
