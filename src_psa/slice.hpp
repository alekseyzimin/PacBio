/******************************************
Copyright University of Maryland 2015
******************************************/
#ifndef _SLICE_H_
#define _SLICE_H_

#include <thread>
#include <atomic>

template<typename T, typename B>
class slice_for {
  std::atomic<T> m_slice;
  B&             m_barrier;

public:
  slice_for(B& barrier) : m_barrier(barrier) { }

  template<typename R>
  void loop(T start, T end, T step, R block) {
    if(m_barrier.wait())
      m_slice = 0;
    m_barrier.wait();
    while(true) {
      const T slice_start = start + step * m_slice++;
      if(slice_start >= end) break;
      const T slice_end = std::min(slice_start + step, end);
      for(T i = slice_start; i < slice_end; ++i)
        block(i);
    }
  }

  template<typename R>
  void chunk(T start, T end, T step, R block) {
    if(m_barrier.wait())
      m_slice = 0;
    m_barrier.wait();
    while(true) {
      const T slice_start = start + step * m_slice++;
      if(slice_start >= end) break;
      const T slice_end = std::min(slice_start + step, end);
      block(slice_start, slice_end);
    }
  }
};

#endif /* _SLICE_H_ */
