/******************************************
Copyright University of Maryland 2015
******************************************/
#ifndef __BARRIER_H__
#define __BARRIER_H__

#include <thread>
#include <mutex>
#include <condition_variable>
#include <stdexcept>

template<typename MutexT>
class barrier {
public:
  barrier(unsigned int count)
    : m_threshold(count)
    , m_count(count)
    , m_generation(0)
  {
    if (count == 0)
      throw std::invalid_argument("count cannot be zero.");
  }

  bool wait() {
    std::unique_lock<MutexT> lock(m_mutex);
    const unsigned int gen = m_generation;

    if(--m_count == 0) {
      ++m_generation;
      m_count = m_threshold;
      m_cond.notify_all();
      return true;
    }

    while(gen == m_generation)
      m_cond.wait(lock);

    return false;
  }

private:
  MutexT                  m_mutex;
  std::condition_variable m_cond;
  const unsigned int      m_threshold;
  volatile unsigned int   m_count;
  volatile unsigned int   m_generation;
};

#endif /* __BARRIER_H__ */
