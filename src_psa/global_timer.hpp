#ifndef __GLOBAL_TIMER_H__
#define __GLOBAL_TIMER_H__

#include <chrono>
#include <string>

class global_timer_type {
#ifdef SHOW_TIMING
  const char*                           m_msg;
  std::chrono::steady_clock::time_point m_start;
#endif

public:
#ifdef SHOW_TIMING
  global_timer_type() : m_msg(nullptr) { };
#else
  global_timer_type() = default;
#endif
  void start(const char* msg);
  void start(const std::string& msg) { return start(msg.c_str()); }
  void stop();

#ifdef SHOW_TIMING
  template<typename Barrier>
  void start(Barrier& b, const char* msg) {
    if(b.wait())
      start(msg);
    b.wait();
  }
  template<typename Barrier>
  void start(Barrier& b, const std::string& msg) { start(b, msg.c_str()); }

  template<typename Barrier>
  void stop(Barrier& b) {
    if(b.wait())
      stop();
    b.wait();
  }
#else
  template<typename Barrier>
  void start(Barrier& b, const char* msg) { }
  template<typename Barrier>
  void start(Barrier& b, const std::string& msg) { }
  template<typename Barrier>
  void stop(Barrier& b) { }
#endif
};

extern global_timer_type global_timer;

#endif /* __GLOBAL_TIMER_H__ */
