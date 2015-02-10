#include <src_psa/global_timer.hpp>

#ifdef SHOW_TIMING
#include <iostream>

double seconds(std::chrono::steady_clock::time_point start, std::chrono::steady_clock::time_point end) {
  return std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
}

void global_timer_type::start(const char* msg) {
  stop();
  m_msg = msg;
  std::cerr << "Starting " << msg << " ..." << std::flush;
  m_start = std::chrono::steady_clock::now();
}

void global_timer_type::stop() {
  auto end = std::chrono::steady_clock::now();
  if(m_msg)
    std::cerr << ' ' << seconds(m_start, end) << std::endl;
  m_msg = nullptr;
}

#else // don't show timing
void global_timer_type::start(const char* msg) { }
void global_timer_type::stop() { }
#endif

global_timer_type global_timer;
