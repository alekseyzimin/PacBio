#ifndef __MULTIPLEXER_H__
#define __MULTIPLEXER_H__

#include <iostream>
#include <sstream>

//#include <pthread.h>
#include <mutex>


class Multiplexer {
public:
  Multiplexer(std::ostream& os, size_t min_buffer = 4096, size_t max_buffer = 4 * 4096) :
    os_(os),
    min_buffer_(std::max((size_t)1024, min_buffer)),
    max_buffer_(std::max(max_buffer, min_buffer))
  { }

  virtual ~Multiplexer() {
    os_.flush();
  }

  class ostream : public std::stringstream {
  public:
    ostream(Multiplexer& m) : m_(m) { }
    ostream(Multiplexer* m) : m_(*m) { }
    ~ostream() { do_flush(); }

    void end_record() {
      const size_t current_size = this->tellp();
      if(current_size >= m_.min_buffer_) {
        if(current_size < m_.max_buffer_) {
          std::unique_lock<std::mutex> lock(m_.mutex_, std::try_to_lock);
          if(lock.owns_lock())
            do_flush_();
        } else {
          do_flush();
        }
      }
    }

    void do_flush() {
      std::unique_lock<std::mutex> lock(m_.mutex_);
      do_flush_();
    }

  protected:
    void do_flush_() {
      if(this->tellp() == 0) return;
      m_.os_ << this->rdbuf();
      // m_.os_ << this->str();
      this->str("");
      this->seekg(0);
    }

    Multiplexer& m_;
  };
  friend class ostream;

private:
  std::ostream& os_;
  const size_t  min_buffer_, max_buffer_;
  std::mutex    mutex_;
};

Multiplexer::ostream& endr(Multiplexer::ostream& os) {
  os.end_record();
  return os;
}

#endif /* __MULTIPLEXER_H__ */
