#ifndef __DEBUG_H__
#define __DEBUG_H__

#ifndef NDEBUG

#include <sstream>
#include <iostream>
#include <vector>

// Convenience: print vectors
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
  os << '[';
  auto it = v.cbegin();
  if(it != v.cend()) {
    os << *it;
    for(++it; it != v.cend(); ++it)
      os << ", " << *it;
  }
  return os << ']';
}

// Debugging class: output message on destructor. Behave like an
// output string stream.
struct dbg {
  std::ostringstream msg_;

  dbg() = default;
  virtual ~dbg() {
    std::cerr << msg_.str() << std::endl;
  }

  dbg& operator<<(dbg& (*pf)(dbg&)) {
    return pf(*this);
  }

  template<typename T>
  dbg& operator<<(const T& x) {
    msg_ << x;
    return *this;
  }
};

// Convenience macros
#define DBG if(1) dbg() << __FILE__ << ':' << __LINE__ << ':' << __FUNCTION__
#define W(m,x) xspace << m << ':' << (x)
#define V(x) xspace << #x << ':' << (x)

// Manipulators for dbg class. xspace add a space if needed
inline std::ostream& xspace(std::ostream& os) {
return os << ' ';
}

inline dbg& xspace(dbg& os) {
  if(os.msg_.tellp() > 0)
    os << ' ';
  return os;
}

#else
// NDEBUG defined
#define V(x)
#define W(m, x)
struct dbg {
  template<typename T>
  dbg& operator<<(const T& x) { return *this; }
};

#endif

#endif /* __DEBUG_H__ */
