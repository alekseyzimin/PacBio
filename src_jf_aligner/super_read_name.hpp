#ifndef __SUPER_READ_NAME_H__
#define __SUPER_READ_NAME_H__

#include <limits>
#include <algorithm>
#include <debug.hpp>

class super_read_name {
  std::vector<long> unitigs_;

public:
  static const unsigned long invalid = std::numeric_limits<unsigned long>::max();

  super_read_name() = default;
  explicit super_read_name(size_t n) : unitigs_(n) { }
  explicit super_read_name(const std::string& name) : unitigs_(parse(name)) { }
  explicit super_read_name(std::vector<long>&& unitigs) : unitigs_(std::move(unitigs)) { }
  // super_read_name(super_read_name&& rhs) : unitigs_(std::move(rhs.unitigs_)) { }
  // super_read_name(const super_read_name& rhs) : unitigs_(rhs.unitigs_) { }

  struct unitig {
    unsigned long id;
    char          ori; // 'F' or 'R'
  };

  const std::vector<long>& unitigs() const { return unitigs_; }
  size_t nb_unitigs() const { return unitigs_.size(); }
  size_t size() const { return nb_unitigs(); }
  std::string unitig_name(long i) const { return std::to_string(abs(i)) + (i > 0 ? 'F' : 'R'); }
  std::string name() const;

  void append(const super_read_name& rhs, size_t skip = 0);

  // Prepend all but 'all_but' elements from the 'rhs' into 'this'
  // starting at position 'offset'. Returns 'offset' minus the number
  // of elements copied. If 'rhs' is too big, nothing is copied.
  size_t prepend(const super_read_name& rhs, size_t all_but = 0, size_t offset = std::numeric_limits<size_t>::max());

  unitig operator[](size_t i) const {
    if(i < nb_unitigs())
      return { (unsigned long)abs(unitigs_[i]), unitigs_[i] > 0 ? 'F' : 'R' };
    return { invalid, '\0' };
  }

  unsigned long unitig_id(size_t i) const { return i < nb_unitigs() ? abs(unitigs_[i]) : invalid; }

  // static char reverse_ori(char ori) {
  //   switch(ori) {
  //   case 'F': return 'R';
  //   case 'R': return 'F';
  //   default: return ori; // Hum?? Bad orientation to begin with.
  //   }
  // }

  void reverse();

  super_read_name get_reverse() const {
    super_read_name res(*this);
    res.reverse();
    return res;
  }

  // Return the length of the longest overlap by unitigs between two
  // super reads, in a dovetail fashion. The return value is the
  // largest integer m such that the last m k-unitigs of *this are
  // equal to the first m k-unitigs of rhs. Note that this
  // relationship is NOT symmetrical, i.e. if *this overlaps with rhs
  // then it comes before rhs. If the return is 0, there is no
  // overlap.
  int overlap(const super_read_name& rhs) const;

protected:
  static std::vector<long> parse(std::string name);

  friend std::ostream& operator<<(std::ostream& os, const super_read_name& sr);
};

std::ostream& operator<<(std::ostream& os, const super_read_name& sr);
#endif /* __SUPER_READ_NAME_H__ */
