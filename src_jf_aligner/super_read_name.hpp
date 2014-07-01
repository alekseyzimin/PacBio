#ifndef __SUPER_READ_NAME_H__
#define __SUPER_READ_NAME_H__

#include <limits>
#include <algorithm>
#include <debug.hpp>

class super_read_name {
  const std::vector<long> unitigs_;

public:
  static const unsigned long invalid = std::numeric_limits<unsigned long>::max();

  super_read_name(const std::string& name) : unitigs_(parse(name)) { }
  super_read_name(std::vector<long>&& unitigs) : unitigs_(std::move(unitigs)) { }

  struct unitig {
    unsigned long id;
    char          ori; // 'F' or 'R'
  };

  size_t nb_unitigs() const { return unitigs_.size(); }
  size_t size() const { return nb_unitigs(); }
  std::string unitig_name(long i) const { return std::to_string(abs(i)) + (i > 0 ? 'F' : 'R'); }
  std::string name() const {
    std::string res;
    auto it = unitigs_.cbegin();
    if(it != unitigs_.cend()) {
      res += unitig_name(*it);
      for(++it; it != unitigs_.cend(); ++it)
        (res += '_') += unitig_name(*it);
    }
    return res;
  }

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

  super_read_name reverse() const {
    const size_t n = unitigs_.size();
    std::vector<long> us(n, 0);
    for(size_t i = 0; i < us.size() / 2; ++i) {
      us[i]         = -unitigs_[n - i - 1];
      us[n - i - 1] = -unitigs_[i];
    }
    if(n % 2 == 1)
      us[n / 2] = -unitigs_[n / 2];
    return super_read_name(std::move(us));
  }

  // Return the length of the longest overlap by unitigs between two
  // super reads, in a dovetail fashion. The return value is the
  // largest integer m such that the last m k-unitigs of *this are
  // equal to the first m k-unitigs of rhs. Note that this
  // relationship is NOT symmetrical, i.e. if *this overlaps with rhs
  // then it comes before rhs. If the return is 0, there is no
  // overlap.
  int overlap(const super_read_name& rhs) const {
    if(rhs.unitigs_.empty()) return 0;
    const long fu = rhs.unitigs_.front();

    for(auto it = std::find(unitigs_.cbegin(), unitigs_.cend(), fu);
        it != unitigs_.cend();
        it = std::find(it + 1, unitigs_.cend(), fu)) {
      const int olen = unitigs_.cend() - it;
      if(olen <= rhs.unitigs_.size() && std::equal(it, unitigs_.cend(), rhs.unitigs_.cbegin()))
        return olen;
    }
    return 0;
  }

protected:
  static std::vector<long> parse(std::string name) {
    std::vector<long> res;
    if(!name.empty()) {
      size_t pn = 0;
      for(size_t n = name.find_first_of('_'); n != std::string::npos; pn = n + 1, n = name.find_first_of('_', pn)) {
        long id = std::stoul(name.c_str() + pn);
        res.push_back(name[n - 1] == 'R' ? -id : id);
      }
      long id = std::stoul(name.c_str() + pn);
      res.push_back(name[name.size() - 1] == 'R' ? -id : id);
    }
    return res;
  }
};


#endif /* __SUPER_READ_NAME_H__ */
