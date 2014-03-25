#ifndef __SUPER_READ_NAME_H__
#define __SUPER_READ_NAME_H__

#include <limits>

class super_read_name {
  const std::string         name_;
  const std::vector<size_t> unitigs_;

public:
  static const unsigned long invalid = std::numeric_limits<unsigned long>::max();

  super_read_name(const std::string& name) : name_(name), unitigs_(parse(name)) { }

  struct unitig {
    unsigned long id;
    char          ori; // 'F' or 'R'
  };

  size_t nb_unitigs() const { return unitigs_.size() - 1; }
  size_t size() const { return nb_unitigs(); }
  const std::string& name() const { return name_; }

  unitig operator[](size_t i) const {
    if(i >= nb_unitigs()) return { std::numeric_limits<unsigned long>::max(), '\0' };
    return { std::stoul(name_.substr(unitigs_[i], unitigs_[i + 1] - unitigs_[i] - 2)), name_[unitigs_[i + 1] - 2] };
  }

  unsigned long unitig_id(size_t i) const {
    if(i < nb_unitigs())
      return std::stoul(name_.substr(unitigs_[i], unitigs_[i + 1] - unitigs_[i] - 2));
    return std::numeric_limits<unsigned long>::max();
  }

  static char reverse_ori(char ori) {
    switch(ori) {
    case 'F': return 'R';
    case 'R': return 'F';
    default: return ori; // Hum?? Bad orientation to begin with.
    }
  }

  std::string reverse() const {
    std::string res;
    if(nb_unitigs() > 0) {
      for(size_t i = nb_unitigs(); i >= 1; --i) {
        res += name_.substr(unitigs_[i - 1], unitigs_[i] - unitigs_[i - 1] - 2);
        res += reverse_ori(name_[unitigs_[i] - 2]);
        if(i > 1) res += '_';
      }
    }
    return res;
  }

protected:
  static std::vector<size_t> parse(std::string name) {
    std::vector<size_t> res;
    res.push_back(0);
    if(!name.empty()) {
      for(size_t n = name.find_first_of('_');  n != std::string::npos; n = name.find_first_of('_', n + 1))
        res.push_back(n + 1);
      res.push_back(name.size() + 1);
    }
    return res;
  }
};


#endif /* __SUPER_READ_NAME_H__ */
