#include <src_jf_aligner/super_read_name.hpp>

std::string super_read_name::name() const {
  std::string res;
  auto it = unitigs_.cbegin();
  if(it != unitigs_.cend()) {
    res += unitig_name(*it);
    for(++it; it != unitigs_.cend(); ++it)
      (res += '_') += unitig_name(*it);
  }
  return res;
}

void super_read_name::append(const super_read_name& rhs, size_t skip) {
  if(skip >= rhs.size()) return;
  const auto osize = unitigs_.size();
  unitigs_.resize(unitigs_.size() + rhs.size() - skip);
  std::copy(rhs.unitigs_.cbegin() + skip, rhs.unitigs_.cend(), unitigs_.begin() + osize);
}

size_t super_read_name::prepend(const super_read_name& rhs, size_t all_but, size_t offset) {
  offset = std::min(offset, unitigs_.size());
  if(all_but >= rhs.size()) return offset;
  const size_t to_copy = rhs.size() - all_but;
  if(to_copy > offset) return offset;
  const size_t new_offset = offset - to_copy;
  std::copy_n(rhs.unitigs_.cbegin(), to_copy, unitigs_.begin() + new_offset);
  return new_offset;
}

void super_read_name::reverse() {
  const size_t n = unitigs_.size();
  for(size_t i = 0; i < unitigs_.size() / 2; ++i) {
    long tmp            = unitigs_[i];
    unitigs_[i]         = -unitigs_[n - i - 1];
    unitigs_[n - i - 1] = -tmp;
  }
  if(n % 2 == 1)
    unitigs_[n / 2] = -unitigs_[n / 2];
}

int super_read_name::overlap(const super_read_name& rhs) const {
  if(rhs.unitigs_.empty()) return 0;
  const long fu = rhs.unitigs_.front();

  for(auto it = std::find(unitigs_.cbegin(), unitigs_.cend(), fu);
      it != unitigs_.cend();
      it = std::find(it + 1, unitigs_.cend(), fu)) {
    const int olen = unitigs_.cend() - it;
    if((size_t)olen <= rhs.unitigs_.size() && std::equal(it, unitigs_.cend(), rhs.unitigs_.cbegin()))
      return olen;
  }
  return 0;
}

std::vector<long> super_read_name::parse(std::string name) {
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

std::ostream& operator<<(std::ostream& os, const super_read_name& sr) {
  auto it = sr.unitigs_.cbegin();
  if(it != sr.unitigs_.cend()) {
    os << abs(*it) << (*it > 0 ? 'F' : 'R');
    for(++it; it != sr.unitigs_.cend(); ++it)
      os << '_' << abs(*it) << (*it > 0 ? 'F' : 'R');
  }
  return os;
}
