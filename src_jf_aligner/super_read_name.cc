#include <src_jf_aligner/super_read_name.hpp>
#include <stdexcept>

const uint32_t super_read_name::invalid_id;

std::string super_read_name::name() const {
  std::string res;
  auto it = unitigs_.cbegin();
  if(it != unitigs_.cend()) {
    res += it->name();
    for(++it; it != unitigs_.cend(); ++it)
      (res += '_') += it->name();
  }
  return res;
}

void super_read_name::append(const super_read_name& rhs, size_t skip) {
  if(skip >= rhs.size()) return;
  const auto osize = unitigs_.size();
  unitigs_.resize(unitigs_.size() + rhs.size() - skip);
  std::copy(rhs.unitigs_.cbegin() + skip, rhs.unitigs_.cend(), unitigs_.begin() + osize);
}

size_t super_read_name::prepend(size_t offset, const super_read_name& rhs, size_t first, size_t last) {
  if(first > last || first >= rhs.size()) return offset;
  const size_t to_copy = std::min(last, rhs.size() - 1) - first + 1;
  if(to_copy > offset) return offset;
  const size_t new_offset = offset - to_copy;
  std::copy_n(rhs.unitigs_.cbegin() + first, to_copy, unitigs_.begin() + new_offset);
  return new_offset;
}

void super_read_name::reverse() {
  const size_t n = unitigs_.size();
  for(size_t i = 0; i < unitigs_.size() / 2; ++i) {
    auto tmp            = unitigs_[i];
    unitigs_[i]         = unitigs_[n - i - 1].reversed();
    unitigs_[n - i - 1] = tmp.reversed();
  }
  if(n % 2 == 1)
    unitigs_[n / 2].reverse();
}

int super_read_name::overlap(const super_read_name& rhs) const {
  if(rhs.unitigs_.empty()) return 0;
  const u_id_ori *plhs,*prhs;
  uint32_t slhs=unitigs_.size(),srhs=rhs.unitigs_.size(),start_offset;
  if(slhs<2 || srhs<2) return 0;
  plhs=unitigs_.data();
  prhs=rhs.unitigs_.data();
  start_offset=(slhs>=srhs)?(slhs-srhs+1):1;
  uint32_t res=0;
  for(uint32_t i=start_offset;i<slhs;i++){
	if(prhs[0]==plhs[i]){
	//found the candidate
	res=1;
	for(uint32_t j=i+1,k=1;j<slhs;j++,k++){
		if(plhs[j]==prhs[k]){
		res++;
		}else{
		res=0;
		break;
		}
	}
	if(res>0)
		return res;
	}
  }
  
/*  const auto fu = rhs.unitigs_.front();

  for(auto it = std::find(unitigs_.cbegin(), unitigs_.cend(), fu);
      it != unitigs_.cend();
      it = std::find(it + 1, unitigs_.cend(), fu)) {
    const int olen = unitigs_.cend() - it;
    if((size_t)olen <= rhs.unitigs_.size() && std::equal(it, unitigs_.cend(), rhs.unitigs_.cbegin()))
      return olen;
  }*/
  return 0;
}

super_read_name::unitigs_list super_read_name::parse(const std::string& name) {
  unitigs_list res;
  try {
    if(!name.empty()) {
      size_t pn = 0;
      for(size_t n = name.find_first_of('_'); n != std::string::npos; pn = n + 1, n = name.find_first_of('_', pn)) {
        uint32_t id = std::stoul(name.c_str() + pn);
        res.push_back(u_id_ori(id, name[n - 1]));
      }
      uint32_t id = std::stoul(name.c_str() + pn);
      res.push_back(u_id_ori(id, name[name.size() - 1]));
    }
  } catch(std::invalid_argument) {
    res.clear();
  }
  return res;
}

std::ostream& operator<<(std::ostream& os, const super_read_name& sr) {
  auto it = sr.unitigs_.cbegin();
  if(it != sr.unitigs_.cend()) {
    os << it->id() << it->ori();
    for(++it; it != sr.unitigs_.cend(); ++it)
      os << '_' << it->id() << it->ori();
  }
  return os;
}

static char rev_comp_(char c) {
  switch(c) {
  case 'a': case 'A': return 'T';
  case 'c': case 'C': return 'G';
  case 'g': case 'G': return 'C';
  case 't': case 'T': return 'A';
  default: return 'N';
  }
}

static void print_unitig(std::ostream& os, bool rev_comp, const std::string& seq, int offset) {
  if((size_t)offset < seq.size()) {
    if(rev_comp) {
      for(auto it = seq.crbegin() + offset; it != seq.crend(); ++it)
        os << rev_comp_(*it);
    } else {
      os << (seq.c_str() + offset);
    }
  }
}

void super_read_name::print_sequence(std::ostream& os, const std::vector<std::string>& unitigs_sequences,
                                     int k_len, size_t start_unitig, ssize_t nb_unitigs) const {
  auto       it  = unitigs_.cbegin() + std::min(start_unitig, unitigs_.size());
  const auto end = (nb_unitigs == -1) ? unitigs_.cend() : unitigs_.cbegin() + std::min(start_unitig + nb_unitigs, unitigs_.size());
  if(it != end) {
    print_unitig(os, it->ori_, unitigs_sequences.at(it->id()), 0);
    for(++it; it != end; ++it)
      print_unitig(os, it->ori_, unitigs_sequences.at(it->id()), k_len - 1);
  }
}
