#ifndef _SUPERREAD_PARSER_HPP_
#define _SUPERREAD_PARSER_HPP_

#include <vector>
#include <jellyfish/mer_dna.hpp>

#include <src_jf_aligner/jf_aligner.hpp>
#include <src_psa/compact_dna.hpp>
#include <src_psa/compact_index.hpp>
typedef int32_t saint_t;
#include <src_psa/psa.hpp>


// Bare minimum struct to behave like a pointer into a mer_dna, as far
// as mer_sa_imp::str_to_mer and SA::search are concerned.
template<typename T>
struct mer_dna_ptr {
  typedef typename jellyfish::mer_dna_ns::mer_base<T> base_type;

  const T& m_m;
  unsigned int m_off;

  explicit mer_dna_ptr(const T& m, unsigned int off = 0) : m_m(m), m_off(off) { }
  mer_dna_ptr(const mer_dna_ptr& rhs) : m_m(rhs.m_m), m_off(rhs.m_off) { }
  mer_dna_ptr& operator+=(unsigned int x) {
    m_off += x;
    return *this;
  }
  mer_dna_ptr operator+(unsigned int x) const {
    mer_dna_ptr res(*this);
    return res += x;
  }

  uint64_t str_to_mer(unsigned int mer_size) const {
    return m_m.get_bits(2 * (m_m.k() - mer_size - m_off), 2 * mer_size);
  }
};

namespace mer_sa_imp {
template<>
inline uint64_t str_to_mer<mer_dna_ptr<jellyfish::mer_dna_ns::mer_base_static<uint64_t, 0>>>(mer_dna_ptr<jellyfish::mer_dna_ns::mer_base_static<uint64_t, 0>> x, unsigned int mer_size) {
  return x.str_to_mer(mer_size);
}
template<>
inline uint64_t str_to_mer<mer_dna_ptr<jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>>>(mer_dna_ptr<jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>> x, unsigned int mer_size) {
  return x.str_to_mer(mer_size);
}
}

struct sequence_psa {
  struct offset_type {
    size_t header;
    size_t sequence;
  };
  typedef std::vector<offset_type> offsets_type;
  typedef mer_sa_imp::SA<compact_dna::const_iterator, compact_index<uint64_t>::iterator> SA;

  struct it_elt {
    const frag_lists::frag_info* frag;
    int                          offset;
  };
  static const int m_search_bits = 20;


  class pos_iterator : public std::iterator<std::input_iterator_tag, const it_elt> {
    const sequence_psa* m_psa;
    size_t              m_fwd_index, m_fwd_end;
    size_t              m_bwd_index, m_bwd_end;
    it_elt              m_elt;
    size_t              m_len;

  public:
    pos_iterator() : m_psa(nullptr), m_fwd_index(0), m_bwd_index(0) { }
    pos_iterator(const sequence_psa* psa, size_t fwd_index, size_t fwd_end, size_t bwd_index, size_t bwd_end, size_t len)
      : m_psa(psa)
      , m_fwd_index(fwd_index)
      , m_fwd_end(fwd_end)
      , m_bwd_index(bwd_index)
      , m_bwd_end(bwd_end)
      , m_len(len)
    { ++*this; }
    pos_iterator(const pos_iterator& rhs)
      : m_psa(rhs.m_psa)
      , m_fwd_index(rhs.m_fwd_index)
      , m_fwd_end(rhs.m_fwd_end)
      , m_bwd_index(rhs.m_bwd_index)
      , m_bwd_end(rhs.m_bwd_end)
      , m_elt(rhs.m_elt)
      , m_len(rhs.m_len)
    { }
    pos_iterator& operator=(const pos_iterator& rhs) {
      m_psa       = rhs.m_psa;
      m_fwd_index = rhs.m_fwd_index;
      m_fwd_end   = rhs.m_fwd_end;
      m_bwd_index = rhs.m_bwd_index;
      m_bwd_end   = rhs.m_bwd_end;
      m_elt       = rhs.m_elt;
      m_len       = rhs.m_len;
      return *this;
    }
    bool operator==(const pos_iterator& rhs) const {
      return m_fwd_index == rhs.m_fwd_index && m_bwd_index == rhs.m_bwd_index && m_psa == rhs.m_psa;
    }
    bool operator!=(const pos_iterator& rhs) const { return !operator==(rhs); }
    const it_elt& operator*() const { return m_elt; }
    const it_elt* operator->() const { return &m_elt; }
    pos_iterator& operator++() {
      bool fwd;
      while((fwd = m_fwd_index != m_fwd_end) || m_bwd_index != m_bwd_end) {
        const size_t x = (*m_psa->m_sa)[fwd ? m_fwd_index++ : m_bwd_index++];
        const size_t search = x >> m_search_bits;
        if(search >= m_psa->m_header_search.size()) continue;
        const auto start = m_psa->m_offsets.cbegin() + m_psa->m_header_search[search];
        const bool at_end = search >= m_psa->m_header_search.size() - 1;
        const auto end =
          at_end || m_psa->m_header_search[search + 1] >= m_psa->m_offsets.size() - 1
          ? m_psa->m_offsets.cend()
          : m_psa->m_offsets.cbegin() + m_psa->m_header_search[search + 1] + 1;
        // assert(start->sequence <= x);
        // assert(end == m_psa->m_offsets.cend() || end->sequence > x);
        const auto next  =
          start == end ? start : std::lower_bound(start, end, x, [](const offset_type& j, size_t i) { return j.sequence <= i; });
        const auto limit =
          next == m_psa->m_offsets.cend() ? m_psa->sequence_size() : next->sequence;
        if(x + m_len > limit) continue;
        const auto res = next - 1;
        m_elt.frag     = &m_psa->m_headers[res->header];
        m_elt.offset   = x - res->sequence + 1;
        if(!fwd) m_elt.offset = -m_elt.offset;
        return *this;
      }
      // Reach the end
      m_psa   = 0;
      m_fwd_index = m_fwd_end = 0;
      m_bwd_index = m_bwd_end = 0;
      return *this;
    }
    pos_iterator operator++(int) { pos_iterator res(*this); ++*this; return res; }
  };

  std::vector<frag_lists::frag_info> m_headers;
  std::vector<uint64_t>              m_sequence;
  std::vector<size_t>                m_header_search;
  offsets_type                       m_offsets;
  std::unique_ptr<PSA<compact_dna::const_iterator> > m_sa;

  sequence_psa() {
    m_offsets.push_back({ (size_t)0, (size_t)0 });
  }

  void append_fasta(std::istream& is);
  void append_fasta(const char* path) {
    std::ifstream       is(path);
    if(!is.good())
      throw std::runtime_error(std::string("Can't open file ") + path);
    append_fasta(is);
  }

  size_t sequence_size() const { return m_offsets.back().sequence; }
  size_t nb_mers() const { return sequence_size() - (m_sa->min_size() - 1) * m_headers.size(); }
  size_t nb_sequences() const { return m_headers.size(); }

  void compute_psa(unsigned int min_size, unsigned int max_size,
                   unsigned int threads = std::thread::hardware_concurrency()) {
    std::cerr << "compute_psa " << m_headers.size() << ' ' << sequence_size() << '\n';
    m_sa.reset(new PSA<compact_dna::const_iterator>(compact_dna::const_iterator_at(m_sequence.data()),
                                                    sequence_size(), min_size, max_size, threads));
  }
  bool check_suffixes(std::ostream& out = std::cout) const;
  bool check_psa() const {
    return m_sa->check() && m_sa->check_suffixes();
  }

  template<typename MerType>
  std::pair<pos_iterator, pos_iterator> equal_range(const MerType& m) const {
    return equal_range(m, m.get_reverse_complement());
  }

  template<typename MerType>
  std::pair<pos_iterator, pos_iterator> equal_range(const MerType& m, const MerType& rm,
                                                    size_t max = std::numeric_limits<size_t>::max()) const {
    auto fwd_res = m_sa->search(mer_dna_ptr<MerType>(m), m.k());
    if(fwd_res.first >= max)
      return std::make_pair(pos_iterator(), pos_iterator());
    auto bwd_res = m_sa->search(mer_dna_ptr<MerType>(rm), rm.k());
    if(fwd_res.first + bwd_res.first >= max)
      return std::make_pair(pos_iterator(), pos_iterator());
    if(rm < m) std::swap(fwd_res, bwd_res);
    return std::make_pair(pos_iterator(this,
                                       fwd_res.second, fwd_res.second + fwd_res.first,
                                       bwd_res.second, bwd_res.second + bwd_res.first,
                                       m.k()),
                          pos_iterator());
  }

  template<typename MerType>
  std::pair<pos_iterator, size_t> find_pos_size(const MerType& m, const MerType& rm) const {
    auto fwd_res = m_sa->search(mer_dna_ptr<MerType>(m), m.k());
    auto bwd_res = m_sa->search(mer_dna_ptr<MerType>(rm), rm.k());
    if(rm < m) std::swap(fwd_res, bwd_res);
    return std::make_pair(pos_iterator(this,
                                       fwd_res.second, fwd_res.second + fwd_res.first,
                                       bwd_res.second, bwd_res.second + bwd_res.first,
                                       m.k()),
                          fwd_res.first + bwd_res.first);
  }

  const pos_iterator pos_end() const { return pos_iterator(); }

  template<typename Iterator>
  static sequence_psa read_files(Iterator start, Iterator end) {
    sequence_psa res;
    for( ; start != end; ++start)
      res.append_fasta(*start);
    return res;
  }

  template<typename File>
  static sequence_psa read_file(File& path) {
    sequence_psa res;
    res.append_fasta(path);
    return res;
  }
};

template<typename File>
sequence_psa superread_parse(File& file, unsigned int min_size = 10, unsigned int max_size = 17) {
  auto res = sequence_psa::read_file(file);
  res.compute_psa(min_size, max_size);
  return res;
}

template<typename FileIterator>
sequence_psa superread_parse(FileIterator first, FileIterator last, unsigned int min_size = 10, unsigned int max_size = 17) {
  auto res = sequence_psa::read_files(first, last);
  res.compute_psa(min_size, max_size);
  return res;
}


#endif /* _SUPERREAD_PARSER_HPP_ */
