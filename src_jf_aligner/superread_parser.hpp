#ifndef _SUPERREAD_PARSER_HPP_
#define _SUPERREAD_PARSER_HPP_

#include <vector>
#include <jellyfish/mer_dna.hpp>

#include <src_jf_aligner/jf_aligner.hpp>
#include <src_psa/compact_dna.hpp>
#include <src_psa/compact_index.hpp>
typedef int32_t saint_t;
#include <src_psa/mer_sa_imp.hpp>

// jellyfish mer_dna are already compact dna.
// struct mer_dna_off {
//   const mer_dna& m_m;
//   unsigned int m_off;

//   mer_dna_off(const mer_dna& m) : m_m(m), m_off(0) { }
//   mer_dna_off(const mer_dna_off& rhs) : m_m(rhs.m_m), m_off(rhs.m_off) { }
//   mer_dna_off& operator+=(unsigned int x) {
//     m_off += x;
//     return *this;
//   }
//   mer_dna_off operator+(unsigned int x) const {
//     mer_dna_off res(*this);
//     return res += x;
//   }
// };
// template<>
// inline uint64_t mer_sa_imp::str_to_mer<mer_dna_off>(mer_dna_off x, unsigned int mer_size) {
//   return x.m_m.get_bits(x.m_off, mer_size);
// }

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
    size_t              m_index, m_end;
    it_elt              m_elt;
    size_t              m_len;

  public:
    pos_iterator() : m_psa(nullptr), m_index(0) { }
    pos_iterator(const sequence_psa& psa, size_t index, size_t end, size_t len)
      : m_psa(&psa)
      , m_index(index)
      , m_end(end)
      , m_len(len)
    { ++*this; }
    pos_iterator(const pos_iterator& rhs)
      : m_psa(rhs.m_psa)
      , m_index(rhs.m_index)
      , m_end(rhs.m_end)
      , m_elt(rhs.m_elt)
      , m_len(rhs.m_len)
    { }
    pos_iterator& operator=(const pos_iterator& rhs) {
      m_psa      = rhs.m_psa;
      m_index    = rhs.m_index;
      m_end      = rhs.m_end;
      m_elt      = rhs.m_elt;
      m_len      = rhs.m_len;
      return *this;
    }
    bool operator==(const pos_iterator& rhs) const { return m_index == rhs.m_index && m_psa == rhs.m_psa; }
    bool operator!=(const pos_iterator& rhs) const { return m_index != rhs.m_index || m_psa != rhs.m_psa; }
    const it_elt& operator*() const { return m_elt; }
    const it_elt* operator->() const { return &m_elt; }
    pos_iterator& operator++() {
      std::cout << "op++ m_index:" << m_index << " m_end:" << m_end << "\n";
      for( ; m_index != m_end; ++m_index) {
        const size_t x = (*m_psa->m_sa)[m_index];
        const size_t search = x >> m_search_bits;
        if(search >= m_psa->m_header_search.size()) continue;
        const auto start = m_psa->m_offsets.cbegin() + m_psa->m_header_search[search];
        const auto end = start + (search == m_psa->m_header_search.size() - 1
                                  ? m_psa->m_header_search.size() - 1
                                  : m_psa->m_header_search[search + 1] + 1);
        const auto next = std::lower_bound(start, end, x, [](const offset_type& j, size_t i) { return j.sequence <= i; });
        std::cout << "loop m_index:" << m_index << " m_end:" << m_end << " x:" << x << " m_len:" << m_len << " next:" << next->sequence << "\n";
        if(x + m_len > next->sequence) continue;
        const auto res = next - 1;
        m_elt.frag = &m_psa->m_headers[res->header];
        m_elt.offset = x - res->sequence;
        return *this;
      }
      // Reach the end
      m_psa   = 0;
      m_index = 0;
      return *this;
    }
    pos_iterator operator++(int) { pos_iterator res(*this); ++*this; return res; }
  };

  std::vector<frag_lists::frag_info>       m_headers;
  std::vector<uint64_t>                    m_sequence;
  std::vector<size_t>                      m_header_search;
  offsets_type                             m_offsets;
  std::vector<uint64_t>                    m_counts;
  std::unique_ptr<compact_index<uint64_t>> m_sa;
  unsigned int                             m_min_size, m_max_size;

  //  explicit sequence_psa(int search_bits = 20) : m_search_bits(search_bits) {
  sequence_psa() {
    m_offsets.push_back({ (size_t)0, (size_t)0 });
    m_header_search.push_back(0);
  }

  void append_fasta(std::istream& is);
  void append_fasta(const char* path) {
    std::ifstream       is(path);
    if(!is.good())
      throw std::runtime_error(std::string("Can't open file ") + path);
    append_fasta(is);
  }
  size_t sequence_size() const { return m_offsets.back().sequence; }
  size_t nb_mers() const { return sequence_size() - (m_min_size - 1) * m_headers.size(); }

  void compute_psa(unsigned int min_size, unsigned int max_size,
                   unsigned int threads = std::thread::hardware_concurrency());
  bool check_psa() const {
    std::cout << "check " << sequence_size() << ' ' << nb_mers() << std::endl;
    return SA::check(compact_dna::iterator_at(m_sequence.data()), sequence_size(),
                     m_sa->cbegin(), nb_mers(), m_counts.data(), m_min_size, m_max_size);
  }

  template<typename MerType>
  std::pair<pos_iterator, pos_iterator> equal_range(const MerType& m) {
    std::cout << "search for " << m.k() << ' ' << m << ' ' << std::hex << m.word(0) << std::endl;
    std::cout << mer_sa_imp::str_to_mer(m.to_str().c_str(), m.k()) << std::endl;
    auto res = SA::search(compact_dna::const_iterator_at(m_sequence.data()), sequence_size(),
                          m_sa->begin(), nb_mers(), m_counts.data(),
                          m_min_size, m_max_size,
                          m.to_str().c_str(), m.k());
                          //
                          //                          mer_dna_off(m), m.k());
    std::cout << res.first << ' ' << res.second << std::endl;
    return std::make_pair(pos_iterator(*this, res.second, res.second + res.first, m.k()), pos_iterator());
  }

  template<typename Iterator>
  static sequence_psa read_files(Iterator start, Iterator end) {
    sequence_psa res;
    for( ; start != end; ++start) {
      //      global_timer.start(std::string("Reading file ") + *start);
      res.append_fasta(*start);
    }
    //    global_timer.stop();
    return res;
  }
};

// void superread_parse(int threads, mer_pos_hash_type& hash, frag_lists& names,
//                      file_vector::const_iterator begin, file_vector::const_iterator end,
//                      bool compress = false);

// void superread_parse(int threads, mer_pos_hash_type& hash, frag_lists& names, const char* file, bool compress = false);

#endif /* _SUPERREAD_PARSER_HPP_ */
