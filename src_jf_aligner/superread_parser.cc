#include <thread>

#include <src_jf_aligner/superread_parser.hpp>
#include <src_psa/barrier.hpp>
#include <src_psa/slice.hpp>
#include <src_psa/compact_index.hpp>
#include <src_psa/mer_sa_imp.hpp>

void sequence_psa::append_fasta(std::istream& is) {
  int c = is.peek();
  if(is.peek() != '>')
    throw std::runtime_error("Not in fasta format");

  size_t       seq_offset    = m_offsets.back().sequence;
  size_t       search_offset = m_header_search.back();
  std::string  line;
  std::string  header_line;
  for( ; c != EOF; c = is.peek()) {
    // read header
    std::getline(is, header_line);

    // read sequence
    size_t old_seq_offset = seq_offset;
    for(c = is.peek(); c != '>' && c != EOF; c = is.peek()) {
      std::getline(is, line);
      if(m_sequence.size() * sizeof(uint64_t) * 4 < line.size() + seq_offset)
        m_sequence.resize(std::max(1 + (size_t)(line.size() + seq_offset) / (sizeof(uint64_t) * 4), m_sequence.size() * 2));
      compact_dna::copy_from_str(compact_dna::iterator(m_sequence.data(), 2, 0) + seq_offset, line);
      // TODO: Should we check that everything was copied?
      seq_offset += line.size();
    }

    if(seq_offset > old_seq_offset) {
      m_headers.push_back(frag_lists::frag_info(seq_offset - old_seq_offset, header_line.c_str() + 1));
      const size_t next_search_offset = seq_offset >> m_search_bits;
      for(size_t i = search_offset + 1; i <= next_search_offset; ++i)
        m_header_search.push_back(m_offsets.size() - 1);
      search_offset = next_search_offset;
      m_offsets.push_back({ m_headers.size(), seq_offset });
    }
  }
}

static void
create_thread(int thid, int nb_threads, barrier<std::mutex>* thread_barrier, slice_for<size_t, barrier<std::mutex>>* slicer,
              compact_dna::const_iterator T, compact_index<uint64_t>::iterator SA,
              const sequence_psa::offsets_type& offsets,
              size_t n, uint64_t* mer_counts, unsigned int mer_size, unsigned int max_size) {
  const size_t nb_counts    = (1 << (2 * mer_size)) + 1;
  const size_t counts_mask  = ~(size_t)0 >> (8 * sizeof(size_t) - 2 * mer_size);
  const size_t counts_start = (nb_counts - 1) * thid / nb_threads;
  const size_t counts_end   = (nb_counts - 1) * (thid + 1) / nb_threads;
  //  const size_t n_step       = std::min((size_t)(1024 * 1024), std::max((size_t)1, n / nb_threads));

  // Add 1 to length of zeroing if thid == nb_threads because
  // counts_end is at most nb_counts - 1, and need to zero up to
  // nb_counts.
  std::fill_n(mer_counts + counts_start, counts_end - counts_start + (thid == nb_threads), (uint64_t)0);
  thread_barrier->wait();

  {
    std::unique_ptr<uint64_t[]> tmp_counts(new uint64_t[nb_counts - 1]);
    std::fill_n(tmp_counts.get(), nb_counts - 1, (uint64_t)0);
    slicer->loop(0, offsets.size() - 1, 10, [=](size_t i) {
        sequence_psa::SA::count_mers(T + offsets[i].sequence,
                                     offsets[i + 1].sequence - offsets[i].sequence,
                                     mer_counts, mer_size);
      });
    for(size_t i = 0, j = counts_start; i < nb_counts - 1; ++i, j = ((j + 1) & counts_mask))
      __sync_fetch_and_add(mer_counts + j, tmp_counts[j]);
  }

  if(thread_barrier->wait())
    sequence_psa::SA::partial_sums(mer_counts, nb_counts);

  typedef typename compactsufsort_imp::parallel_iterator_traits<compact_index<uint64_t>::iterator>::type PSAIDPTR;
  PSAIDPTR PSA(SA);
  slicer->loop(0, offsets.size() - 1, 10, [=](size_t i) {
    sequence_psa::SA::fill_mers(T, offsets[i].sequence, offsets[i + 1].sequence, PSA, mer_counts, mer_size);
    });

  slicer->loop(0, nb_counts - 1, 100, [=](size_t i) {
      sequence_psa::SA::sort_one_mer(T, n, PSA, mer_counts[i], mer_counts[i + 1], mer_size, max_size);
    });
}


void sequence_psa::compute_psa(unsigned int min_size, unsigned int max_size, unsigned int threads) {
  m_min_size = min_size;
  m_max_size = max_size;

  barrier<std::mutex>      SA_barrier(threads);
  slice_for<size_t, barrier<std::mutex>> slicer(SA_barrier);
  std::vector<std::thread> thread_handles;
  m_counts.resize((1 << (2 * min_size)) + 1);
  m_sa.reset(new compact_index<uint64_t>(nb_mers(), compact_index<uint64_t>::required_bits(sequence_size())));

  for(unsigned int i = 0; i < threads; ++i)
    thread_handles.push_back(std::thread(create_thread, i, threads, &SA_barrier, &slicer,
                                         compact_dna::const_iterator_at(m_sequence.data()),
                                         m_sa->begin(), m_offsets,
                                         sequence_size(), m_counts.data(), min_size, max_size));
  for(unsigned int i = 0; i < threads; ++i)
    thread_handles[i].join();
}

// void superread_parse(int threads, mer_pos_hash_type& hash, frag_lists& names,
//                      file_vector::const_iterator begin, file_vector::const_iterator end,
//                      bool compress) {
//   names.ensure(threads);
//   stream_manager streams(begin, end);
//   superreads_read_mers reader(threads, &hash, 0, names, streams, compress);
//   reader.exec_join(threads);
// }

// void superread_parse(int threads, mer_pos_hash_type& hash, frag_lists& names, const char* file, bool compress) {
//   file_vector files;
//   files.push_back(file);
//   superread_parse(threads, hash, names, files.cbegin(), files.cend(), compress);
// }
