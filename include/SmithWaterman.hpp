#include <sys/types.h>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <queue>
#include <string>
#include <vector>

namespace std {
class thread;
}
struct FastaSequence {
  std::string description;
  std::string sequence;
};

class SmithWaterman {
 public:
  SmithWaterman(const std::string& query_seq_path,
                const std::string& target_seq_path);

  std::vector<size_t> solve();

  int validate(const std::string& ref_path);

  void report() const;

  static constexpr int64_t match_score = 5;      // Do not modify.
  static constexpr int64_t mismatch_score = -3;  // Do not modify.
  static constexpr int64_t gap_score = -4;       // Do not modify.
  // The highest scores
  std::vector<size_t> max_scores;

 private:

  // The query/target sequences
  std::vector<FastaSequence> query_seqs;
  std::vector<FastaSequence> target_seqs;

  // The reference scores.
  std::vector<size_t> refs;

  std::queue<std::pair<size_t,size_t> > problems;
  std::queue<std::thread*> threads;

  // Number of query/target sequences
  size_t query_seqs_size;
  size_t target_seqs_size;

  // Length of each sequence
  //std::vector<std::vector<int64_t>> query_seqs_lens;
  //std::vector<std::vector<int64_t>> target_seqs_lens;

  // The indices of the highest score for each query-target pair.
  //std::vector<size_t> max_positions;

  void read_seq(const std::string& seq_path, std::vector<FastaSequence>& seqs);

  void read_ref(const std::string& ref_path, std::vector<size_t>& refs);
};
void pair_align(FastaSequence* query_seq, FastaSequence* target_seq, size_t* sw);