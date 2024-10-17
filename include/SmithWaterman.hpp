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

  std::vector<uint32_t> solve();

  int validate(const std::string& ref_path);

  void report() const;

  static constexpr int32_t match_score = 5;      // Do not modify.
  static constexpr int32_t mismatch_score = -3;  // Do not modify.
  static constexpr int32_t gap_score = -4;       // Do not modify.
  // The highest scores
  std::vector<uint32_t> max_scores;

 private:

  // The query/target sequences
  std::vector<FastaSequence> query_seqs;
  std::vector<FastaSequence> target_seqs;

  // The reference scores.
  std::vector<uint32_t> refs;

  std::vector<std::pair<uint32_t,uint32_t> > problems;

  // Number of query/target sequences
  uint32_t query_seqs_size;
  uint32_t target_seqs_size;

  // Length of each sequence
  //std::vector<std::vector<int32_t>> query_seqs_lens;
  //std::vector<std::vector<int32_t>> target_seqs_lens;

  // The indices of the highest score for each query-target pair.
  //std::vector<uint32_t> max_positions;

  void read_seq(const std::string& seq_path, std::vector<FastaSequence>& seqs);

  void read_ref(const std::string& ref_path, std::vector<uint32_t>& refs);
};
void pair_align(FastaSequence* query_seq, FastaSequence* target_seq, uint32_t* sw);