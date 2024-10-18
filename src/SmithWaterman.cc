#include "SmithWaterman.hpp"

#include <immintrin.h>
#include <omp.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>

SmithWaterman::SmithWaterman(const std::string& query_seq_path,
                             const std::string& target_seq_path) {
  read_seq(query_seq_path, query_seqs);
  read_seq(target_seq_path, target_seqs);

  query_seqs_size = query_seqs.size();
  assert(query_seqs_size >= 1);
  target_seqs_size = target_seqs.size();
  assert(target_seqs_size >= 1);
}

std::vector<uint32_t> SmithWaterman::solve() {
  while (!problems.empty()) {
    problems.pop_back();
  }
  // Iterate through the query sequences
  for (uint32_t i = 0; i < query_seqs.size(); i++) {
    // Iterate throuth the target sequences
    for (uint32_t j = 0; j < target_seqs.size(); j++) {
      std::pair<uint32_t, uint32_t> problem = {i, j};
      problems.push_back(problem);
    }
  }

  max_scores.resize(query_seqs.size() * target_seqs.size());

  omp_set_num_threads(problems.size());
#pragma omp parallel
  {
    uint32_t thread_number = omp_get_thread_num();
    if (thread_number < problems.size()) {
      pair_align(&query_seqs[problems[thread_number].first],
                 &target_seqs[problems[thread_number].second],
                 max_scores.begin().base() + thread_number);
    }
  }
  return max_scores;
}

void SmithWaterman::report() const {
  for (uint32_t i = 0; i < query_seqs_size; i++) {
    auto max_score_ptr =
        std::max_element(max_scores.cbegin() + i * target_seqs_size,
                         max_scores.cbegin() + (i + 1) * target_seqs_size);
    uint32_t max_score_idx = std::distance(max_scores.cbegin(), max_score_ptr);
    std::cout << "Current query sequence: " << query_seqs.at(i).description
              << std::endl;
    std::cout << "The most similar sequence: "
              << target_seqs.at(max_score_idx % target_seqs_size).description
              << std::endl;
    std::cout << "The Simiarity Score: " << *max_score_ptr << std::endl
              << std::endl;
  }
}

void pair_align(FastaSequence* query_seq, FastaSequence* target_seq,
                uint32_t* score) {
  uint32_t target_seq_length = target_seq->sequence.size();
  uint32_t query_seq_length = query_seq->sequence.size();
  // similarity matrix(scoring matrix)
  std::vector<uint32_t> H;
  std::vector<int32_t> up;
  std::vector<int32_t> upleft;
  // Resize the similarity_matrix
  H.resize(2 * (target_seq_length + 1), 0);
  up.resize(target_seq_length + 1, 0);
  upleft.resize(target_seq_length + 1, 0);
  // Store the highest score in each pairwise-alignment process.
  // Default to 0.
  uint32_t max_score = 0;
  // Pairwise-Alignment between the two sequences
  __m512i H_512, up_512,
      gap_score_512 = _mm512_set1_epi32(SmithWaterman::gap_score);
  __m512i mis_score_512 = _mm512_set1_epi32(SmithWaterman::mismatch_score);
  __mmask16 strmask;
  __m128i query_seq_data;
  for (int32_t i = 1; i <= query_seq_length; i++) {
    query_seq_data = _mm_set1_epi8(query_seq->sequence.at(i - 1));
    for (int32_t j = 1; j + 15 <= target_seq_length; j += 16) {
      // From the upper element
      H_512 = _mm512_loadu_epi32(H.begin().base() + j);
      up_512 = _mm512_add_epi32(H_512, gap_score_512);
      _mm512_storeu_epi32(up.begin().base() + j, up_512);

      // From the upleft element
      strmask = _mm_cmpeq_epi8_mask(
          query_seq_data,
          _mm_loadu_epi8((target_seq->sequence.data() + j - 1)));
      __m512i upleft_score = _mm512_mask_set1_epi32(mis_score_512, strmask,
                                                    SmithWaterman::match_score);
      H_512 = _mm512_loadu_epi32(H.begin().base() + j - 1);
      _mm512_storeu_epi32(upleft.begin().base() + j,
                          _mm512_add_epi32(H_512, upleft_score));
    }
    for (int32_t j = target_seq_length - target_seq_length % 16 + 1;
         j <= target_seq_length; j++) {
      up[j] = H[j] + SmithWaterman::gap_score;
      upleft[j] =
          ((query_seq->sequence.at(i - 1) == target_seq->sequence.at(j - 1))
               ? SmithWaterman::match_score
               : SmithWaterman::mismatch_score) +
          H[j - 1];
    }
    for (int32_t j = 1; j <= target_seq_length; j++) {
      int32_t index = target_seq_length + 1 + j;

      // From the left element
      int32_t left = H[index - 1] + SmithWaterman::gap_score;

      int32_t max = std::max({up[j], left, upleft[j], 0});

      H[index] = max;

      if (max > max_score) {
        max_score = max;
      }
    }
    for (int j = 1; j <= target_seq_length; j++)
      H[j] = H[target_seq_length + 1 + j];
  }
  *score = max_score;
}

int SmithWaterman::validate(const std::string& ref_path) {
  read_ref(ref_path, refs);
  if (refs == max_scores) {
    std::cout << "Result correct!!!" << std::endl;
    report();
    return 0;
  } else {
    std::cout << "Result not match!!!" << std::endl;
    std::cout << "Reference Scores: ";
    std::copy(refs.cbegin(), refs.cend(),
              std::ostream_iterator<uint32_t>(std::cout, " "));
    std::cout << std::endl;

    std::cout << "Calculated Scores: ";
    std::copy(max_scores.cbegin(), max_scores.cend(),
              std::ostream_iterator<uint32_t>(std::cout, " "));
    std::cout << std::endl;

    return -1;
  }
}

void SmithWaterman::read_seq(const std::string& seq_path,
                             std::vector<FastaSequence>& seqs) {
  // Open file
  std::ifstream seq_file{seq_path};
  if (!seq_file.is_open()) {
    std::cerr << "Error opening file: " << seq_path << std::endl;
    throw std::runtime_error("Failed to open sequence file");
  }

  std::string line;
  FastaSequence curr_seq;

  // Read seqs
  while (std::getline(seq_file, line)) {
    if (line[0] == '>') {
      if (!curr_seq.sequence.empty()) {
        seqs.push_back(curr_seq);
        curr_seq = FastaSequence();
      }
      curr_seq.description = line.substr(1);
    } else {
      curr_seq.sequence += line;
    }
  }
  if (!curr_seq.sequence.empty()) {
    seqs.push_back(curr_seq);
  }

  // Close file
  seq_file.close();
}

void SmithWaterman::read_ref(const std::string& ref_path,
                             std::vector<uint32_t>& refs) {
  std::ifstream ref_file{ref_path};
  if (!ref_file.is_open()) {
    std::cerr << "Error opening the reference file: " << ref_path << std::endl;
    throw std::runtime_error("Failed to open reference file");
  }

  std::string line;
  while (std::getline(ref_file, line)) {
    std::istringstream iss(line);
    uint32_t score;
    while (iss >> score) {
      refs.push_back(score);
    }
  }
  ref_file.close();
}
