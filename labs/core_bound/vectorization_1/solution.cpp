#include "solution.hpp"
#include <algorithm>
#include <cassert>
#include <type_traits>

using score_t = std::int16_t;
using symbol_t = std::uint8_t;
using sequence_column_t = std::array<symbol_t, sequence_count_v>;
using sequence_group_t = std::array<sequence_column_t, sequence_size_v>;

sequence_group_t make_sequence_group(const std::vector<sequence_t>& sequences) {
  assert(size(sequences) == sequence_count_v);
  sequence_group_t result;
  for (std::size_t i=0; i<size(sequences); ++i) {
    for (std::size_t j=0; j<sequence_size_v; ++j) {
      result[j][i] = sequences[i][j];
    }
  }
  return result;
}

result_t compute_alignment(const sequence_group_t& sequences1,
                           const sequence_group_t& sequences2)
{
  constexpr score_t gap_open{-11};
  constexpr score_t gap_extension{-1};
  constexpr score_t match{6};
  constexpr score_t mismatch{-4};

  using score_column_t = std::array<score_t, sequence_count_v>;
  using score_mat_t = std::array<score_column_t, sequence_size_v + 1>;

  score_mat_t scores;
  score_mat_t horizontal_gaps;

  for (std::size_t i=0; i<sequence_count_v; ++i) {
    horizontal_gaps[0][i] = gap_open;
  }

  for (std::size_t i=1; i<size(scores); ++i) {
    for (std::size_t j=0; j<sequence_count_v; ++j) {
      scores[i][j] = gap_open + (i - 1) * gap_extension;
      horizontal_gaps[i][j] = 2 * gap_open + (i - 1) * gap_extension;
    }
  }

  score_column_t last_vertical_gaps;
  score_column_t last_diagonal_scores;

  for (std::size_t col=1; col<=sequence_size_v; ++col) {
    for (std::size_t seq=0; seq<sequence_count_v; ++seq) {
      last_diagonal_scores[seq] = scores[0][seq];
      scores[0][seq] = horizontal_gaps[0][seq];
      last_vertical_gaps[seq] = horizontal_gaps[0][seq] + gap_open;
      horizontal_gaps[0][seq] += gap_extension;
    }

    for (std::size_t row=1; row<=sequence_size_v; ++row) {
      for (std::size_t seq=0; seq<sequence_count_v; ++seq) {
        // Compute next score from diagonal direction with match/mismatch.
        score_t best_cell_score =
            last_diagonal_scores[seq] +
            (sequences1[row - 1][seq] == sequences2[col - 1][seq] ? match : mismatch);
        // Determine best score from diagonal, vertical, or horizontal
        // direction.
        best_cell_score = std::max(best_cell_score, last_vertical_gaps[seq]);
        best_cell_score = std::max(best_cell_score, horizontal_gaps[row][seq]);
        // Cache next diagonal value and store optimum in score_column.
        last_diagonal_scores[seq] = scores[row][seq];
        scores[row][seq] = best_cell_score;
        // Compute the next values for vertical and horizontal gap.
        best_cell_score += gap_open;
        last_vertical_gaps[seq] += gap_extension;
        horizontal_gaps[row][seq] += gap_extension;
        // Store optimum between gap open and gap extension.
        last_vertical_gaps[seq] = std::max(last_vertical_gaps[seq], best_cell_score);
        horizontal_gaps[row][seq] =
            std::max(horizontal_gaps[row][seq], best_cell_score);
      }
    }
  }

  return scores.back();
}

// The alignment algorithm which computes the alignment of the given sequence
// pairs.
result_t compute_alignment(std::vector<sequence_t> const &sequences1,
                           std::vector<sequence_t> const &sequences2) {
#ifdef SOLUTION
  return compute_alignment(
    make_sequence_group(sequences1),
    make_sequence_group(sequences2));
#endif

  result_t result{};

  for (size_t sequence_idx = 0; sequence_idx < sequences1.size();
       ++sequence_idx) {
    using score_t = int16_t;
    using column_t = std::array<score_t, sequence_size_v + 1>;

    sequence_t const &sequence1 = sequences1[sequence_idx];
    sequence_t const &sequence2 = sequences2[sequence_idx];

    /*
     * Initialise score values.
     */
    score_t gap_open{-11};
    score_t gap_extension{-1};
    score_t match{6};
    score_t mismatch{-4};

    /*
     * Setup the matrix.
     * Note we can compute the entire matrix with just one column in memory,
     * since we are only interested in the last value of the last column in the
     * score matrix.
     */
    column_t score_column{};
    column_t horizontal_gap_column{};
    score_t last_vertical_gap{};

    /*
     * Initialise the first column of the matrix.
     */
    horizontal_gap_column[0] = gap_open;
    last_vertical_gap = gap_open;

    for (size_t i = 1; i < score_column.size(); ++i) {
      score_column[i] = last_vertical_gap;
      horizontal_gap_column[i] = last_vertical_gap + gap_open;
      last_vertical_gap += gap_extension;
    }

    /*
     * Compute the main recursion to fill the matrix.
     */
    for (unsigned col = 1; col <= sequence2.size(); ++col) {
      score_t last_diagonal_score =
          score_column[0]; // Cache last diagonal score to compute this cell.
      score_column[0] = horizontal_gap_column[0];
      last_vertical_gap = horizontal_gap_column[0] + gap_open;
      horizontal_gap_column[0] += gap_extension;

      for (unsigned row = 1; row <= sequence1.size(); ++row) {
        // Compute next score from diagonal direction with match/mismatch.
        score_t best_cell_score =
            last_diagonal_score +
            (sequence1[row - 1] == sequence2[col - 1] ? match : mismatch);
        // Determine best score from diagonal, vertical, or horizontal
        // direction.
        best_cell_score = std::max(best_cell_score, last_vertical_gap);
        best_cell_score = std::max(best_cell_score, horizontal_gap_column[row]);
        // Cache next diagonal value and store optimum in score_column.
        last_diagonal_score = score_column[row];
        score_column[row] = best_cell_score;
        // Compute the next values for vertical and horizontal gap.
        best_cell_score += gap_open;
        last_vertical_gap += gap_extension;
        horizontal_gap_column[row] += gap_extension;
        // Store optimum between gap open and gap extension.
        last_vertical_gap = std::max(last_vertical_gap, best_cell_score);
        horizontal_gap_column[row] =
            std::max(horizontal_gap_column[row], best_cell_score);
      }
    }

    // Report the best score.
    result[sequence_idx] = score_column.back();
  }

  return result;
}
