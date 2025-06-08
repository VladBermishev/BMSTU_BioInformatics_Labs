from typing import Callable, Tuple
import argparse
import sys

PRINT_MAX_LINE_LENGTH = 80
DEBUG = False
GLOBAL_MINIMUM = -10**10


def score_fun(a: str,  b: str, match_score: int = 5, mismatch_score: int = -4) -> int:
    return match_score if a == b else mismatch_score


def needleman_wunsch(seq1: str,
                     seq2: str,
                     score: Callable[[str, str], int] = score_fun,
                     gap_penalty: int = -10,
                     visible_range: int = None) -> Tuple[int, str, str]:

    """Given two sequences, aligns them using the Needleman-Wunsch algorithm.

    This function takes two sequences and optionally a scoring function and a
    gap penalty value as arguments.
    The function returns a tuple containing the optimal alignment score and the
    aligned sequences, e.g. (10, 'ACCGT', 'AC-GT').

    Args:
        seq1: The first sequence, e.g. 'CCGT'
        seq2: The second sequence, e.g. 'ACGT'
        score: The scoring function, e.g. score_fun('A', 'A') returns 5
        gap_penalty: The gap penalty value, e.g. -10
        visible_range: The width of limited diagonal (<= len(seq1)), e.g. 5

    Returns:
        score: The optimal alignment score, e.g. 10
        aligned_seq1: The first aligned sequence, e.g. 'ACCGT'
        aligned_seq2: The second aligned sequence, e.g. 'AC-GT'
    """

    """ Initialize the score matrix.
        If visible_range is defined we won't evaluate whole matrix, so it's important to set everything to GLOBAL_MINIMUM 
        for restored aligned sequences to be correct
    """
    default_value = GLOBAL_MINIMUM if visible_range else 0
    score_matrix = [[default_value for _ in range(len(seq2) + 1)] for _ in range(len(seq1) + 1)]
    # Setting seeds: D(i,0) = i*gap, D(0,j) = j*gap
    for i in range(len(seq1) + 1): score_matrix[i][0] = i * gap_penalty
    for i in range(len(seq2) + 1): score_matrix[0][i] = i * gap_penalty
    """ Evaluating score_matrix,
        if visible range is defined and |i - j| > visible_range than score_matrix[i][j] shouldn't be evaluated
    """
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            if visible_range and abs(i - j) > visible_range: continue
            score_max = max(score_matrix[i - 1][j - 1] + score(seq1[i - 1], seq2[j - 1]),
                            score_matrix[i - 1][j] + gap_penalty,
                            score_matrix[i][j - 1] + gap_penalty)
            score_matrix[i][j] = score_max
    # Restoring aligned sequences: finding max among diag, up, left, then go there until i == 0 and j == 0
    i, j = len(seq1), len(seq2)
    aligned_seq1, aligned_seq2 = '', ''
    while i != 0 or j != 0:
        diag, up, left = (
        score_matrix[i - 1][j - 1] + score(seq1[i - 1], seq2[j - 1]) if i >= 1 and j >= 1 else GLOBAL_MINIMUM,
        score_matrix[i - 1][j] + gap_penalty if i >= 1 else GLOBAL_MINIMUM,
        score_matrix[i][j - 1] + gap_penalty if j >= 1 else GLOBAL_MINIMUM)
        i_step, j_step = 0, 0
        push_seq1, push_seq2 = '-', '-'
        if diag >= up and diag >= left:
            push_seq1, push_seq2 = seq1[i - 1], seq2[j - 1]
            i_step, j_step = -1, -1
        elif up >= diag and up >= left:
            push_seq1 = seq1[i - 1]
            i_step = -1
        elif left >= up and left >= diag:
            push_seq2 = seq2[j - 1]
            j_step = -1
        i, j = i + i_step, j + j_step
        aligned_seq1 = push_seq1 + aligned_seq1
        aligned_seq2 = push_seq2 + aligned_seq2
    return score_matrix[-1][-1], aligned_seq1, aligned_seq2

def needleman_wunsch_k(seq1: str,
                     seq2: str,
                     score: Callable[[str, str], int] = score_fun,
                     gap_penalty: int = -10,
                     visible_range: int = 1) -> Tuple[int, str, str]:

    """Given two sequences, aligns them using the improved k-banded Needleman-Wunsch algorithm.

    This function takes two sequences with same length, optionally a scoring function, a gap penalty value,
     visible_range representing k as arguments.
    The function returns a tuple containing the optimal alignment score and the
    aligned sequences, e.g. (10, 'ACCGT', 'AC-GT').

    Args:
        Sequences must be the same length
        seq1: The first sequence, e.g. 'CCGT'
        seq2: The second sequence, e.g. 'ACGT'
        score: The scoring function, e.g. score_fun('A', 'A') returns 5
        gap_penalty: The gap penalty value, e.g. -10
        visible_range: The width of limited diagonal (<= len(seq1)), e.g. 5
    Returns:
        score: The optimal alignment score, e.g. 10
        aligned_seq1: The first aligned sequence, e.g. 'ACCGT'
        aligned_seq2: The second aligned sequence, e.g. 'AC-GT'
    """
    assert len(seq1) == len(seq2)
    return needleman_wunsch(seq1, seq2, score, gap_penalty, visible_range)

