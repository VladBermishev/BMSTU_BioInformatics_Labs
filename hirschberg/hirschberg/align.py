from typing import Callable, Tuple

DEBUG = False

def score_fun(a: str, 
              b: str,
              match_score: int = 5, 
              mismatch_score: int = -4) -> int:
    return match_score if a == b else mismatch_score


def needleman_wunsch(seq1: str, seq2: str, score_fun: Callable = score_fun, gap_score: int = -5):
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

        Returns:
            score: The optimal alignment score, e.g. 10
            aligned_seq1: The first aligned sequence, e.g. 'ACCGT'
            aligned_seq2: The second aligned sequence, e.g. 'AC-GT'
        """

    m, n = len(seq1) + 1, len(seq2) + 1
    matrix = [[0] * n for _ in range(m)]

    for i in range(m):
        matrix[i][0] = i * gap_score
    for j in range(n):
        matrix[0][j] = j * gap_score

    for i in range(1, m):
        for j in range(1, n):
            matrix[i][j] = max(matrix[i - 1][j - 1] + score_fun(seq1[i - 1], seq2[j - 1]),
                               matrix[i - 1][j] + gap_score,
                               matrix[i][j - 1] + gap_score)
    if DEBUG:
        print_array(matrix)

    score = matrix[-1][-1]
    i, j = m - 1, n - 1
    aln1 = ""
    aln2 = ""
    while i > 0 or j > 0:
        a, b = '-', '-'
        # (A, B)
        if i > 0 and j > 0 and matrix[i][j] == matrix[i - 1][j - 1] + score_fun(seq1[i - 1], seq2[j - 1]):
            a = seq1[i - 1]
            b = seq2[j - 1]
            i -= 1
            j -= 1

        # (A, -)
        elif i > 0 and matrix[i][j] == matrix[i - 1][j] + gap_score:
            a = seq1[i - 1]
            i -= 1

        # (-, A)
        elif j > 0 and matrix[i][j] == matrix[i][j - 1] + gap_score:
            b = seq2[j - 1]
            j -= 1

        aln1 += a
        aln2 += b
    return aln1[::-1], aln2[::-1], score

def nw_score_evaluate(seq1: str, seq2: str, score: Callable = score_fun, gap_score: int = -5):
    N, M = len(seq1) + 1, len(seq2) + 1
    score_rows = [[0 for j in range(len(seq2) + 1)] for i in range(2)]
    for j in range(M):
        score_rows[0][j] = j*gap_score
    for i in range(2):
        score_rows[i][0] = i*gap_score

    for i in range(1, N):
        score_rows[(i-1) % 2][0] = (i-1)*gap_score
        score_rows[i % 2][0] = i*gap_score
        for j in range(1, M):
            score_rows[i % 2][j] = max(score_rows[(i-1) % 2][j-1] + score(seq1[i-1], seq2[j-1]),
                                   score_rows[(i-1) % 2][j] + gap_score,
                                   score_rows[i % 2][j-1] + gap_score)
    return score_rows[(N+1) % 2]


def hirschberg(seq1: str, 
               seq2: str, 
               score: Callable = score_fun,
               gap_score: int = -5) -> Tuple[str, str, int]:
    '''
    Inputs:
    seq1 - first sequence
    seq2 - second sequence
    score_fun - function that returns score for two symbols
    gap_score - score for gap in final alignment

    Outputs:
    aln1 - first sequence in alignment
    aln2 - second sequence in alignment
    score - score of alignment
    '''
    # Dealing with end of recursion when |sequences| < 2
    if len(seq1) == 0:
        return '-' * len(seq2), seq2, gap_score*len(seq2)
    elif len(seq2) == 0:
        return seq1, '-' * len(seq1), gap_score*len(seq1)
    elif len(seq1) == 1 or len(seq2) == 1:
        return needleman_wunsch(seq1, seq2, score, gap_score)
    else:
        i = len(seq1) // 2
        # Evaluating needleman_wunsch scores for upper and lower halfs
        s_up = nw_score_evaluate(seq1[:i], seq2, score, gap_score)
        s_down = nw_score_evaluate(seq1[len(seq1):i-1:-1], seq2[::-1], score, gap_score)
        # Finding pivot index for columns & continue evaluating on smaller matrices
        s = [s_up[i] + s_down[len(seq2) - i] for i in range(len(seq2) + 1)]
        j = s.index(max(s))
        res_up = hirschberg(seq1[:i], seq2[:j], score, gap_score)
        res_down = hirschberg(seq1[i:len(seq1)], seq2[j:len(seq2)], score, gap_score)
        return res_up[0] + res_down[0], res_up[1] + res_down[1], res_up[2] + res_down[2]

def print_array(matrix: list):
    for row in matrix:
        for element in row:
            print(f"{element:6}", end="")
        print()



if __name__ == "__main__":    
    #aln1, aln2, score = hirschberg("ATCT", "ACT", gap_score=-5)
    aln1, aln2, score = needleman_wunsch("ATCT", "ACT", gap_score=-5)
    
    assert len(aln1) == len(aln2)
    print(aln1)
    print(aln2)
    print(score)
