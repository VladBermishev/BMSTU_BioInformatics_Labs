from typing import Callable, Tuple

DEBUG = False

def score_fun(a: str, 
              b: str,
              match_score: int = 5, 
              mismatch_score: int = -4) -> int:
    return match_score if a == b else mismatch_score

def needleman_wunsch_affine(seq1: str, 
                            seq2: str, 
                            score_fun: Callable = score_fun, 
                            gap_open: int = -10, 
                            gap_extend: int = -1) -> Tuple[str, str, int]:
    '''
    Inputs:
    seq1 - first sequence
    seq2 - second sequence
    score_fun - function that takes two characters and returns score
    gap_open - gap open penalty
    gap_extend - gap extend penalty
    Outputs:
    aln1 - first aligned sequence
    aln2 - second aligned sequence
    score - score of the alignment
    '''
    n, m = len(seq1) + 1, len(seq2) + 1
    #infinity = 2 * gap_open + (n + m - 2) * gap_extend + 1
    infinity = float('-inf')

    # 1. Initialize matrices
    match_matrix = [[0 for _ in range(m)] for _ in range(n)]
    insertion_matrix = [[0 for _ in range(m)] for _ in range(n)]
    deletion_matrix = [[0 for _ in range(m)] for _ in range(n)]
    for i in range(m):
        match_matrix[0][i] = infinity
        insertion_matrix[0][i] = infinity
        deletion_matrix[0][i] = gap_open + (i - 1) * gap_extend
    for i in range(n):
        match_matrix[i][0] = infinity
        insertion_matrix[i][0] = gap_open + (i - 1) * gap_extend
        deletion_matrix[i][0] = infinity
    match_matrix[0][0] = 0
    # 2. Fill matrices
    # We assume that consecutive gaps on different sequences are not allowed
    for i in range(1, n):
        for j in range(1, m):
            match_matrix[i][j] = max(match_matrix[i-1][j-1],
                                     insertion_matrix[i-1][j-1],
                                     deletion_matrix[i-1][j-1]) + score_fun(seq1[i-1], seq2[j-1])
            insertion_matrix[i][j] = max(insertion_matrix[i][j-1] + gap_extend,
                                         match_matrix[i][j-1] + gap_open)
            deletion_matrix[i][j] = max(deletion_matrix[i-1][j] + gap_extend,
                                        match_matrix[i-1][j] + gap_open)
    # 3. Traceback
    score = max(match_matrix[-1][-1], insertion_matrix[-1][-1], deletion_matrix[-1][-1])
    aln1, aln2 = '', ''
    i, j = len(seq1), len(seq2)
    current_matrix = "match" if score == match_matrix[-1][-1] else\
        "insertion" if score == insertion_matrix[-1][-1] else "deletion"
    while i > 0 or j > 0:
        i_step, j_step = 0, 0
        push_seq1, push_seq2 = '-', '-'
        if current_matrix == "match":
            current_matrix = "match" if match_matrix[i][j] == match_matrix[i-1][j-1] + score_fun(seq1[i-1], seq2[j-1]) else\
                "insertion" if match_matrix[i][j] == insertion_matrix[i-1][j-1] + score_fun(seq1[i-1], seq2[j-1]) else "deletion"
            push_seq1, push_seq2 = seq1[i-1], seq2[j-1]
            i_step, j_step = -1, -1
        elif current_matrix == "insertion":
            current_matrix = "insertion" if j <= 0 or insertion_matrix[i][j] == insertion_matrix[i][j-1] + gap_extend else "match"
            push_seq1, push_seq2 = ('-', seq2[j - 1]) if j > 0 else (seq1[i - 1], '-')
            i_step, j_step = (0, -1) if j > 0 else (-1, 0)
        else:
            current_matrix = "deletion" if i <= 0 or deletion_matrix[i][j] == deletion_matrix[i-1][j] + gap_extend else "match"
            push_seq1, push_seq2 = (seq1[i - 1], '-') if i > 0 else ('-', seq2[j - 1])
            i_step, j_step = (-1, 0) if i > 0 else (0, -1)
        i, j = i + i_step, j + j_step
        aln1 += push_seq1
        aln2 += push_seq2
    
    return aln1[::-1], aln2[::-1], score

def print_array(matrix: list):
    for row in matrix:
        for element in row:
            print(f"{element:6}", end="")
        print()

def main():
    aln1, aln2, score = needleman_wunsch_affine("ACGT", "TAGT", gap_open=-10, gap_extend=-1) 
    print(f'str 1: {aln1}')
    print(f'str 2: {aln2}')
    print(f'score: {score}')
    


if __name__ == "__main__":
    main()