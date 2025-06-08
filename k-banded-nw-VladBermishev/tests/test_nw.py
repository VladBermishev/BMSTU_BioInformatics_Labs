import datetime
import src.nw as align

def test_nw_1():
    """Identical sequences, match=5, mismatch=-4, gap=-10, range=1
        results should be equal, but banded algorithm should be faster
    """
    seq1 = 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'
    seq2 = 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'
    start = datetime.datetime.now()
    result1 = align.needleman_wunsch(seq1,
                                     seq2,
                                     score=lambda x, y: 5 if x == y else -4,
                                     gap_penalty=-10)
    elapsed_time1 = datetime.datetime.now() - start
    start = datetime.datetime.now()
    result2 = align.needleman_wunsch_k(seq1,
                                     seq2,
                                     score=lambda x, y: 5 if x == y else -4,
                                     gap_penalty=-10,
                                     visible_range=1)
    elapsed_time2 = datetime.datetime.now() - start
    assert elapsed_time1 > elapsed_time2
    assert result1 == result2

def test_nw_2():
    """4 gaps, match=5, mismatch=-4, gap=-10
        in first case with visible_range = 1, gives suboptimal solution
        in second case with visible_range = 2, gives correct solution
    """
    seq1 = 'ACTGGTCAACTGGTCAACTGGTCAACTGGTCA'
    seq2 = 'TTACTGGTCAACTGGTCAACTTCAACTGGTCA'
    score_1, aligned_seq1_1, aligned_seq2_1 = align.needleman_wunsch_k(seq1,
                                                               seq2,
                                                               score=lambda x, y: 5 if x == y else -4,
                                                               gap_penalty=-10,
                                                               visible_range=1)
    score_2, aligned_seq1_2, aligned_seq2_2 = align.needleman_wunsch_k(seq1,
                                                                 seq2,
                                                                 score=lambda x, y: 5 if x == y else -4,
                                                                 gap_penalty=-10,
                                                                 visible_range=2)
    assert score_1 < score_2
    assert aligned_seq1_1 != aligned_seq1_2
    assert aligned_seq2_1 != aligned_seq2_2

def test_nw_3():
    """Full mismatch, match=5, mismatch=-4, gap=-10, visible_range=2"""
    seq1 = 'AAAAAA'
    seq2 = 'CCCCCC'
    score, aligned_seq1, aligned_seq2 = align.needleman_wunsch_k(seq1,
                                                               seq2,
                                                               score=lambda x, y: 5 if x == y else -4,
                                                               gap_penalty=-10,
                                                               visible_range=2)
    assert score == -4*len(seq1)
    assert aligned_seq1 == 'AAAAAA'
    assert aligned_seq2 == 'CCCCCC'


test_nw_1()
test_nw_2()
test_nw_3()
