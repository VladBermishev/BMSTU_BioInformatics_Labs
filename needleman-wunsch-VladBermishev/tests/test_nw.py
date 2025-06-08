import src.nw as align


def test_nw_1():
    score, aln1, aln2 = align.needleman_wunsch("ACGT", "ACGT")
    assert len(aln1) == len(aln2)
    assert aln1 == "ACGT"
    assert score == 20


def test_nw_2():
    score, aln1, aln2 = align.needleman_wunsch("ACG", "ACGT")
    assert len(aln1) == len(aln2)
    assert aln1 == "ACG-"
    assert aln2 == "ACGT"
    assert score == 5


def test_nw_3():
    score, aln1, aln2 = align.needleman_wunsch("ACGT", "ACG")
    assert len(aln1) == len(aln2)
    assert aln1 == "ACGT"
    assert aln2 == "ACG-"
    assert score == 5


def test_nw_4():
    score, aln1, aln2 = align.needleman_wunsch("ACAGT", "ACGT")
    assert len(aln1) == len(aln2)
    assert aln1 == "ACAGT"
    assert aln2 == "AC-GT"
    assert score == 10


def test_nw_5():
    score, aln1, aln2 = align.needleman_wunsch("ACGT", "ACAGT")
    assert len(aln1) == len(aln2)
    assert aln1 == "AC-GT"
    assert aln2 == "ACAGT"
    assert score == 10


def test_nw_6():
    score, aln1, aln2 = align.needleman_wunsch("CAGT", "ACAGT")
    assert len(aln1) == len(aln2)
    assert aln1 == "-CAGT"
    assert aln2 == "ACAGT"
    assert score == 10


def test_nw_7():
    score, aln1, aln2 = align.needleman_wunsch("ACAGT", "CAGT")
    assert len(aln1) == len(aln2)
    assert aln1 == "ACAGT"
    assert aln2 == "-CAGT"
    assert score == 10


def test_nw_8():
    score, aln1, aln2 = align.needleman_wunsch("ACGT", "A")
    assert len(aln1) == len(aln2)
    assert aln1 == "ACGT"
    assert aln2 == "A---"
    assert score == -25


def test_nw_9():
    score, aln1, aln2 = align.needleman_wunsch("ACGT", "")
    assert len(aln1) == len(aln2)
    assert aln1 == "ACGT"
    assert aln2 == "----"
    assert score == -40


def test_nw_10():
    score, aln1, aln2 = align.needleman_wunsch("A", "ACGT")
    assert len(aln1) == len(aln2)
    assert aln1 == "A---"
    assert aln2 == "ACGT"
    assert score == -25


def test_nw_11():
    score, aln1, aln2 = align.needleman_wunsch("", "ACGT")
    assert len(aln1) == len(aln2)
    assert aln1 == "----"
    assert aln2 == "ACGT"
    assert score == -40


def test_nw_12():
    score, aln1, aln2 = align.needleman_wunsch("", "")
    assert aln1 == ""
    assert aln2 == ""
    assert score == 0


def test_nw_13():
    score, aln1, aln2 = align.needleman_wunsch("TACGT", "ATGT")
    assert len(aln1) == len(aln2)
    assert aln1 == "TACGT"
    assert aln2 == "-ATGT"
    assert score == 1


def test_nw_14():
    score, aln1, aln2 = align.needleman_wunsch("TACGT", "ACTGT")
    assert len(aln1) == len(aln2)
    assert aln1 == "TAC-GT"
    assert aln2 == "-ACTGT"
    assert score == 0


def test_nw_15():
    score, aln1, aln2 = align.needleman_wunsch("ACGT", "TAGTA", gap_penalty=-5)
    assert len(aln1) == len(aln2)
    assert aln1 == "-ACGT-"
    assert aln2 == "TA-GTA"
    assert score == 0


def test_nw_16():
    score, aln1, aln2 = align.needleman_wunsch("TAGTA", "ACGT", gap_penalty=-5)
    assert len(aln1) == len(aln2)
    assert aln1 == "TA-GTA"
    assert aln2 == "-ACGT-"
    assert score == 0


def test_nw_17():
    score, aln1, aln2 = align.needleman_wunsch("ACGT", "TAGT", gap_penalty=0)
    assert len(aln1) == len(aln2)
    assert aln1 == "-ACGT"
    assert aln2 == "TA-GT"
    assert score == 15


def test_nw_18():
    score, aln1, aln2 = align.needleman_wunsch("TAGT", "ACGT", gap_penalty=10)
    assert len(aln1) == len(aln2)
    assert len(aln1) == 8
    assert score == 80


def test_nw_19():
    score, aln1, aln2 = align.needleman_wunsch("GGAGCCAAGGTGAAGTTGTAGCAGTGTGTCC",
                                         "GACTTGTGGAACCTCTGTCCTCCGAGCTCTC", gap_penalty=-5)
    assert len(aln1) == len(aln2)
    assert len(aln1) == 36
    assert score == 8


def test_nw_20():
    score, aln1, aln2 = align.needleman_wunsch("AAAAAAATTTTTTT", "TTTTTTTAAAAAAA", gap_penalty=-5)
    assert len(aln1) == len(aln2)
    assert len(aln1) == 21
    assert score == -35

def test_nw_21():
    """Full match, except 2 gaps and 1 mismatch, match=5, mismatch=-4, gap=-10"""
    seq1 = 'ACTGGTCAACTGGTCAACTGGTCAACTGGTCA'
    seq2 = 'TACTGGTCAACTGGTCAACTGTCAACTGGTCA'
    score, aligned_seq1, aligned_seq2 = align.needleman_wunsch(seq1,
                                                               seq2,
                                                               score=lambda x, y: 5 if x == y else -4,
                                                               gap_penalty=-10)
    assert score == 135
    assert aligned_seq1 == '-ACTGGTCAACTGGTCAACTGGTCAACTGGTCA'
    assert aligned_seq2 == 'TACTGGTCAACTGGTCAACT-GTCAACTGGTCA'

def test_nw_22():
    """Full mismatch, match=5, mismatch=-4, gap=-10"""
    seq1 = 'AAAAAA'
    seq2 = 'CCCCCC'
    score, aligned_seq1, aligned_seq2 = align.needleman_wunsch(seq1,
                                                               seq2,
                                                               score=lambda x, y: 5 if x == y else -4,
                                                               gap_penalty=-10)
    assert score == -4*len(seq1)
    assert aligned_seq1 == 'AAAAAA'
    assert aligned_seq2 == 'CCCCCC'