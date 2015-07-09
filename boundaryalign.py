from __future__ import division, print_function, unicode_literals
from operator import itemgetter

import numpy as np
import pandas
from util import by_chrom, chrom_sorted


GAP = 666 # gaps are evil!
DIAG = 0
UP = 1
LEFT = 2


def backtrack(B, S):
    i, j = S.shape[0] - 1, S.shape[1] - 1
    path = [(i,j, -1)]
    while not (i == 0 and j == 0):
        if B[i,j] == DIAG:   # match
            i -= 1
            j -= 1
            path.append((i,j, DIAG))
        elif B[i,j] == UP:   # gap in y
            i -= 1
            path.append((i,j, UP))
        elif B[i,j] == LEFT: # gap in x
            j -= 1
            path.append((i,j, LEFT))
    return path[::-1]


def boundary_align(chrom, seq1, seq2, gap_cost=20000):
    m = len(seq1)
    n = len(seq2)
    
    DP = np.zeros((m+1, n+1))
    for i in xrange(1, m+1):
        DP[i, 0] = DP[i-1, 0] + gap_cost
    for j in xrange(1, n+1):
        DP[0, j] = DP[0, j-1] + gap_cost
    
    B = np.zeros((m+1, n+1), dtype=int)
    B[0, :] = LEFT
    B[:, 0] = UP
    
    for i in xrange(1, m+1):
        for j in xrange(1, n+1):
            choices = (
                (DIAG, DP[i-1, j-1] + abs(seq1[i-1] - seq2[j-1])),
                (UP,   DP[i-1, j] + gap_cost),
                (LEFT, DP[i, j-1] + gap_cost)
            )
            B[i,j], DP[i,j] = min(choices, key=itemgetter(1))
            
    cost = DP[m, n]
    path = backtrack(B, DP)
    
    return np.array(path), cost
