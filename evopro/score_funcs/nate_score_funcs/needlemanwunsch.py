#!/usr/bin/env python
"""
The Needleman-Wunsch Algorithm
==============================

This is a dynamic programming algorithm for finding the optimal alignment of
two strings.

Example
-------

    >>> x = "GATTACA"
    >>> y = "GCATGCU"
    >>> print(nw(x, y))
    G-ATTACA
    GCA-TGCU


LICENSE

This is free and unencumbered software released into the public domain.
Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <http://unlicense.org/>
"""

import numpy as np

def nw(x, y, match = 1, mismatch = 1, gap = 1):
    nx = len(x)
    ny = len(y)
    # Optimal score at each possible pair of characters.
    F = np.zeros((nx + 1, ny + 1))
    F[:,0] = np.linspace(0, -nx * gap, nx + 1)
    F[0,:] = np.linspace(0, -ny * gap, ny + 1)
    # Pointers to trace through an optimal aligment.
    P = np.zeros((nx + 1, ny + 1))
    P[:,0] = 3
    P[0,:] = 4
    # Temporary scores.
    t = np.zeros(3)
    for i in range(nx):
        for j in range(ny):
            if x[i] == y[j]:
                t[0] = F[i,j] + match
            else:
                t[0] = F[i,j] - mismatch
            t[1] = F[i,j+1] - gap
            t[2] = F[i+1,j] - gap
            tmax = np.max(t)
            F[i+1,j+1] = tmax
            if t[0] == tmax:
                P[i+1,j+1] += 2
            if t[1] == tmax:
                P[i+1,j+1] += 3
            if t[2] == tmax:
                P[i+1,j+1] += 4
    # Trace through an optimal alignment.
    i = nx
    j = ny
    rx = []
    ry = []
    while i > 0 or j > 0:
        if P[i,j] in [2, 5, 6, 9]:
            rx.append(x[i-1])
            ry.append(y[j-1])
            i -= 1
            j -= 1
        elif P[i,j] in [3, 5, 7, 9]:
            rx.append(x[i-1])
            ry.append('-')
            i -= 1
        elif P[i,j] in [4, 6, 7, 9]:
            rx.append('-')
            ry.append(y[j-1])
            j -= 1
    # Reverse the strings.
    rx = ''.join(rx)[::-1]
    ry = ''.join(ry)[::-1]
    return '\n'.join([rx, ry])

#x = "GATTACA"
#y = "GCATGCU"
#print(nw(x, y))
# G-ATTACA
# GCA-TGCU

#np.random.seed(42)
#x = np.random.choice(['A', 'T', 'G', 'C'], 50)
#y = np.random.choice(['A', 'T', 'G', 'C'], 50)

#print(nw(x, y, gap = 0))
# ----G-C--AGGCAAGTGGGGCACCCGTATCCT-T-T-C-C-AACTTACAAGGGT-C-CC-----CGT-T
# GTGCGCCAGAGG-AAGT----CA--C-T-T--TATATCCGCG--C--AC---GGTACTCCTTTTTC-TA-

#print(nw(x, y, gap = 1))
# GCAG-GCAAGTGG--GGCAC-CCGTATCCTTTC-CAAC-TTACAAGGGTCC-CCGT-T-
# G-TGCGCCAGAGGAAGTCACTTTATATCC--GCGC-ACGGTAC-----TCCTTTTTCTA

#print(nw(x, y, gap = 2))
# GCAGGCAAGTGG--GGCAC-CCGTATCCTTTCCAACTTACAAGGGTCCCCGTT
# GTGCGCCAGAGGAAGTCACTTTATATCC-GCGCACGGTAC-TCCTTTTTC-TA

