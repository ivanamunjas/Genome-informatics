# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 13:42:51 2020

@author: Ivana Pesic
"""

import burrows_wheeler

bw = burrows_wheeler.SuffixArrayBurrowsWheeler()

# calculate the first occurence of letter in sorted string
def calc_first_occ(s):
    # s - is the bwt transformed string
    A = {}  # letter count
    for i, c in enumerate(s):
        if A.get(c):
            A[c] += 1
        else:
            A[c] = 1

    # sort the letters
    letters = sorted(A.keys())

    # first index of letter
    occ = {}

    idx = 0
    for c in letters:
        occ[c] = idx
        idx += A[c]
    del idx, A

    return occ

# count the number of letters for each tally matrix step and return the list of the counts
def calc_checkpoints(s, step):
    A = {}  # letter count
    C = []  # checkpoints
    for i, c in enumerate(s):
        if i % step == 0:
            C.append(A.copy())
        if A.get(c):
            A[c] += 1
        else:
            A[c] = 1
    return C

# count the number of letters for each suffix array step and return list of the counts
def calc_sa_checkpoints(sa, steps):
    sa_temp = []
    for i in range(len(sa)):
        if (sa[i] % steps) != 0:
            sa_temp.append(None)
        else:
            sa_temp.append(sa[i])
    return sa_temp


def count_letter_with_checkpoints(C, step, s, idx, letter):
    """ 
    Arguments:
    C      - the list of checkpoints
    step   - the step of the checkpoints
    s      - transformed string
    idx    - count upto this position
    letter - count for this letter
    """

    # find the nearest checkpoint for idx
    check = int((idx + (step / 2)) / step)
    if check >= len(C):
        check = len(C) - 1
    pos = check * step

    # count of the letter s[idx] upto pos (not included)
    count = C[check].get(letter)
    if count == None:
        count = 0

    # range between pos and idx
    if pos < idx:
        r = range(pos, idx)
    else:
        r = range(idx, pos)

    # count of letters between pos, idx
    k = 0
    for i in r:
        if letter == s[i]:
            k += 1

    # calculate the letter count upto idx (not included)
    if pos < idx:
        count += k
    else:
        count -= k

    return count


class FMSimpleIndex(object):
    def __init__(self, data):
        self.data = bw.transform(data)
        self.offset = {}
        self._build(data)

    # build index
    def _build(self, data):
        self.occ = calc_first_occ(self.data)

    # get occurance of specific letter
    def _occ(self, qc):
        c = self.occ.get(qc)
        if c == None:
            return 0
        return c

    # count occurances of qc up to idx
    def _count(self, idx, qc):
        if not qc in self.occ.keys(): return 0
        c = 0
        for i in range(idx):
            if self.data[i] == qc:
                c += 1
        return c

    # get LF mapping for letter qc at position idx
    def _lf(self, idx, qc):
        o = self._occ(qc)
        c = self._count(idx, qc)
        return o + c

    def suffix(self, i):
        count = 0
        while self.sa[i] == None:
            i = self._lf(i, self.data[i])
            count = count + 1
        return self.sa[i] + count

    # find the first and the last suffix position
    def bounds(self, q):
        top = 0
        bot = len(self.data)
        for i, qc in enumerate(q[::-1]):
            top = self._lf(top, qc)
            bot = self._lf(bot, qc)
            if top == bot: return (-1, -1)
        return (top, bot)

    def search(self, q):
        # find the suffixes for the query
        top, bot = self.bounds(q)
        matches = []
        # find the location of the suffixes by walking the reverse text from that position
        # with lf mapping
        for i in range(top, bot):
            pos = self.suffix(i)
            matches.append(pos)
        return sorted(matches)

    # count occurances
    def count(self, q):
        top, bot = self.bounds(q)
        return bot - top

# FM index with checkpointing
class FMCheckpointing(FMSimpleIndex):

    def __init__(self, data, sa_step = 10, tally_step = 10):
        bwt_ref, idx, self.sa = bw.bwt(data)
        encoded = "".join(data[x] for x in bwt_ref)
        self.data = encoded
        self.offset = {}
        self.step = tally_step
        self.step_sa = sa_step
        self._build()
    
    # building index
    def _build(self):
        self.occ = calc_first_occ(self.data)
        self.C = calc_checkpoints(self.data, self.step)
        self.sa = calc_sa_checkpoints(self.sa, self.step_sa)
 
    # count occurances of specific letter
    def _count(self, idx, qc):
        count = count_letter_with_checkpoints(self.C, self.step, self.data, idx, qc)
        return count