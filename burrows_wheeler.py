# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 13:41:05 2020

@author: Ivana Pesic
"""

from operator import itemgetter
import collections

class BurrowsWheeler():
    
    # end marker
    EOS = "\0"

    # bwt and suffix array
    def suffix_array(self, arr):
        arr_size = len(arr)
        arr_int = {v: k for k, v in enumerate(sorted(set(arr)))}
        arr = [arr_int[x] for x in arr]
        arr.append(-1)
        suf = [[i, arr[i], arr[i + 1]] for i in range(arr_size)]
        # suf = radixSort(suf)
        suf.sort(key=itemgetter(1, 2))
        idx = [0] * arr_size
        k = 2
        while k < arr_size:
            r = 0
            prev_r = suf[0][1]
            for i in range(arr_size):
                if suf[i][1] != prev_r or suf[i - 1][2] != suf[i][2]:
                    r += 1
                prev_r = suf[i][1]
                suf[i][1] = r
                idx[suf[i][0]] = i
            for i in range(arr_size):
                next_idx = suf[i][0] + k
                suf[i][2] = suf[idx[next_idx]][1] if next_idx < arr_size else -1
            suf.sort(key=itemgetter(1, 2))
            k <<= 1
        return [x[0] for x in suf]

    def bwt(self,input_string):
        data_ref = self.suffix_array(input_string)
        return (x - 1 for x in data_ref), data_ref.index(0), data_ref

class SuffixArrayBurrowsWheeler(BurrowsWheeler):

    def transform(self, s):

        assert self.EOS not in s, "Input string cannot contain null character (%s)" % self.EOS

        # add end of text marker
        s += self.EOS

        # table of suffixes
        rotations = [s[i:] for i in range(len(s))]

        sa = {}
        for i, q in enumerate(rotations):
            sa[q] = i
        od = collections.OrderedDict(sorted(sa.items()))

        # sort the suffixes
        rotations.sort()

        # get the length of ordered suffixes
        k = len(rotations)

        r = [0] * k
        for i in range(k):
            l = len(rotations[i])
            if l == k:
                r[i] = self.EOS
            else:
                r[i] = s[-l - 1]
        r = ''.join(r)

        return r, list(od.values())