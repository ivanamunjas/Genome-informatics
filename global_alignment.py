# -*- coding: utf-8 -*-
"""
Created on Sun Nov 22 11:11:07 2020

@author: Ivana Pesic
"""

import numpy

class GlobalAlignment():
    
    def __init__(self, match, mismatch, gap):
        self.match = match
        self.mismatch = mismatch  
        self.gap = gap

    def scoring_matrix(self, a, b):
        if a == b: return self.match
        if a == '_' or b == '_' : return self.gap
        return self.mismatch
    
    def global_alignment(self, x, y, s):
        D = numpy.zeros((len(x) + 1, len(y) + 1), dtype=int)
        
        for i in range(1, len(x) + 1):
            D[i,0] = D[i-1,0] + s(x[i-1], '_')  
        for j in range(1, len(y)+1):
            D[0,j] = D[0,j-1] + s('_', y[j-1])
        
        for i in range(1, len(x) + 1):
            for j in range(1, len(y) + 1):
                D[i,j] = max(D[i-1,j]   + s(x[i-1], '_'),
                             D[i,j-1]   + s('_', y[j-1]), 
                             D[i-1,j-1] + s(x[i-1], y[j-1]))
                
        # function returns table and global alignment score
        #alignment score is in cell (n,m) of the matrix
        return D, D[len(x),len(y)]
    
    
    def traceback(self, x, y, V, s):
        
        # initializing starting position cell(n,m)
        i=len(x)
        j=len(y)
        
        # initializing strings we use to represent alignments in x, y, edit transcript and global alignment
        ax, ay, am, tr = '', '', '', ''
        
        # exit condition is when we reach cell (0,0)
        while i > 0 or j > 0:
            
            # calculating diagonal, horizontal and vertical scores for current cell
            d, v, h = -100, -100, -100
            
            if i > 0 and j > 0:
                delta = 1 if x[i-1] == y[j-1] else 0
                d = V[i-1,j-1] + s(x[i-1], y[j-1])  # diagonal movement   
            if i > 0: v = V[i-1,j] + s(x[i-1], '_')  # vertical movement
            if j > 0: h = V[i,j-1] + s('_', y[j-1])  # horizontal movement
                
            # backtracing to next (previous) cell
            if d >= v and d >= h:
                ax += x[i-1]
                ay += y[j-1]
                if delta == 1:
                    tr += 'M'
                    am += '|'
                else:
                    tr += 'R'
                    am += ' '
                i -= 1
                j -= 1
            elif v >= h:
                ax += x[i-1]
                ay += '_'
                tr += 'D'
                am += ' '
                i -= 1
            else:
                ay += y[j-1]
                ax += '_'
                tr += 'I'
                am += ' '
                j -= 1
                
        alignment='\n'.join([ax[::-1], am[::-1], ay[::-1]])
        return alignment, tr[::-1]
    
    def return_parameters(self, x, y):
        D, alignment_score = self.global_alignment(x, y, self.scoring_matrix)
        alignment, transcript = self.traceback(x, y, D, self.scoring_matrix)
        return alignment_score, transcript
        
    
    
