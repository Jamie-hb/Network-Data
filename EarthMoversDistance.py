import numpy as np
import ot

def d_LCSS(I1, I2):
    # computes the Lowest Common SubSequence distance between I1 amd I2
    # I1, I2 : one-dimensional array of numbers 
    m = len(I1)
    n = len(I2)
    
    d = np.zeros((m+1,n+1))
    
    d[:,0] = list(range(0,m+1))
    d[0,:] = list(range(0,n+1))
    
    for j in range(n):
        for i in range(m):
            if I1[i] == I2[j]:
                substitutioncost = 0
            else:
                substitutioncost = 2
            
            d[i+1,j+1] = min(d[i,j+1] + 1, d[i+1,j] + 1, d[i,j] + substitutioncost)
    
    return(d[m,n])

def dbar_LCSS(I1, I2):
    # returns the Normalised (via the Steinhaus/biotope transform) Lowest
    # Common SubSequence distance between I1 and I2
    # I1, I2 : one-dimensional arrays of numbers
    m = len(I1)
    n = len(I2)
    d = d_LCSS(I1, I2)
    delta = (m + n - d)/2 # d_LCSS = m + n - 2δ => δ = ...
    
    return(1 - (delta / (m + n - delta)))

def INS_2_EMD(S1, S2):
    # Input : two sequences of interaction networks S1, S2, each stored as a list of lists, 
    # with each (sub-)list being an single interaction network 
    # Output : a list of 4 elements
    # 1) the unique entries of S1, stored as a list, with the order of the list being
    #    the order in which the interaction networks appear in S1
    # 2) the unique entries of S2, stored as a list, with the order of the list being
    #    the order in which the interaction networks appear in S2
    # 3) list of weights for S1, with entry i being the amount of times the i'th row of (1) appears in S1
    # 4) list of weights for S2, with entry i being the amount of times the i'th row of (2) appears in S2
    
    L1 = len(S1)
    L2 = len(S2)
    
    res1 = list()
    res2 = list()
    res3 = list()
    res4 = list()
    
    for i in range(L1):
        t = S1[i]
        if t in res1:
            res3[res1.index(t)] += 1
        else:
            res1.append(t)
            res3.append(1)
    
    for j in range(L2):
        t = S2[j]
        if t in res2:
            res4[res2.index(t)] += 1
        else:
            res2.append(t)
            res4.append(1)
    
    return([res1, res2, res3, res4])

    
def EMD(S1, S2, d):
    # Takes two sequences of interaction networks, S1 and S2, and returns the
    # Earth Mover's Distance between them, using d as a metric on the space of interaction 
    # networks.
    
    t = INS_2_EMD(S1, S2)
    m = len(t[0])
    n = len(t[1])
    D = np.zeros((m,n))
    
    for i in range(m):
        for j in range(n):
            r = d(t[0][i], t[1][j])
            D[i,j] = r
    
    w1 = t[2]
    w2 = t[3]        
        
    F = ot.emd(w1, w2, D)
    
    return( np.sum(F * D) / np.sum(F) )
    
    








