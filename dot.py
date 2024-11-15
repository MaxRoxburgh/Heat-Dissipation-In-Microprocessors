import numpy as np
def dot(m, v):
    V = [i for i in v.A1]
    M = [i for i in m.A1]
    n = len(V)
    
    new = []
    for i in range(n):
        row_multiple = sum([V[j]*M[j+n*i] for j in range(n)])
        new.append(row_multiple)
    
    return np.matrix(new).transpose()
