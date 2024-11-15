from Plotting import trim, np
import matplotlib.pyplot as plt
from numba import jit, cuda

def matrix_maker(heat_map, trim_it=True):
    if trim_it:
        heat_map = trim(heat_map)
    
    rows, cols = len(heat_map), len(heat_map[0].A1) 
    #print(cols,rows)
    total_points = cols*rows # length of the vector
    #print(total_points)
    
    matrix = np.identity(total_points)*(-4)
    
    for i in range(total_points): # row
        for j in range(total_points): # col
            if abs(i-j) == 1:
                matrix[i,j] = 1
                
            if i < cols:
                if j == i+cols:
                    matrix[i,j] = 1
            elif i < (total_points-cols):
                if j == i+cols or j == i-cols:
                    matrix[i,j] = 1
            else:
                if j == i-cols:
                    matrix[i,j] = 1
                    
                
        if i > 0:
            if i<(total_points-cols/2):
                if i%cols==0:
                    matrix[i,i-1] = 0
                    
                if (i+1)%cols==0:
                    matrix[i,i+1] = 0
                    
    return matrix

def matrix_maker_2(heat_map, trim_it=True):
    rows, cols = np.shape(heat_map)
    rows, cols = rows-2, cols-2
    M1 = matrix_maker(np.matrix(heat_map), trim_it)
    n = len(M1)
    
    for i in range(n):
        if i//cols >=2 and (rows-i//cols) > 2:
            if i%cols >=2 and (cols-i%cols) > 2:
                M1[i,i] = -5
                M1[i,i+1] = 4/3
                M1[i,i-1] = 4/3
                M1[i,i+2] = -1/12
                M1[i,i-2] = -1/12
                M1[i,i+(cols)] = 4/3
                M1[i,i-(cols)] = 4/3
                M1[i,i+(cols)*2] = -1/12
                M1[i,i-(cols)*2] = -1/12

    return M1

# @jit(target_backend='cuda')   
def gaus_jordan(M, negative_diag=True):
    if negative_diag:
        # this is for our special matricies that we are using in this project
        M = -M
        
    
    n = len(M)    
    I = np.identity(n)

    
    for i in range(n):
        for j in range(n):
            if j!=i:
                scaling = M[j,i]/M[i,i]
                M[j] = M[j] - scaling*M[i]
                I[j] = I[j] - scaling*I[i]
    
    
    for i in range(n):
        scaling = M[i,i]
        M[i] = M[i]/scaling
        I[i] = I[i]/scaling
    
    
    if negative_diag:
        return -I

    
def BC_vector_creator(T_0):
    rows, cols = np.shape(T_0)
    T_0 = T_0.copy()
    
    BC_vec = []
    
    temp = T_0[0].A1
    top_row = temp[1:-1]
    top_row[0] += T_0[1].A1[0].copy()
    top_row[-1] += T_0[1].A1[-1].copy()

    for i in top_row:
        BC_vec.append(i)
    
    for i in range(rows-4):
        mid_row = np.array([0. for i in range(cols-2)])
        mid_row[0] = np.array(T_0[2+i].A1,dtype='float')[0]
        mid_row[-1] = T_0[2+i].A1[-1]
        
        for i in mid_row:
            BC_vec.append(i)

    
    
    temp = T_0[-1].A1
    bottom_row = temp[1:-1]
    bottom_row[0] += T_0[-2].A1[0]
    bottom_row[-1] += T_0[-2].A1[-1]
    
    for i in bottom_row:
            BC_vec.append(i)
            
    return -np.matrix(BC_vec).transpose()

# @jit(target_backend='cuda')   
def update_T_matrix_method(T_0, inv_M, q=0, h=1, inside = False, K_si=150, T_a=293):
    rows, cols = np.shape(T_0)    
    if q == 0:
        T_inside = np.dot(inv_M, BC_vector_creator(T_0)).reshape(rows-2, cols-2)
    else:
        T_inside = np.dot(inv_M, BC_vector_creator(T_0) - q*h**2/ (K_si*T_a)).reshape(rows-2, cols-2)
   
    
    T_k = T_0.copy()
    for i in range(1, rows-1):
        for j in range(1,cols-1):
            T_k[i,j] = T_inside[i-1,j-1]
            
    if inside:
        return T_k, T_inside
    else:
        return T_k
    
def find_D_Tl_Tr(M):
    n = len(M)
    D = np.matrix([[0. for j in range(n)] for i in range(n)] )
    Tl = np.matrix([[0. for j in range(n)] for i in range(n)] )
    Tr = np.matrix([[0. for j in range(n)] for i in range(n)] )
    for i in range(n):
        for j in range(n):
            if i == j:
                D[i,j] = M[i,j]
            if i > j:
                Tl[i,j] = M[i,j] 
            if i < j:
                Tr[i,j] = M[i,j] 
    
    return D, Tl, Tr

def update_T_k_jacobi_matrix(T_k, T_0, M, q=0, h=0, K_si=150, T_a=293):
    BC_vec = BC_vector_creator(np.matrix(T_0.copy()))
    rows, cols = np.shape(T_0)
    T_inside = 1/4 * (-BC_vec + np.dot(M, np.matrix(trim(T_k).A1).transpose()) + q*h**2/ (K_si*T_a))
    T_inside = np.matrix(T_inside).reshape(rows-2,cols-2)
    
    T_k_P1 = T_0.copy()
    for i in range(1, rows-1):
        for j in range(1,cols-1):
            T_k_P1[i,j] = T_inside[i-1,j-1]
    
    return T_k_P1

def update_T_k_GS_matrix(T_k, T_0, M, Tu, q=0, h=0, K_si=150, T_a=293):
    '''
    M in this function is (D+T_l)^-1 aka the inverse of D+T_l
    
    '''
    BC_vec = BC_vector_creator(T_0)
    rows, cols = np.shape(T_0)
    if q != 0:
        q_vec = np.matrix(np.array([q for i in range((rows-2)*(cols-2))])).transpose()
        T_inside = np.dot(M,BC_vec) - np.dot(M, np.dot(Tu, np.matrix(trim(T_k).A1).transpose())) - np.dot(M, q_vec)*h**2 / (K_si*T_a) 
    else:
        T_inside = np.dot(M,BC_vec) - np.dot(M, np.dot(Tu, np.matrix(trim(T_k).A1).transpose()))
    
    T_inside = T_inside.reshape(rows-2,cols-2)
    T_k_P1 = T_0.copy()
    for i in range(1, rows-1):
        for j in range(1,cols-1):
            T_k_P1[i,j] = T_inside[i-1,j-1]
    
    return T_k_P1

# a = np.ones((7,7), dtype='float')
# b = matrix_maker_2(a)


