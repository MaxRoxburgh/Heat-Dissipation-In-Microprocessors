import numpy as np

def initialise_boundary(top, left, right, bottom, T_ambient, corners=[0,0,0,0], Get_T_K=True):
    '''
    Takes in the boundary values and fills all the middle points in with either 0 for T_0 
    or the ambient (starting) temperature to T_k
    adds the corners in either averaged with the points around them or defined prior
    '''
    cols = len(top)+2
    rows = len(left)+2
    
    middle_zeros = np.array([0 for i in top])
    
    if corners == [0,0,0,0]:
        # Average the corners 
        initial_top = np.append(np.append([(top[0]+left[0])/2],top), [(top[-1]+right[0])/2])
        initial_bottom = np.append(np.append([(bottom[0]+left[-1])/2],bottom), [(bottom[-1]+right[-1])/2])
    else:
        # puts in predefined corners
        initial_top = np.append(np.append([corners[0]],top), [corners[1]])
        initial_bottom = np.append(np.append([corners[2]],bottom), [corners[3]])
    
    
    
    middle = np.array(np.append([np.append(np.append([left[i]], middle_zeros), [right[i]]) for i in range(rows-2)],[]))
    
    initial_array = np.append(np.append(initial_top, middle), initial_bottom)
    initial_matrix_T_0 = np.matrix(initial_array.reshape(rows, cols), dtype='float')
    
    if not Get_T_K:
        return initial_matrix_T_0
    
    else:
        middle_ambient = np.array([T_ambient for i in top])
        middle_ambient = np.array(np.append([np.append(np.append([left[i]], middle_ambient), [right[i]]) for i in range(rows-2)],[]))
        initial_array = np.append(np.append(initial_top, middle_ambient), initial_bottom)
        initial_matrix_T_k = np.matrix(initial_array.reshape(rows, cols), dtype='float')
        return initial_matrix_T_0, initial_matrix_T_k
    
def heat_sinc_initialisation(n_fins, ambient, spacing=2, scaling=3):
    '''
    Creates initial boundary and heat map with initial ambient temperature
    this is based on scaling (passed in in mm) 
    scaling is the number of points per mm
    also calculates the solution to poisson equation of array size called M_inv
    also calculates the number of points in the heat sinc body called vol_sinc
    '''
    from MatrixFunctions import matrix_maker, gaus_jordan
    heat_sinc_length = ((n_fins-1)*(spacing+1) + 1)*scaling
    
    top = np.array([ambient for i in range(heat_sinc_length)])
    bottom = top.copy()
    left = np.array([ambient for i in range(4*scaling)])
    right = left.copy()
    
    T_0, T_k = initialise_boundary(top, left, right, bottom, ambient,[ambient,ambient,ambient,ambient])
    
    M = matrix_maker(T_k)
    M_inv = gaus_jordan(M)
    vol_sinc = heat_sinc_length*4*scaling
    
    return T_0, T_k, M_inv, vol_sinc

def full_system_initialisation(ambient, n_fins, scaling=3, get_heat_sinc=True, new_fin_length=30, spacing = 2):
    '''
    initialises the temperatures of boundaries and internal grids for each section of the system
    n_fins passed in is the number wanted to create, in our case we only create half of them as
    we use symetry to copy the elements in our code
    '''
    # set up the processor
    
    vol_pr = 14*scaling**2
    
    top = np.array([ambient for i in range(14*scaling)])
    bottom = top.copy()
    left = np.array([ambient for i in range(scaling)])
    right = left.copy()
    
    T_0_pr, T_k_pr = initialise_boundary(top, left, right, bottom, ambient,[ambient,ambient,ambient,ambient])
    
    # set up the case
    
    vol_cs = 20*2*scaling**2
    
    top = np.array([ambient for i in range(20*scaling)])
    bottom = top.copy()
    left = np.array([ambient for i in range(2*scaling)])
    right = left.copy()
    
    T_0_cs, T_k_cs = initialise_boundary(top, left, right, bottom, ambient,[ambient,ambient,ambient,ambient])
    
    
    # set up the main body of the heat sinc (34x4)
    
    if get_heat_sinc:
        
        sinc_length = ((n_fins*2-1)*(1+spacing) + 1)*scaling
        vol_sinc = sinc_length*4*scaling
        
        top = np.array([ambient for i in range(int(sinc_length))])
        bottom = top.copy()
        left = np.array([ambient for i in range(4*scaling)])
        right = left.copy()
        
        T_0_hs, T_k_hs = initialise_boundary(top, left, right, bottom, ambient,[ambient,ambient,ambient,ambient])
    
    # set up fins for heat sinc
    
    fin_length = new_fin_length*scaling
    
    top = np.array([ambient for i in range(scaling)])
    bottom = top.copy()
    left = np.array([ambient for i in range(fin_length)])
    right = left.copy()
    
    
    vol_fins = fin_length*scaling**2 * n_fins
    fin_0, fin_k = initialise_boundary(top, left, right, bottom, ambient,[ambient,ambient,ambient,ambient])
    T_0_fn = np.array([fin_0 for i in range(int(n_fins))])
    T_k_fn = np.array([fin_k for i in range(int(n_fins))])
    
    if not get_heat_sinc:
        return T_0_pr, T_0_cs, T_0_fn, T_k_pr, T_k_cs, T_k_fn, vol_pr, vol_cs, vol_fins
    
    return T_0_pr, T_0_cs, T_0_hs, T_0_fn, T_k_pr, T_k_cs, T_k_hs, T_k_fn, vol_pr, vol_cs, vol_sinc, vol_fins