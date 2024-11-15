from Plotting import trim, np

def gradient_nat_convection(T_k_i, T_ambient, k, h):
    value = - h_nat(T_k_i, T_ambient)/ k * (T_k_i - 293/T_ambient)
    return value

def grad_forced_convection(T_k_i, T_ambient, k):
    value = - (11.4 + 20*5.7)/ k * (T_k_i - 293/T_ambient)
    return value

def h_nat(T_s_hat, T_a):
    return 1.31*(T_s_hat*T_a - 293)**(1/3)

def new_imaginary_points(grad, previous_points, h):
    return previous_points + 2*h*grad  

def average_temp(T_k, return_size = False):
    M = trim(T_k)
    rows, cols = len(M), len(M[0].A1)
    
    number_of_points = rows*cols
    total_temp = sum(M.A1)
    avg_temp = total_temp/number_of_points
    if return_size:
        return avg_temp, number_of_points
    else:
        return avg_temp
    
def flip_fins(T):
    n = len(T)
    T_fliped = np.array([np.fliplr(T[-(1+i)]) for i in range(n)])
    return T_fliped

def mirror_array(array, list_of_fins=False):
    if not list_of_fins:
        Array_fist_half = array
        Array_second_half = np.fliplr(array)
        return np.append(Array_fist_half, Array_second_half)
    else:
        Array = [i for i in array]
        Array_second_half = flip_fins(array)
        for i in Array_second_half:
            Array.append(i)
        return np.array(Array)
        


    
    
    
    
    
    
    
    
    
    