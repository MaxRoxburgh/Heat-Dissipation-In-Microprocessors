import matplotlib.pyplot as plt
import numpy as np


def plot(T, cbar = False, vbar = True, invis = False, show=True):
    fig, ax = plt.subplots()
    im = ax.imshow(T, cmap='plasma')
    if invis:
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)
    if cbar:
        #fig.suptitle('Temperature map')
        if vbar:
            cbar = ax.figure.colorbar(im, orientation='vertical')
            cbar.set_label('Temperature Scale', loc='center')
        else:
            cbar = ax.figure.colorbar(im, orientation='horizontal')
            cbar.ax.set_xlabel('Temperature Scale')
    if show:
        plt.show()  
        
def trim(M):
    # takes the edge rows off a matrix
    # this is used as the T_k matricies that are generated contain the imaginary points from the neuman BC
    rows = len(M)
    cols = len(M[0].A1)
    
    M_list = M.A1[cols:-cols]
    M = np.matrix(M_list.reshape(rows-2, cols)).transpose()
    M_list = M.A1[rows-2:-(rows-2)]
    M = np.matrix(M_list.reshape(cols-2, rows-2)).transpose()
    
    return M    

def combine(M1, M2, Temp, Trim = True):
    # This function is to combine the two matricies, M1 on top and M2 on bottom
    # It assumes that the matricies contain the boundary conditions and therefore trims them both
    # It will add k to anywhere that is dead-space
    # Assumes that M1 has more cols than M2
    # print(M1)
    
    if Trim:
        M1 = trim(M1)
        M2 = trim(M2)
    
    M1_rows, M1_cols = len(M1), len(M1[0].A1)
    M2_rows, M2_cols = len(M2), len(M2[0].A1)
    
    total_rows = M1_rows + M2_rows
    total_cols = M1_cols
    
    if M1_cols < M2_cols:
        raise Exception('M1 should be the bigger matrix, but it has less columns')
    
    cols_diff = M1_cols - M2_cols
    
    total_array = M1.A1
    
    for row in M2:
        ks = [Temp for i in range(int(cols_diff/2))]
        middle = np.append(np.append(ks, row),ks)
        total_array = np.append(total_array, middle)
        
    return np.matrix(total_array).reshape(total_rows, total_cols)

def combine_fins(fins, spacing=20, ambient=0):
    test_fin = trim(np.matrix(fins[0]))
    rows, cols = len(test_fin), len(test_fin[0].A1)
    number_of_fins = len(fins)
    
    total_matrix = np.array([])
    for r in range(rows):
        for f in range(number_of_fins-1):
            fin = trim(np.matrix(fins[f]))
            total_matrix = np.append(np.append(total_matrix, fin[r]),[ambient for i in range(spacing)])
        fin = trim(np.matrix(fins[-1]))
        total_matrix = np.append(total_matrix, fin[r])
    
    total_matrix = np.matrix(total_matrix).reshape(rows, (cols*number_of_fins+spacing*(number_of_fins-1)))
    
    return total_matrix  

def combine_fins_wo_trim(fins, spacing=20, ambient=0):
    test_fin = np.matrix(fins[0])
    rows, cols = len(test_fin), len(test_fin[0].A1)
    number_of_fins = len(fins)
    
    total_matrix = np.array([])
    for r in range(rows):
        for f in range(number_of_fins-1):
            fin = np.matrix(fins[f])
            total_matrix = np.append(np.append(total_matrix, fin[r]),[ambient for i in range(spacing)])
        fin = np.matrix(fins[-1])
        total_matrix = np.append(total_matrix, fin[r])
    
    total_matrix = np.matrix(total_matrix).reshape(rows, (cols*number_of_fins+spacing*(number_of_fins-1)))
    
    return total_matrix  

def plot_whole_system(T_k_cs, T_k_pr, T_k_fn, T_k_hs, ambient, spacing=20, plotting=True):
    pr_cs = combine(T_k_cs, T_k_pr, ambient)
    fins_combined = combine_fins(T_k_fn, spacing, ambient)
    hs_wFins = combine(fins_combined, trim(T_k_hs), ambient, False)
    total_system = combine(hs_wFins, pr_cs, ambient, False)
    # print(min(total_system.A1))
    if plotting:
        plot(total_system, True, True, True)
    else:
        return total_system
    
def take_outside(M):
    from Initializing import initialise_boundary
    top = M[0].A1[1:-1]
    bottom = M[-1].A1[1:-1]
    M = M.transpose()
    left = M[0].A1[1:-1]
    right = M[-1].A1[1:-1]
    
    M = initialise_boundary(top, left, right, bottom, 0, Get_T_K=False)
    return M
    
def temp_calc(T_k_pr,T_k_cs,T_k_hs,T_k_fn, avT_k_1, vol_fins, vol_cs, vol_pr, vol_sinc, i, every=1000):
    from Calculators import average_temp, mirror_array
    
    average_T_pr = average_temp(T_k_pr)
    average_T_cs = average_temp(T_k_cs)
    average_T_fn = sum([average_temp(np.matrix(T_k)) for T_k in T_k_fn])/len(T_k_fn)
    average_T_hs = average_temp(T_k_cs)
    
    average_fin_hs = (2*vol_fins*average_T_fn + vol_sinc*average_T_hs)/(2*vol_fins+vol_sinc)
    
    print('T:     pr     |     cs     |      hs      |\n', average_T_pr, average_T_cs, average_fin_hs)
    
    ambient = min(T_k_fn[0][0])
    plot_whole_system(T_k_cs, T_k_pr, mirror_array(T_k_fn, True), T_k_hs, ambient, 6)
    
    avT_k = (2*vol_fins*average_T_fn + vol_sinc*average_T_hs + vol_cs*average_T_cs + vol_pr* average_T_pr)/(2*vol_fins+vol_sinc+vol_cs+vol_pr)
    diff = avT_k-avT_k_1
    print(f'\nAverage T, difference from {every} runs ago: run {i}\n',avT_k, diff)
    return avT_k, diff
  
def save(T_k_pr_update, T_k_cs_update, T_k_hs_update, T_k_fn_update, run_title, path="Optimised_T_k"):
    import pandas as pd
    df = pd.DataFrame(T_k_pr_update)
    df.to_csv(f"{path}/T_k_pr_{run_title}", index=False, header=False)
    df = pd.DataFrame(T_k_cs_update)
    df.to_csv(f"{path}/T_k_cs_{run_title}", index=False, header=False)
    df = pd.DataFrame(T_k_hs_update)
    df.to_csv(f"{path}/T_k_hs_{run_title}", index=False, header=False)
    for i, T_k in enumerate(T_k_fn_update):
        df = pd.DataFrame(T_k)
        df.to_csv(f"{path}/T_k_fn_{i}_{run_title}", index=False, header=False)

def plot_3d_matrix(matrix, elev, azim, ax, number, fig):
    rows, cols = matrix.shape
    x = np.arange(cols)
    y = np.arange(rows)
    x, y = np.meshgrid(x, y)
    z = matrix

    # fig = plt.figure()
    ax = fig.add_subplot(number, projection='3d')
    ax.plot_surface(x, y, z, cmap='plasma')

    # ax.set_xlabel('X-axis')
    # ax.set_ylabel('Y-axis')
    ax.set_zlabel('Temperature [K]')
    ax.view_init(elev=elev, azim=azim)
    ax.axis('off')
#
    plt.show()    
    
def add_boundary(M, T, n):
    r,c = M.shape
    n-=1    
    M_new = np.matrix(np.ones((r+2,c+2), dtype='float'))*T
    for i in range(1,r-1):
        for j in range(1,c-1):
            M_new[i,j]=M[i,j]
            
    if n == 1:
        return M
    else:
        return add_boundary(np.matrix(M_new), T, n)
    