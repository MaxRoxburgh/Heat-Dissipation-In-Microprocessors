from Updates import update_system_n_times, update_pr_cs_n_times
from Initializing import full_system_initialisation, np
from Calculators import mirror_array, average_temp
from Plotting import plot_whole_system, combine
from Updates import update_T_matrix_method, update_boundary_convection
import pandas as pd
import time
import matplotlib.pyplot as plt

scaling=6 # refers to points per mm that are being simulated

q = 0.5*1e9 # W/m^3
T_a = 293 # K
K_si = 150 # W/mK
K_ceramic = 230 # W/mK
K_alu = 250 # W/mK
h = (1/scaling)*1e-3 # m    
above_T = 21 # Initial Temp above
bellow_T = 18 # Initial Temp below
n_fins = int(12)  

T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, vol_pr, vol_cs, vol_sinc, vol_fins = full_system_initialisation(above_T, n_fins=n_fins/2, scaling=scaling)
T_0_pr_bellow, T_0_cs_bellow, T_0_hs_bellow, T_0_fn_bellow, T_k_pr_bellow, T_k_cs_bellow, T_k_hs_bellow, T_k_fn_bellow, vol_pr, vol_cs, vol_sinc, vol_fins = full_system_initialisation(bellow_T, n_fins=n_fins/2, scaling=scaling)

try:
    M_inv_pr = np.matrix(pd.read_csv(f'M_inverses/M_inv_pr_s{scaling}', header=None))
    M_inv_cs = np.matrix(pd.read_csv(f'M_inverses/M_inv_cs_s{scaling}', header=None))
    M_inv_hs = np.matrix(pd.read_csv(f'M_inverses/M_inv_hs_34_s{scaling}', header=None))
    M_inv_fn = np.matrix(pd.read_csv(f'M_inverses/M_inv_fn_30_s{scaling}', header=None))
except:
    '''
    NOTE:
        
        if you wish for this to be quicker, swap gaus_jordan for np.linalg.inv
        
        
        
    '''
    from MatrixFunctions import matrix_maker_2, gaus_jordan
    M_pr = matrix_maker_2(T_0_pr_above)
    M_cs = matrix_maker_2(T_0_pr_above)
    M_hs = matrix_maker_2(T_0_pr_above)
    M_fn = matrix_maker_2(np.matrix(T_0_fn_above[0]))
    
    M_inv_pr = gaus_jordan(M_pr)
    M_inv_cs = gaus_jordan(M_cs)
    M_inv_hs = gaus_jordan(M_hs)
    M_inv_fn = gaus_jordan(M_fn)

print('Done with the first section')

above_T = 21 # Initial Temp above
bellow_T = 18 # Initial Temp below
n_fins = int(12)

T_above_list = []
above_T_1 = above_T
diff_above = []
T_bellow_lsit = []
bellow_T_1 = bellow_T
diff_bellow = []
switch = True
nat=True

scalings=[3,4,5,6,7,8,9,10,11]
converge_time_pr = []


for scaling in scalings:
    h = (1/scaling)*1e-3 # m      
    print('\n\n----------------------------------------------\n\n')
    print(f"\n\nRunning for scaling {scaling}")
    
    M_inv_pr = np.matrix(pd.read_csv(f'M_inverses/M_inv_pr_s{scaling}', header=None))
    
    print('Done with the first section')
    
    above_T = 34
    bellow_T = 20
    n_fins = int(12)
    
    T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, vol_pr, vol_cs, vol_sinc, vol_fins = full_system_initialisation(above_T, n_fins=n_fins/2, scaling=scaling)
    T_0_pr_bellow, T_0_cs_bellow, T_0_hs_bellow, T_0_fn_bellow, T_k_pr_bellow, T_k_cs_bellow, T_k_hs_bellow, T_k_fn_bellow, vol_pr, vol_cs, vol_sinc, vol_fins = full_system_initialisation(bellow_T, n_fins=n_fins/2, scaling=scaling)
    
    print('Done with the second section')
    
       
    T_above_list = []
    above_T_1 = above_T
    diff_above = []
    T_bellow_lsit = []
    bellow_T_1 = bellow_T
    diff_bellow = []
    nat=True
    
    t_before = time.time()
    
    while (above_T - bellow_T) > 0.1:
        for n in range(1000):
            T_0_pr_above, T_k_pr_above = update_boundary_convection(T_k_pr_above, 150, h=h), update_T_matrix_method(T_0_pr_above, M_inv_pr, q=q,h=h)
            T_0_pr_bellow, T_k_pr_bellow = update_boundary_convection(T_k_pr_bellow, 150, h=h), update_T_matrix_method(T_0_pr_bellow, M_inv_pr, q=q, h=h)
        
        above_T = average_temp(T_k_pr_above)
        bellow_T = average_temp(T_k_pr_bellow)
        T_above_list.append(above_T)
        diff_above.append(above_T - above_T_1)
        above_T_1 = above_T
        T_bellow_lsit.append(bellow_T)
        diff_bellow.append(bellow_T - bellow_T_1)
        bellow_T_1 = bellow_T
    
    t_after = time.time()
    converge_time_pr.append(t_after-t_before)
    print(f'Time to converge: {int((t_after-t_before)//60)} mins {round((t_after-t_before)%60)} seconds')
    print(f'\n\nAbove T:   {above_T}\nBellow T:  {bellow_T}')
    print('Done with the third section')
    
    ###########################################
    
    #              SAVES DATA                 #
    
    ###########################################
    path = 'Method_testing' 
    
    df = pd.DataFrame(T_k_pr_above)
    df.to_csv(f"{path}/T_k_pr_above_s{scaling}", index=False, header=False)
    
    df = pd.DataFrame(T_k_pr_above)
    df.to_csv(f"{path}/T_k_pr_bellow_s{scaling}", index=False, header=False)
    
    df = pd.DataFrame(np.matrix([T_above_list, diff_above, T_bellow_lsit, diff_bellow]).transpose())
    df.columns = ["T_above", "diff_above", "T_bellow", "diff_bellow"]
    df.to_csv(f'{path}/data_s{scaling}', index=False)

    print('\n\n            SAVED\n\n')
    
    
#%% Another section that runs tests for different lengths of fins


length_of_fins = [32,34,36,38,40,42,44]
nat=True # If you wish to do forced convection, change nat to False
above_T = 3.9
converge_time_length = []

for l in length_of_fins:
    
    print('\n\n----------------------------------------------\n\n')
    print(f"\n\nRunning for; {l} mm fins; scaling {scaling}")
    
    above_T = 3.9
    bellow_T = 2.7
    
    T_above_list = []
    above_T_1 = above_T
    diff_above = []
    T_bellow_lsit = []
    bellow_T_1 = bellow_T
    diff_bellow = []
    
    
    T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, vol_pr, vol_cs, vol_sinc, vol_fins = full_system_initialisation(above_T, n_fins=4, scaling=scaling, new_fin_length=l)
    T_0_pr_bellow, T_0_cs_bellow, T_0_hs_bellow, T_0_fn_bellow, T_k_pr_bellow, T_k_cs_bellow, T_k_hs_bellow, T_k_fn_bellow, vol_pr, vol_cs, vol_sinc, vol_fins = full_system_initialisation(bellow_T, n_fins=4, scaling=scaling, new_fin_length=l)
    
    M = matrix_maker_2(np.matrix(T_0_fn_above[0]))
    M_inv_fn = gaus_jordan(M)
    
    t_before = time.time()
    
    for holding in range(30):
       T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, above_T = update_system_n_times(800, scaling, T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, M_inv_pr, M_inv_cs, M_inv_hs, M_inv_fn, nat)
       T_k_pr_bellow, T_k_cs_bellow, T_k_hs_bellow, T_k_fn_bellow, T_0_pr_bellow, T_0_cs_bellow, T_0_hs_bellow, T_0_fn_bellow, bellow_T = update_system_n_times(800, scaling, T_k_pr_bellow, T_k_cs_bellow, T_k_hs_bellow, T_k_fn_bellow, T_0_pr_bellow, T_0_cs_bellow, T_0_hs_bellow, T_0_fn_bellow, M_inv_pr, M_inv_cs, M_inv_hs, M_inv_fn, nat)
      
       print(f'\n\nAbove T:   {above_T*293}\nBellow T:  {bellow_T*293}') # prints after each iteration
       
       plot_whole_system(T_k_cs_above, T_k_pr_above, mirror_array(T_k_fn_above), T_k_hs_above, min(T_k_cs_bellow))
      
       T_above_list.append(above_T)
       diff_above.append(above_T - above_T_1)
       above_T_1 = above_T
       T_bellow_lsit.append(bellow_T)
       diff_bellow.append(bellow_T - bellow_T_1)
       bellow_T_1 = bellow_T
    
    
    t_after = time.time()
    converge_time_length.append(t_after-t_before)
    print(f'Time to converge: {int((t_after-t_before)//60)} mins {round((t_after-t_before)%60)} seconds')
    print(f'\n\nAbove T:   {above_T*293}\nBellow T:  {bellow_T*293}')
    print('Done with the third section')
    print('Saving...')
        
    path = 'whole_system_FDO4'
    
    df = pd.DataFrame(T_k_pr_above)
    df.to_csv(f"{path}/T_k_pr_above_s{scaling}_fin_length_{l}", index=False, header=False)
    df = pd.DataFrame(T_k_cs_above)
    df.to_csv(f"{path}/T_k_cs_above_s{scaling}_fin_length_{l}", index=False, header=False)
    df = pd.DataFrame(T_k_hs_above)
    df.to_csv(f"{path}/T_k_hs_above_s{scaling}_fin_length_{l}", index=False, header=False)
    for index, fin in enumerate(T_k_fn_above):
        df = pd.DataFrame(fin)
        df.to_csv(f"{path}/T_k_fn_above_s{scaling}_fin_length_{l}_{index}", index=False, header=False)
        
    df = pd.DataFrame(T_k_pr_bellow)
    df.to_csv(f"{path}/T_k_pr_bellow_s{scaling}_fin_length_{l}", index=False, header=False)
    df = pd.DataFrame(T_k_cs_bellow)
    df.to_csv(f"{path}/T_k_cs_bellow_s{scaling}_fin_length_{l}", index=False, header=False)
    df = pd.DataFrame(T_k_hs_bellow)
    df.to_csv(f"{path}/T_k_hs_bellow_s{scaling}_fin_length_{l}", index=False, header=False)
    for index, fin in enumerate(T_k_fn_bellow):
        df = pd.DataFrame(fin)
        df.to_csv(f"{path}/T_k_fn_bellow_s{scaling}_fin_length_{l}_{index}", index=False, header=False)
    
    df = pd.DataFrame(np.matrix([T_above_list, diff_above, T_bellow_lsit, diff_bellow]).transpose())
    df.columns = ["T_above", "diff_above", "T_bellow", "diff_bellow"]
    df.to_csv(f'{path}/data_s{scaling}_fin_length_{l}', index=False)

    print('\n\n            SAVED\n\n')
    