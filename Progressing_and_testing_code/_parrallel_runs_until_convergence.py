from Updates import update_system_n_times, update_pr_cs_n_times
from Initializing import full_system_initialisation, np
from Calculators import mirror_array
from Plotting import plot_whole_system, combine
import pandas as pd
import time
import matplotlib.pyplot as plt

def change(M, addition, is_list=False):
   
    if is_list:
        M_list = M.copy()
        M_new_list = []
        for m in M_list:
            r,c=np.shape(m)
            vec = np.matrix(m).A1
            for i in range(len(vec)):
                if vec[i] != 0:
                    vec[i] += addition
            M_new_list.append(np.matrix(vec).reshape(r,c))
       
        return np.array(M_new_list)
           
   
    else:
        m = M.copy()
        r,c=np.shape(m)
        vec = np.matrix(m).A1
        for i in range(len(vec)):
            if vec[i] != 0:
                vec[i] += addition
       
        return np.matrix(vec).reshape(r,c)


#%%
scaling=6

q = 0.5*1e9 # W/m^3
T_a = 293 # K
K_si = 150 # W/mK
K_ceramic = 230 # W/mK
K_alu = 250 # W/mK
h = (1/scaling)*1e-3 # m      

M_inv_pr = np.matrix(pd.read_csv(f'M_inverses/M_inv_pr_s{scaling}', header=None))
M_inv_cs = np.matrix(pd.read_csv(f'M_inverses/M_inv_cs_s{scaling}', header=None))
M_inv_hs = np.matrix(pd.read_csv(f'M_inverses/M_inv_hs_34_s{scaling}', header=None))
M_inv_fn = np.matrix(pd.read_csv(f'M_inverses/M_inv_fn_30_s{scaling}', header=None))

print('Done with the first section')

above_T = 21
bellow_T = 18
n_fins = int(12)


T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, vol_pr, vol_cs, vol_sinc, vol_fins = full_system_initialisation(above_T, n_fins=n_fins/2, scaling=scaling)
T_0_pr_bellow, T_0_cs_bellow, T_0_hs_bellow, T_0_fn_bellow, T_k_pr_bellow, T_k_cs_bellow, T_k_hs_bellow, T_k_fn_bellow, vol_pr, vol_cs, vol_sinc, vol_fins = full_system_initialisation(bellow_T, n_fins=n_fins/2, scaling=scaling)


T_above_list = []
above_T_1 = above_T
diff_above = []
T_bellow_lsit = []
bellow_T_1 = bellow_T
diff_bellow = []
switch = True
nat=True
#%%
# add = -2.
# T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, above_T = change(T_k_pr_above, add), change(T_k_cs_above, add), change(T_k_hs_above, add), change(T_k_fn_above, add, True), change(T_0_pr_above, add), change(T_0_cs_above, add), change(T_0_hs_above, add), change(T_0_fn_above, add, True), above_T+add


while (above_T - bellow_T) > 0.1:
       switch = not switch
       T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, above_T = update_system_n_times(500, scaling, T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, M_inv_pr, M_inv_cs, M_inv_hs, M_inv_fn, nat)
       T_k_pr_bellow, T_k_cs_bellow, T_k_hs_bellow, T_k_fn_bellow, T_0_pr_bellow, T_0_cs_bellow, T_0_hs_bellow, T_0_fn_bellow, bellow_T = update_system_n_times(500, scaling, T_k_pr_bellow, T_k_cs_bellow, T_k_hs_bellow, T_k_fn_bellow, T_0_pr_bellow, T_0_cs_bellow, T_0_hs_bellow, T_0_fn_bellow, M_inv_pr, M_inv_cs, M_inv_hs, M_inv_fn, nat)
      
       # print(f'\n\nAbove T:   {above_T}\nBellow T:  {bellow_T}')
      
       T_above_list.append(above_T)
       diff_above.append(above_T - above_T_1)
       above_T_1 = above_T
       T_bellow_lsit.append(bellow_T)
       diff_bellow.append(bellow_T - bellow_T_1)
       bellow_T_1 = bellow_T
       
        # if switch:
        #     plot_whole_system(T_k_cs_above, T_k_pr_above, mirror_array(T_k_fn_above, True), T_k_hs_above, above_T, spacing=2*scaling)
        #     plot_whole_system(T_k_cs_bellow, T_k_pr_bellow, mirror_array(T_k_fn_bellow, True), T_k_hs_bellow, bellow_T, spacing=2*scaling)
        #     plt.plot(T_above_list[-5:], abs(np.array(diff_above)[-5:]))
        #     plt.plot(T_bellow_lsit[-5:], diff_bellow[-5:])
        #     plt.yscale('log')
print('Done with the third section')

#%%
q = 0.5*1e9 # W/m^3
T_a = 293 # K
K_si = 150 # W/mK
K_ceramic = 230 # W/mK
K_alu = 250 # W/mK

scalings=[3,4,5,6,7,8,9]
converge_time = []


for scaling in scalings:
    h = (1/scaling)*1e-3 # m      
    print('\n\n----------------------------------------------\n\n')
    print(f"\n\nRunning for scaling {scaling}")
    
    M_inv_pr = np.matrix(pd.read_csv(f'M_inverses/M_inv_pr_s{scaling}_FDO4', header=None))
    M_inv_cs = np.matrix(pd.read_csv(f'M_inverses/M_inv_cs_s{scaling}_FDO4', header=None))
    # M_inv_hs = np.matrix(pd.read_csv(f'M_inverses/M_inv_hs_34_s{scaling}', header=None))
    # M_inv_fn = np.matrix(pd.read_csv(f'M_inverses/M_inv_fn_30_s{scaling}', header=None))
    
    print('Done with the first section')
    
    above_T = 16+scaling
    bellow_T = 12+scaling
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
    # switch = True
    nat=True
    
    t_before = time.time()
    
    while (above_T - bellow_T) > 0.1:
        T_k_pr_above, T_k_cs_above, T_0_pr_above, T_0_cs_above, above_T = update_pr_cs_n_times(400*scaling, scaling, T_k_pr_above, T_k_cs_above, T_0_pr_above, T_0_cs_above, M_inv_pr, M_inv_cs, nat)
        T_k_pr_bellow, T_k_cs_bellow, T_0_pr_bellow, T_0_cs_bellow, bellow_T = update_pr_cs_n_times(400*scaling, scaling, T_k_pr_bellow, T_k_cs_bellow, T_0_pr_bellow, T_0_cs_bellow, M_inv_pr, M_inv_cs, nat)
       
        # print(f'\n\nAbove T:   {above_T}\nBellow T:  {bellow_T}')
       
        T_above_list.append(above_T)
        diff_above.append(above_T - above_T_1)
        above_T_1 = above_T
        T_bellow_lsit.append(bellow_T)
        diff_bellow.append(bellow_T - bellow_T_1)
        bellow_T_1 = bellow_T
    
    t_after = time.time()
    converge_time.append(t_after-t_before)
    print(f'Time to converge: {int((t_after-t_before)//60)} mins {round((t_after-t_before)%60)} seconds')
    print(f'\n\nAbove T:   {above_T}\nBellow T:  {bellow_T}')
    print('Done with the third section')
    
    path = f'pr_cs_data_FDO4/scaling_{scaling}'
    
    df = pd.DataFrame(T_k_pr_above)
    df.to_csv(f"{path}/T_k_pr_above_s{scaling}", index=False, header=False)
    df = pd.DataFrame(T_k_cs_above)
    df.to_csv(f"{path}/T_k_cs_above_s{scaling}", index=False, header=False)
    
    df = pd.DataFrame(T_k_pr_above)
    df.to_csv(f"{path}/T_k_pr_bellow_s{scaling}", index=False, header=False)
    df = pd.DataFrame(T_k_cs_above)
    df.to_csv(f"{path}/T_k_cs_bellow_s{scaling}", index=False, header=False)
    
    df = pd.DataFrame(np.matrix([T_above_list, diff_above, T_bellow_lsit, diff_bellow]).transpose())
    df.columns = ["T_above", "diff_above", "T_bellow", "diff_bellow"]
    df.to_csv(f'{path}/data_s{scaling}', index=False)

    print('\n\n            SAVED\n\n')


### 9.337260521064849 for converged value s10
# converge times for scalings = [3,4,5,6,7,8,9] - normal FD scheme
# array([ 208.17692184,  415.98871136,  708.64798927, 1162.86647415, 1894.4648818969727, 2932.9399399757385, 4367.7890248298645])

# convergence


#%%
from Updates import update_T_matrix_method, update_boundary_convection
from Calculators import average_temp

q = 0.5*1e9 # W/m^3
T_a = 293 # K
K_si = 150 # W/mK
K_ceramic = 230 # W/mK
K_alu = 250 # W/mK

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
    # switch = True
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
    
    path = f'pr_data/scaling_{scaling}'
    
    df = pd.DataFrame(T_k_pr_above)
    df.to_csv(f"{path}/T_k_pr_above_s{scaling}", index=False, header=False)
    
    df = pd.DataFrame(T_k_pr_above)
    df.to_csv(f"{path}/T_k_pr_bellow_s{scaling}", index=False, header=False)
    
    df = pd.DataFrame(np.matrix([T_above_list, diff_above, T_bellow_lsit, diff_bellow]).transpose())
    df.columns = ["T_above", "diff_above", "T_bellow", "diff_bellow"]
    df.to_csv(f'{path}/data_s{scaling}', index=False)

    print('\n\n            SAVED\n\n')

#converge_time_pr = [90.00240731239319,154.2232186794281, 239.1755564212799, 362.6052348613739, 537.7145259380341, 699.7213654518127, 974.6098198890686, 1369.6851847171783, 1828.2720837593079]
#%%
converge_time_pr_FDO4 = []
scalings=[3,4,5,6,7,8,9]
# scalings = [4]
for scaling in scalings:
    h = (1/scaling)*1e-3 # m      
    print('\n\n----------------------------------------------\n\n')
    print(f"\n\nRunning for scaling {scaling}")
    
    M_inv_pr = np.matrix(pd.read_csv(f'M_inverses/M_inv_pr_s{scaling}_FDO4', header=None))
    
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
    # switch = True
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
    converge_time_pr_FDO4.append(t_after-t_before)
    print(f'Time to converge: {int((t_after-t_before)//60)} mins {round((t_after-t_before)%60)} seconds')
    print(f'\n\nAbove T:   {above_T}\nBellow T:  {bellow_T}')
    print('Done with the third section')
    
    path = f'pr_data_FDO4/scaling_{scaling}'
    
    df = pd.DataFrame(T_k_pr_above)
    df.to_csv(f"{path}/T_k_pr_above_s{scaling}", index=False, header=False)
    
    df = pd.DataFrame(T_k_pr_above)
    df.to_csv(f"{path}/T_k_pr_bellow_s{scaling}", index=False, header=False)
    
    df = pd.DataFrame(np.matrix([T_above_list, diff_above, T_bellow_lsit, diff_bellow]).transpose())
    df.columns = ["T_above", "diff_above", "T_bellow", "diff_bellow"]
    df.to_csv(f'{path}/data_s{scaling}', index=False)

    print('\n\n            SAVED\n\n')

# converge time for FDO4
#array([  91.46322227,  161.5895288 ,  246.50573492,  367.20968032, 524.53289032,  743.82717538, 1020.28084302])

#%% n_fins test
# from Plotting import trim

scaling = 6

q = 0.5*1e9 # W/m^3
T_a = 293 # K
K_si = 150 # W/mK
K_ceramic = 230 # W/mK
K_alu = 250 # W/mK
h = (1/scaling)*1e-3 # m    


print('Done with the first section')

nat=False
n_fins_list = [8,10,12,14,16]
# n_fins_list = [12,14,16]
M_inv_pr = np.matrix(pd.read_csv(f'M_inverses/M_inv_pr_s{scaling}_FDO4', header=None))
M_inv_cs = np.matrix(pd.read_csv(f'M_inverses/M_inv_cs_s{scaling}_FDO4', header=None))
M_inv_fn = np.matrix(pd.read_csv(f'M_inverses/M_inv_fn_30_s{scaling}_FDO4', header=None))

converge_time = []

for n in n_fins_list:
    hs_length = (n-1)*3+1
    # print(hs_length*6)
    
    print('\n\n----------------------------------------------\n\n')
    print(f"\n\nRunning for; {n} fins; scaling {scaling}")
    

    # print(trim(T_0_hs_above).shape, trim(T_0_hs_bellow).shape)
    
    M_inv_hs = np.matrix(pd.read_csv(f'M_inverses/M_inv_hs_{hs_length}_s{scaling}_FDO4', header=None))
    
    above_T = 1.8
    bellow_T = 1.1
    
    # if n == 8:
    #     above_T = 3.95
    #     bellow_T = 3.8
    # elif n ==10:
    #     above_T = 3.5
    #     bellow_T = 3.4
    # elif n ==12:
    #     above_T = 3.2
    #     bellow_T = 3.1
    # elif n == 14:
    #     above_T = 2.95
    #     bellow_T = 2.88
    # else:
    #     above_T = 2.72
    #     bellow_T = 2.79
        
    T_above_list = []
    T_pr_above_list = []
    above_T_1 = above_T
    diff_above = []
    T_bellow_lsit = []
    T_pr_below_list = []
    bellow_T_1 = bellow_T
    diff_bellow = []
    t_before = time.time()
    
    T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, vol_pr, vol_cs, vol_sinc, vol_fins = full_system_initialisation(above_T, n_fins=n/2, scaling=scaling)
    T_0_pr_bellow, T_0_cs_bellow, T_0_hs_bellow, T_0_fn_bellow, T_k_pr_bellow, T_k_cs_bellow, T_k_hs_bellow, T_k_fn_bellow, vol_pr, vol_cs, vol_sinc, vol_fins = full_system_initialisation(bellow_T, n_fins=n/2, scaling=scaling)
    
    # add = -0.5
    # T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, above_T = change(T_k_pr_above, add), change(T_k_cs_above, add), change(T_k_hs_above, add), change(T_k_fn_above, add, True), change(T_0_pr_above, add), change(T_0_cs_above, add), change(T_0_hs_above, add), change(T_0_fn_above, add, True), above_T+add
    
    t_before = time.time()
    
    for m in range(25):
       T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, above_T = update_system_n_times(500, scaling, T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, M_inv_pr, M_inv_cs, M_inv_hs, M_inv_fn, nat)
       T_k_pr_bellow, T_k_cs_bellow, T_k_hs_bellow, T_k_fn_bellow, T_0_pr_bellow, T_0_cs_bellow, T_0_hs_bellow, T_0_fn_bellow, bellow_T = update_system_n_times(500, scaling, T_k_pr_bellow, T_k_cs_bellow, T_k_hs_bellow, T_k_fn_bellow, T_0_pr_bellow, T_0_cs_bellow, T_0_hs_bellow, T_0_fn_bellow, M_inv_pr, M_inv_cs, M_inv_hs, M_inv_fn, nat)
      
       # print(f'\n\nAbove T:   {above_T}\nBellow T:  {bellow_T}')
      
       T_pr_above = average_temp(T_k_pr_above)
       T_pr_below = average_temp(T_k_pr_bellow)  
      
       T_above_list.append(above_T)
       T_pr_above_list.append(T_pr_above)
       diff_above.append(above_T - above_T_1)
       above_T_1 = above_T
       T_bellow_lsit.append(bellow_T)
       T_pr_below_list.append(T_pr_below)
       diff_bellow.append(bellow_T - bellow_T_1)
       bellow_T_1 = bellow_T
       
    
    
    t_after = time.time()
    converge_time.append(t_after-t_before)
    print(f'Time to converge: {int((t_after-t_before)//60)} mins {round((t_after-t_before)%60)} seconds')
    print(f'\n\nAbove T:   {above_T}\nBellow T:  {bellow_T}')
    print('Done with the third section')
    print('Saving...')
        
    path = 'whole_system_FDO4/n_fins_test_forced'
    
    df = pd.DataFrame(T_k_pr_above)
    df.to_csv(f"{path}/T_k_pr_above_s{scaling}_nfin_{n}", index=False, header=False)
    df = pd.DataFrame(T_k_cs_above)
    # df.to_csv(f"{path}/T_k_cs_above_s{scaling}_nfin_{n}", index=False, header=False)
    # df = pd.DataFrame(T_k_hs_above)
    # df.to_csv(f"{path}/T_k_hs_above_s{scaling}_nfin_{n}", index=False, header=False)
    # for index, fin in enumerate(T_k_fn_above):
    #     df = pd.DataFrame(fin)
    #     df.to_csv(f"{path}/T_k_fn_above_s{scaling}_nfin_{n}_{index}", index=False, header=False)
        
    df = pd.DataFrame(T_k_pr_bellow)
    df.to_csv(f"{path}/T_k_pr_bellow_s{scaling}_nfin_{n}", index=False, header=False)
    df = pd.DataFrame(T_k_cs_bellow)
    # df.to_csv(f"{path}/T_k_cs_bellow_s{scaling}_nfin_{n}", index=False, header=False)
    # df = pd.DataFrame(T_k_hs_bellow)
    # df.to_csv(f"{path}/T_k_hs_bellow_s{scaling}_nfin_{n}", index=False, header=False)
    # for index, fin in enumerate(T_k_fn_bellow):
    #     df = pd.DataFrame(fin)
    #     df.to_csv(f"{path}/T_k_fn_bellow_s{scaling}_nfin_{n}_{index}_2", index=False, header=False)
    
    df = pd.DataFrame(np.matrix([T_above_list, diff_above, T_bellow_lsit, diff_bellow, T_pr_above_list, T_pr_below_list]).transpose())
    df.columns = ["T_above", "diff_above", "T_bellow", "diff_bellow", 'T_pr_above', 'T_pr_below']
    df.to_csv(f'{path}/data_s{scaling}_nfins_{n}', index=False)

    print('\n\n            SAVED\n\n')
    

#%%   
    
scaling = 6

q = 0.5*1e9 # W/m^3
T_a = 293 # K
K_si = 150 # W/mK
K_ceramic = 230 # W/mK
K_alu = 250 # W/mK
h = (1/scaling)*1e-3 # m      

print('Done with the first section')

nat=True
length_of_fins = [32,34,36,38,40,42,44]

M_inv_pr = np.matrix(pd.read_csv(f'M_inverses/M_inv_pr_s{scaling}_FDO4', header=None))
M_inv_cs = np.matrix(pd.read_csv(f'M_inverses/M_inv_cs_s{scaling}_FDO4', header=None))
M_inv_hs = np.matrix(pd.read_csv(f'M_inverses/M_inv_hs_22_s{scaling}_FDO4', header=None))
above_T = 3.9

converge_time_length = []

for l in length_of_fins:
    # print(hs_length*6)
    
    print('\n\n----------------------------------------------\n\n')
    print(f"\n\nRunning for; {l} mm fins; scaling {scaling}")
    

    # print(trim(T_0_hs_above).shape, trim(T_0_hs_bellow).shape)
    
    M_inv_fn = np.matrix(pd.read_csv(f'M_inverses/M_inv_fn_{l}_s{scaling}_FDO4', header=None))
    bellow_T = 2.7
    
    # if l == 32:
    #     above_T = 3.83
    #     bellow_T = 3.78
    # elif l ==34:
    #     above_T = 3.71
    #     bellow_T = 3.65
    # elif l ==36:
    #     above_T = 3.61
    #     bellow_T = 3.55
    # elif l == 38:
    #     above_T = 3.53
    #     bellow_T = 3.48
    # elif l == 40:
    #     above_T = 3.42
    #     bellow_T = 3.38
    # elif l == 42:
    #     above_T = 3.61
    #     bellow_T = 3.55
    # elif l == 44:
    #     above_T = 3.275
    #     bellow_T = 3.21
    
    
    T_above_list = []
    above_T_1 = above_T
    diff_above = []
    T_bellow_lsit = []
    bellow_T_1 = bellow_T
    diff_bellow = []
    
    
    T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, vol_pr, vol_cs, vol_sinc, vol_fins = full_system_initialisation(above_T, n_fins=4, scaling=scaling, new_fin_length=l)
    T_0_pr_bellow, T_0_cs_bellow, T_0_hs_bellow, T_0_fn_bellow, T_k_pr_bellow, T_k_cs_bellow, T_k_hs_bellow, T_k_fn_bellow, vol_pr, vol_cs, vol_sinc, vol_fins = full_system_initialisation(bellow_T, n_fins=4, scaling=scaling, new_fin_length=l)
    
    t_before = time.time()
    
    for holding in range(15):
       T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, above_T = update_system_n_times(800, scaling, T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, M_inv_pr, M_inv_cs, M_inv_hs, M_inv_fn, nat)
       T_k_pr_bellow, T_k_cs_bellow, T_k_hs_bellow, T_k_fn_bellow, T_0_pr_bellow, T_0_cs_bellow, T_0_hs_bellow, T_0_fn_bellow, bellow_T = update_system_n_times(800, scaling, T_k_pr_bellow, T_k_cs_bellow, T_k_hs_bellow, T_k_fn_bellow, T_0_pr_bellow, T_0_cs_bellow, T_0_hs_bellow, T_0_fn_bellow, M_inv_pr, M_inv_cs, M_inv_hs, M_inv_fn, nat)
      
       print(f'\n\nAbove T:   {above_T}\nBellow T:  {bellow_T}')
      
       T_above_list.append(above_T)
       diff_above.append(above_T - above_T_1)
       above_T_1 = above_T
       T_bellow_lsit.append(bellow_T)
       diff_bellow.append(bellow_T - bellow_T_1)
       bellow_T_1 = bellow_T
    
    
    t_after = time.time()
    converge_time_length.append(t_after-t_before)
    print(f'Time to converge: {int((t_after-t_before)//60)} mins {round((t_after-t_before)%60)} seconds')
    print(f'\n\nAbove T:   {above_T}\nBellow T:  {bellow_T}')
    print('Done with the third section')
    print('Saving...')
        
    path = 'whole_system_FDO4/fin_length_test'
    
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
    
    
#%%   
from Plotting import take_outside

scaling = 6

q = 0.5*1e9 # W/m^3
T_a = 293 # K
K_si = 150 # W/mK
K_ceramic = 230 # W/mK
K_alu = 250 # W/mK
h = (1/scaling)*1e-3 # m      

print('Done with the first section')

nat=False
length_of_fins = [30,32,34,36,38,40,42,44]
length_of_fins = [30,34]
M_inv_pr = np.matrix(pd.read_csv(f'M_inverses/M_inv_pr_s{scaling}_FDO4', header=None))
M_inv_cs = np.matrix(pd.read_csv(f'M_inverses/M_inv_cs_s{scaling}_FDO4', header=None))
M_inv_hs = np.matrix(pd.read_csv(f'M_inverses/M_inv_hs_23_s{scaling}_FDO4', header=None))
above_T = 1.05

converge_time_length = []

for l in length_of_fins:
    # print(hs_length*6)
    
    print('\n\n----------------------------------------------\n\n')
    print(f"\n\nRunning for; {l} mm fins; scaling {scaling}")
    

    # print(trim(T_0_hs_above).shape, trim(T_0_hs_bellow).shape)
    
    M_inv_fn = np.matrix(pd.read_csv(f'M_inverses/M_inv_fn_{l}_s{scaling}_FDO4', header=None))
    # bellow_T = 2.
    
    if l == 32:
        above_T = 1.22
        bellow_T = 1.205
    elif l ==34:
        above_T = 1.21
        bellow_T = 1.195
    elif l ==36:
        above_T = 1.188
        bellow_T = 1.17
    elif l == 38:
        above_T = 1.178
        bellow_T = 1.175
    elif l == 40:
        above_T = 1.17
        bellow_T = 1.166
    elif l == 42:
        above_T = 1.163
        bellow_T = 1.158
    elif l == 44:
        above_T = 1.157
        bellow_T = 1.153
    elif l == 30:
        above_T = 1.22
        bellow_T = 1.205
    
    T_above_list = []
    above_T_1 = above_T
    diff_above = []
    T_bellow_lsit = []
    bellow_T_1 = bellow_T
    diff_bellow = []
    how_many=20
    # if l!=42:
    #     T_k_pr_above = np.matrix(pd.read_csv(f'whole_system_FDO4/fin_length_test/T_k_pr_above_s6_fin_length_{l}', header=None))
    #     T_k_cs_above = np.matrix(pd.read_csv(f'whole_system_FDO4/fin_length_test/T_k_cs_above_s6_fin_length_{l}', header=None))
    #     T_k_hs_above = np.matrix(pd.read_csv(f'whole_system_FDO4/fin_length_test/T_k_hs_above_s6_fin_length_{l}', header=None))
    #     T_k_fn_above = np.array([np.matrix(pd.read_csv(f'whole_system_FDO4/fin_length_test/T_k_fn_above_s6_fin_length_{l}_{k}', header=None)) for k in range(4)])
    #     T_0_pr_above = take_outside(T_k_pr_above)
    #     T_0_cs_above = take_outside(T_k_cs_above)
    #     T_0_hs_above = take_outside(T_k_hs_above)
    #     T_0_fn_above = np.array([take_outside(np.matrix(T_k_above)) for T_k_above in T_k_fn_above])
        
    #     T_k_pr_bellow = np.matrix(pd.read_csv(f'whole_system_FDO4/fin_length_test/T_k_pr_bellow_s6_fin_length_{l}', header=None))
    #     T_k_cs_bellow = np.matrix(pd.read_csv(f'whole_system_FDO4/fin_length_test/T_k_cs_bellow_s6_fin_length_{l}', header=None))
    #     T_k_hs_bellow = np.matrix(pd.read_csv(f'whole_system_FDO4/fin_length_test/T_k_hs_bellow_s6_fin_length_{l}', header=None))
    #     T_k_fn_bellow = np.array([np.matrix(pd.read_csv(f'whole_system_FDO4/fin_length_test/T_k_fn_bellow_s6_fin_length_{l}_{k}', header=None)) for k in range(4)])
    #     T_0_pr_bellow = take_outside(T_k_pr_bellow)
    #     T_0_cs_bellow = take_outside(T_k_cs_bellow)
    #     T_0_hs_bellow = take_outside(T_k_hs_bellow)
    #     T_0_fn_bellow = np.array([take_outside(np.matrix(T_k_bellow)) for T_k_bellow in T_k_fn_bellow])
    #     how_many = 18
    # else:
    #     how_many = 36 
    T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, vol_pr, vol_cs, vol_sinc, vol_fins = full_system_initialisation(above_T, n_fins=6, scaling=scaling, new_fin_length=l, spacing=1)
    T_0_pr_bellow, T_0_cs_bellow, T_0_hs_bellow, T_0_fn_bellow, T_k_pr_bellow, T_k_cs_bellow, T_k_hs_bellow, T_k_fn_bellow, vol_pr, vol_cs, vol_sinc, vol_fins = full_system_initialisation(bellow_T, n_fins=6, scaling=scaling, new_fin_length=l, spacing=1)
    
    T_above_list = []
    above_T_1 = above_T
    diff_above = []
    T_bellow_lsit = []
    bellow_T_1 = bellow_T
    diff_bellow = []
    t_before = time.time()
    for holding in range(how_many):
       T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, above_T = update_system_n_times(500, scaling, T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, M_inv_pr, M_inv_cs, M_inv_hs, M_inv_fn, nat, spacing=1)
       T_k_pr_bellow, T_k_cs_bellow, T_k_hs_bellow, T_k_fn_bellow, T_0_pr_bellow, T_0_cs_bellow, T_0_hs_bellow, T_0_fn_bellow, bellow_T = update_system_n_times(500, scaling, T_k_pr_bellow, T_k_cs_bellow, T_k_hs_bellow, T_k_fn_bellow, T_0_pr_bellow, T_0_cs_bellow, T_0_hs_bellow, T_0_fn_bellow, M_inv_pr, M_inv_cs, M_inv_hs, M_inv_fn, nat, spacing=1)
      
       # print(f'\n\nAbove T:   {above_T}\nBellow T:  {bellow_T}')
      
       T_above_list.append(above_T)
       diff_above.append(above_T - above_T_1)
       above_T_1 = above_T
       T_bellow_lsit.append(bellow_T)
       diff_bellow.append(bellow_T - bellow_T_1)
       bellow_T_1 = bellow_T
    
    
    t_after = time.time()
    converge_time_length.append(t_after-t_before)
    print(f'Time to converge: {int((t_after-t_before)//60)} mins {round((t_after-t_before)%60)} seconds')
    print(f'\n\nAbove T:   {above_T}\nBellow T:  {bellow_T}')
    print('Done with the third section')
    print('Saving...')
        
    path = 'whole_system_FDO4/fin_length_test_forced'
    
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
    
    
#%%

scaling = 6

q = 0.5*1e9 # W/m^3
T_a = 293 # K
K_si = 150 # W/mK
K_ceramic = 230 # W/mK
K_alu = 250 # W/mK
h = (1/scaling)*1e-3 # m      

print('Done with the first section')

nat=True
l = 30

M_inv_pr = np.matrix(pd.read_csv(f'M_inverses/M_inv_pr_s{scaling}_FDO4', header=None))
M_inv_cs = np.matrix(pd.read_csv(f'M_inverses/M_inv_cs_s{scaling}_FDO4', header=None))
M_inv_fn = np.matrix(pd.read_csv(f'M_inverses/M_inv_fn_30_s{scaling}_FDO4', header=None))
above_T = 3.9

converge_time_length = []
spacing = [1,3,4]

for s in spacing:
    # print(hs_length*6)
    
    print('\n\n----------------------------------------------\n\n')
    print(f"\n\nRunning for; {s} spacing; scaling {scaling}; 12 fins 30mm")
    

    # print(trim(T_0_hs_above).shape, trim(T_0_hs_bellow).shape)
    
    
    M_inv_hs = np.matrix(pd.read_csv(f'M_inverses/M_inv_hs_{11*(s+1)+1}_s{scaling}_FDO4', header=None))
    above_T = 3.5
    bellow_T = 2.8
    

    
    T_above_list = []
    above_T_1 = above_T
    diff_above = []
    T_bellow_lsit = []
    bellow_T_1 = bellow_T
    diff_bellow = []
    
    
    T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, vol_pr, vol_cs, vol_sinc, vol_fins = full_system_initialisation(above_T, n_fins=6, scaling=scaling, new_fin_length=30, spacing=s)
    T_0_pr_bellow, T_0_cs_bellow, T_0_hs_bellow, T_0_fn_bellow, T_k_pr_bellow, T_k_cs_bellow, T_k_hs_bellow, T_k_fn_bellow, vol_pr, vol_cs, vol_sinc, vol_fins = full_system_initialisation(bellow_T, n_fins=6, scaling=scaling, new_fin_length=30, spacing=s)
    
    t_before = time.time()
    
    for holding in range(20):
       T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, above_T = update_system_n_times(800, scaling, T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, M_inv_pr, M_inv_cs, M_inv_hs, M_inv_fn, nat, spacing=s)
       T_k_pr_bellow, T_k_cs_bellow, T_k_hs_bellow, T_k_fn_bellow, T_0_pr_bellow, T_0_cs_bellow, T_0_hs_bellow, T_0_fn_bellow, bellow_T = update_system_n_times(800, scaling, T_k_pr_bellow, T_k_cs_bellow, T_k_hs_bellow, T_k_fn_bellow, T_0_pr_bellow, T_0_cs_bellow, T_0_hs_bellow, T_0_fn_bellow, M_inv_pr, M_inv_cs, M_inv_hs, M_inv_fn, nat, spacing=s)
      
       print(f'\n\nAbove T:   {above_T}\nBellow T:  {bellow_T}')
      
       T_above_list.append(above_T)
       diff_above.append(above_T - above_T_1)
       above_T_1 = above_T
       T_bellow_lsit.append(bellow_T)
       diff_bellow.append(bellow_T - bellow_T_1)
       bellow_T_1 = bellow_T
    
    
    t_after = time.time()
    converge_time_length.append(t_after-t_before)
    print(f'Time to converge: {int((t_after-t_before)//60)} mins {round((t_after-t_before)%60)} seconds')
    print(f'\n\nAbove T:   {above_T}\nBellow T:  {bellow_T}')
    print('Done with the third section')
    print('Saving...')
        
    path = 'whole_system_FDO4/spacing_test'
    
    df = pd.DataFrame(T_k_pr_above)
    df.to_csv(f"{path}/T_k_pr_above_s{scaling}_spacing_{s}", index=False, header=False)
    df = pd.DataFrame(T_k_cs_above)
    df.to_csv(f"{path}/T_k_cs_above_s{scaling}_spacing_{s}", index=False, header=False)
    df = pd.DataFrame(T_k_hs_above)
    df.to_csv(f"{path}/T_k_hs_above_s{scaling}_spacing_{s}", index=False, header=False)
    for index, fin in enumerate(T_k_fn_above):
        df = pd.DataFrame(fin)
        df.to_csv(f"{path}/T_k_fn_above_s{scaling}_spacing_{s}", index=False, header=False)
        
    df = pd.DataFrame(T_k_pr_bellow)
    df.to_csv(f"{path}/T_k_pr_bellow_s{scaling}_spacing_{s}", index=False, header=False)
    df = pd.DataFrame(T_k_cs_bellow)
    df.to_csv(f"{path}/T_k_cs_bellow_s{scaling}_spacing_{s}", index=False, header=False)
    df = pd.DataFrame(T_k_hs_bellow)
    df.to_csv(f"{path}/T_k_hs_bellow_s{scaling}_spacing_{s}", index=False, header=False)
    for index, fin in enumerate(T_k_fn_bellow):
        df = pd.DataFrame(fin)
        df.to_csv(f"{path}/T_k_fn_bellow_s{scaling}_spacing_{s}_{index}", index=False, header=False)
    
    df = pd.DataFrame(np.matrix([T_above_list, diff_above, T_bellow_lsit, diff_bellow]).transpose())
    df.columns = ["T_above", "diff_above", "T_bellow", "diff_bellow"]
    df.to_csv(f'{path}/data_s{scaling}_spacing_{s}', index=False)

    print('\n\n            SAVED\n\n')
        
    
#%%

scaling = 6

q = 0.5*1e9 # W/m^3
T_a = 293 # K
K_si = 150 # W/mK
K_ceramic = 230 # W/mK
K_alu = 250 # W/mK
h = (1/scaling)*1e-3 # m      

print('Done with the first section')

nat=False
length_of_fins = [30,32,34,36,38,40,42,44]
length_of_fins = [30]
length_of_fins = [34,38,42,44]
M_inv_pr = np.matrix(pd.read_csv(f'M_inverses/M_inv_pr_s{scaling}_FDO4', header=None))
M_inv_cs = np.matrix(pd.read_csv(f'M_inverses/M_inv_cs_s{scaling}_FDO4', header=None))
M_inv_hs = np.matrix(pd.read_csv(f'M_inverses/M_inv_hs_23_s{scaling}_FDO4', header=None))
above_T = 1.3
bellow_T = 1.15

converge_time_length = []

for l in length_of_fins:
    # print(hs_length*6)
    
    print('\n\n----------------------------------------------\n\n')
    print(f"\n\nRunning for; {l} mm fins; scaling {scaling}")
    

    # print(trim(T_0_hs_above).shape, trim(T_0_hs_bellow).shape)
    
    M_inv_fn = np.matrix(pd.read_csv(f'M_inverses/M_inv_fn_{l}_s{scaling}_FDO4', header=None))
    # bellow_T = 2.
    
    if l == 32:
        above_T = 1.205
        bellow_T = 1.205-0.05
    elif l ==34:
        above_T = 1.192
        bellow_T = 1.185
    elif l ==36:
        above_T = 1.188
        bellow_T = 1.17
    elif l == 38:
        above_T = 1.175
        bellow_T = 1.17
    elif l == 40:
        above_T = 1.17
        bellow_T = 1.166
    elif l == 42:
        above_T = 1.16
        bellow_T = 1.153
    elif l == 44:
        above_T = 1.155
        bellow_T = 1.149
    elif l == 30:
        above_T = 1.205
        bellow_T = 1.205-0.05
    
    T_above_list = []
    above_T_1 = above_T
    diff_above = []
    T_bellow_lsit = []
    bellow_T_1 = bellow_T
    diff_bellow = []
    how_many=18
    # if l!=42:
    #     T_k_pr_above = np.matrix(pd.read_csv(f'whole_system_FDO4/fin_length_test/T_k_pr_above_s6_fin_length_{l}', header=None))
    #     T_k_cs_above = np.matrix(pd.read_csv(f'whole_system_FDO4/fin_length_test/T_k_cs_above_s6_fin_length_{l}', header=None))
    #     T_k_hs_above = np.matrix(pd.read_csv(f'whole_system_FDO4/fin_length_test/T_k_hs_above_s6_fin_length_{l}', header=None))
    #     T_k_fn_above = np.array([np.matrix(pd.read_csv(f'whole_system_FDO4/fin_length_test/T_k_fn_above_s6_fin_length_{l}_{k}', header=None)) for k in range(4)])
    #     T_0_pr_above = take_outside(T_k_pr_above)
    #     T_0_cs_above = take_outside(T_k_cs_above)
    #     T_0_hs_above = take_outside(T_k_hs_above)
    #     T_0_fn_above = np.array([take_outside(np.matrix(T_k_above)) for T_k_above in T_k_fn_above])
        
    #     T_k_pr_bellow = np.matrix(pd.read_csv(f'whole_system_FDO4/fin_length_test/T_k_pr_bellow_s6_fin_length_{l}', header=None))
    #     T_k_cs_bellow = np.matrix(pd.read_csv(f'whole_system_FDO4/fin_length_test/T_k_cs_bellow_s6_fin_length_{l}', header=None))
    #     T_k_hs_bellow = np.matrix(pd.read_csv(f'whole_system_FDO4/fin_length_test/T_k_hs_bellow_s6_fin_length_{l}', header=None))
    #     T_k_fn_bellow = np.array([np.matrix(pd.read_csv(f'whole_system_FDO4/fin_length_test/T_k_fn_bellow_s6_fin_length_{l}_{k}', header=None)) for k in range(4)])
    #     T_0_pr_bellow = take_outside(T_k_pr_bellow)
    #     T_0_cs_bellow = take_outside(T_k_cs_bellow)
    #     T_0_hs_bellow = take_outside(T_k_hs_bellow)
    #     T_0_fn_bellow = np.array([take_outside(np.matrix(T_k_bellow)) for T_k_bellow in T_k_fn_bellow])
    #     how_many = 18
    # else:
    #     how_many = 36 
    T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, vol_pr, vol_cs, vol_sinc, vol_fins = full_system_initialisation(above_T, n_fins=6, scaling=scaling, new_fin_length=l, spacing=1)
    T_0_pr_bellow, T_0_cs_bellow, T_0_hs_bellow, T_0_fn_bellow, T_k_pr_bellow, T_k_cs_bellow, T_k_hs_bellow, T_k_fn_bellow, vol_pr, vol_cs, vol_sinc, vol_fins = full_system_initialisation(bellow_T, n_fins=6, scaling=scaling, new_fin_length=l, spacing=1)
    
    T_above_list = []
    T_pr_above_list = []
    above_T_1 = above_T
    diff_above = []
    T_bellow_lsit = []
    T_pr_below_list = []
    bellow_T_1 = bellow_T
    diff_bellow = []
    t_before = time.time()
    for holding in range(how_many):
       T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, above_T = update_system_n_times(500, scaling, T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, M_inv_pr, M_inv_cs, M_inv_hs, M_inv_fn, nat, spacing=1)
       T_k_pr_bellow, T_k_cs_bellow, T_k_hs_bellow, T_k_fn_bellow, T_0_pr_bellow, T_0_cs_bellow, T_0_hs_bellow, T_0_fn_bellow, bellow_T = update_system_n_times(500, scaling, T_k_pr_bellow, T_k_cs_bellow, T_k_hs_bellow, T_k_fn_bellow, T_0_pr_bellow, T_0_cs_bellow, T_0_hs_bellow, T_0_fn_bellow, M_inv_pr, M_inv_cs, M_inv_hs, M_inv_fn, nat, spacing=1)
      
       # print(f'\n\nAbove T:   {above_T}\nBellow T:  {bellow_T}')
       T_pr_above = average_temp(T_k_pr_above)
       T_pr_below = average_temp(T_k_pr_bellow)
                
       T_above_list.append(above_T)
       T_pr_above_list.append(T_pr_above)
       diff_above.append(above_T - above_T_1)
       above_T_1 = above_T
       T_bellow_lsit.append(bellow_T)
       T_pr_below_list.append(T_pr_below)
       diff_bellow.append(bellow_T - bellow_T_1)
       bellow_T_1 = bellow_T
    
    
    t_after = time.time()
    converge_time_length.append(t_after-t_before)
    print(f'Time to converge: {int((t_after-t_before)//60)} mins {round((t_after-t_before)%60)} seconds')
    print(f'\n\nAbove T:   {above_T}\nBellow T:  {bellow_T}')
    print('Done with the third section')
    print('Saving...')
        
    path = 'whole_system_FDO4/fin_length_test_forced'
    
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
    
    df = pd.DataFrame(np.matrix([T_above_list, diff_above, T_bellow_lsit, diff_bellow, T_pr_above_list, T_pr_below_list]).transpose())
    df.columns = ["T_above", "diff_above", "T_bellow", "diff_bellow", 'T_pr_above', 'T_pr_below']
    df.to_csv(f'{path}/data_s{scaling}_fin_length_{l}', index=False)

    print('\n\n            SAVED\n\n')
      
#%%
from MatrixFunctions import matrix_maker_2, gaus_jordan

scaling = 6

q = 0.5*1e9 # W/m^3
T_a = 293 # K
K_si = 150 # W/mK
K_ceramic = 230 # W/mK
K_alu = 250 # W/mK
h = (1/scaling)*1e-3 # m      

print('Done with the first section')

nat=False
length_of_fins = [30,32,34,36,38,40,42,44]
length_of_fins = [58]
# length_of_fins = [30,34,38,42,44]
length_of_fins = [26,28]
length_of_fins = [39]
n_fins_list = [20]
M_inv_pr = np.matrix(pd.read_csv(f'M_inverses/M_inv_pr_s{scaling}_FDO4', header=None))
M_inv_cs = np.matrix(pd.read_csv(f'M_inverses/M_inv_cs_s{scaling}_FDO4', header=None))
# M_inv_hs = np.matrix(pd.read_csv(f'M_inverses/M_inv_hs_23_s{scaling}_FDO4', header=None))
above_T = 1.14
bellow_T = 1.05

converge_time_length = []

for l, n_fins in zip(length_of_fins, n_fins_list):
    # print(hs_length*6)
    
    print('\n\n----------------------------------------------\n\n')
    print(f"\n\nRunning for; {l} mm fins; scaling {scaling}")
    

    # print(trim(T_0_hs_above).shape, trim(T_0_hs_bellow).shape)
    
    # M_inv_fn = np.matrix(pd.read_csv(f'M_inverses/M_inv_fn_s{scaling}_FDO4', header=None))
    
    # bellow_T = 2.
    
    # if l == 32:
    #     above_T = 1.205
    #     bellow_T = 1.205-0.05
    # elif l ==34:
    #     above_T = 1.192
    #     bellow_T = 1.85
    # elif l ==36:
    #     above_T = 1.188
    #     bellow_T = 1.17
    # elif l == 38:
    #     above_T = 1.175
    #     bellow_T = 1.17
    # elif l == 40:
    #     above_T = 1.17
    #     bellow_T = 1.166
    # elif l == 42:
    #     above_T = 1.16
    #     bellow_T = 1.153
    # elif l == 44:
    #     above_T = 1.155
    #     bellow_T = 1.149
    # elif l == 30:
    #     above_T = 1.205
    #     bellow_T = 1.205-0.05
    above_T = 1.16
    bellow_T = 1.05
    
    T_above_list = []
    above_T_1 = above_T
    diff_above = []
    T_bellow_lsit = []
    bellow_T_1 = bellow_T
    diff_bellow = []
    how_many=15
    # if l!=42:
    #     T_k_pr_above = np.matrix(pd.read_csv(f'whole_system_FDO4/fin_length_test/T_k_pr_above_s6_fin_length_{l}', header=None))
    #     T_k_cs_above = np.matrix(pd.read_csv(f'whole_system_FDO4/fin_length_test/T_k_cs_above_s6_fin_length_{l}', header=None))
    #     T_k_hs_above = np.matrix(pd.read_csv(f'whole_system_FDO4/fin_length_test/T_k_hs_above_s6_fin_length_{l}', header=None))
    #     T_k_fn_above = np.array([np.matrix(pd.read_csv(f'whole_system_FDO4/fin_length_test/T_k_fn_above_s6_fin_length_{l}_{k}', header=None)) for k in range(4)])
    #     T_0_pr_above = take_outside(T_k_pr_above)
    #     T_0_cs_above = take_outside(T_k_cs_above)
    #     T_0_hs_above = take_outside(T_k_hs_above)
    #     T_0_fn_above = np.array([take_outside(np.matrix(T_k_above)) for T_k_above in T_k_fn_above])
        
    #     T_k_pr_bellow = np.matrix(pd.read_csv(f'whole_system_FDO4/fin_length_test/T_k_pr_bellow_s6_fin_length_{l}', header=None))
    #     T_k_cs_bellow = np.matrix(pd.read_csv(f'whole_system_FDO4/fin_length_test/T_k_cs_bellow_s6_fin_length_{l}', header=None))
    #     T_k_hs_bellow = np.matrix(pd.read_csv(f'whole_system_FDO4/fin_length_test/T_k_hs_bellow_s6_fin_length_{l}', header=None))
    #     T_k_fn_bellow = np.array([np.matrix(pd.read_csv(f'whole_system_FDO4/fin_length_test/T_k_fn_bellow_s6_fin_length_{l}_{k}', header=None)) for k in range(4)])
    #     T_0_pr_bellow = take_outside(T_k_pr_bellow)
    #     T_0_cs_bellow = take_outside(T_k_cs_bellow)
    #     T_0_hs_bellow = take_outside(T_k_hs_bellow)
    #     T_0_fn_bellow = np.array([take_outside(np.matrix(T_k_bellow)) for T_k_bellow in T_k_fn_bellow])
    #     how_many = 18
    # else:
    #     how_many = 36 
    T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, vol_pr, vol_cs, vol_sinc, vol_fins = full_system_initialisation(above_T, n_fins=n_fins/2, scaling=scaling, new_fin_length=l, spacing=1)
    T_0_pr_bellow, T_0_cs_bellow, T_0_hs_bellow, T_0_fn_bellow, T_k_pr_bellow, T_k_cs_bellow, T_k_hs_bellow, T_k_fn_bellow, vol_pr, vol_cs, vol_sinc, vol_fins = full_system_initialisation(bellow_T, n_fins=n_fins/2, scaling=scaling, new_fin_length=l, spacing=1)
    

    M = matrix_maker_2(T_0_fn_above[0])
    M_inv_fn = gaus_jordan(M)
    M = matrix_maker_2(T_0_hs_above)
    M_inv_hs = gaus_jordan(M)
    
    T_above_list = []
    T_pr_above_list = []
    above_T_1 = above_T
    diff_above = []
    T_bellow_lsit = []
    T_pr_below_list = []
    bellow_T_1 = bellow_T
    diff_bellow = []
    t_before = time.time()
    for holding in range(how_many):
       T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, above_T = update_system_n_times(500, scaling, T_k_pr_above, T_k_cs_above, T_k_hs_above, T_k_fn_above, T_0_pr_above, T_0_cs_above, T_0_hs_above, T_0_fn_above, M_inv_pr, M_inv_cs, M_inv_hs, M_inv_fn, nat, spacing=1)
       T_k_pr_bellow, T_k_cs_bellow, T_k_hs_bellow, T_k_fn_bellow, T_0_pr_bellow, T_0_cs_bellow, T_0_hs_bellow, T_0_fn_bellow, bellow_T = update_system_n_times(500, scaling, T_k_pr_bellow, T_k_cs_bellow, T_k_hs_bellow, T_k_fn_bellow, T_0_pr_bellow, T_0_cs_bellow, T_0_hs_bellow, T_0_fn_bellow, M_inv_pr, M_inv_cs, M_inv_hs, M_inv_fn, nat, spacing=1)
      
       print(f'\n\nAbove T:   {above_T}\nBellow T:  {bellow_T}')
       T_pr_above = average_temp(T_k_pr_above)
       T_pr_below = average_temp(T_k_pr_bellow)
                
       T_above_list.append(above_T)
       T_pr_above_list.append(T_pr_above)
       diff_above.append(above_T - above_T_1)
       above_T_1 = above_T
       T_bellow_lsit.append(bellow_T)
       T_pr_below_list.append(T_pr_below)
       diff_bellow.append(bellow_T - bellow_T_1)
       bellow_T_1 = bellow_T
    
    
    t_after = time.time()
    converge_time_length.append(t_after-t_before)
    print(f'Time to converge: {int((t_after-t_before)//60)} mins {round((t_after-t_before)%60)} seconds')
    print(f'\n\nAbove T:   {above_T}\nBellow T:  {bellow_T}')
    print('Done with the third section')
    print('Saving...')
        
    path = 'whole_system_FDO4/fin_length_test_forced'
    
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
    
    df = pd.DataFrame(np.matrix([T_above_list, diff_above, T_bellow_lsit, diff_bellow, T_pr_above_list, T_pr_below_list]).transpose())
    df.columns = ["T_above", "diff_above", "T_bellow", "diff_bellow", 'T_pr_above', 'T_pr_below']
    df.to_csv(f'{path}/data_s{scaling}_fin_length_{l}', index=False)

    print('\n\n            SAVED\n\n')  
    