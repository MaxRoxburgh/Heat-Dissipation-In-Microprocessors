from Updates import calc_update_shared_bdd, heat_sinc_boundary, initialise_boundary
from Updates import update_with_source, update_without_source, np, update_boundary_convection
from Calculators import flip_fins, mirror_array


q = 0.5*1e9 # W/m^3
T_a = 293 # K
K_si = 150 # W/mK
K_ceramic = 230 # W/mK
K_alu = 250 # W/mK
h = (1/3)*1e-3 # m   

n_fins = 6     

import pandas as pd
import numpy as np
from Plotting import take_outside, plot, plot_whole_system

#%%
#T_k_pr = pd.read_csv('Optimised_T_k/T_k_pr_udpate_small_low_unexpec', header=None)
#%%

change = -2.2

T_k_pr = np.matrix(pd.read_csv('Optimised_T_k/T_k_pr_udpate_small_low_unexpec', header=None)) + change
T_k_cs = np.matrix(pd.read_csv('Optimised_T_k/T_k_cs_udpate_small_low_unexpec', header=None)) + change
T_k_hs = np.matrix(pd.read_csv('Optimised_T_k/T_k_hs_udpate_small_low_unexpec', header=None)) + change
T_k_fn = np.array([np.matrix(pd.read_csv(f'Optimised_T_k/T_k_fn_{i}_small_low_unexpec', header=None)) + change for i in range(n_fins)])


T_0_pr = take_outside(T_k_pr)
T_0_cs = take_outside(T_k_cs)
T_0_hs = take_outside(T_k_hs)
T_0_fn = np.array([take_outside(np.matrix(T_k)) for T_k in T_k_fn])
ambient = min(T_k_fn[0][0])
plot_whole_system(T_k_cs, T_k_pr, mirror_array(T_k_fn, True), T_k_hs, ambient, 6)

average_T_pr = average_temp(T_k_pr)
average_T_cs = average_temp(T_k_cs)
average_T_fn = sum([average_temp(np.matrix(T_k)) for T_k in T_k_fn])/len(T_k_fn)
average_T_hs = average_temp(T_k_cs)
avT_k_1 = (2*vol_fins*average_T_fn + vol_sinc*average_T_hs + vol_cs*average_T_cs + vol_pr* average_T_pr)/(2*vol_fins+vol_sinc+vol_cs+vol_pr)
#%% We will find the equelibrium with constant BC then update them to see what direction it tends to
###


from Calculators import average_temp

for i in range(20000):
    
    # update boundaries convection
    #T_0_pr_update, T_0_cs_update = update_boundary_convection(T_k_pr_update, K_si), update_boundary_convection(T_k_cs_update, K_ceramic)
    #T_0_hs_update = update_boundary_convection(T_k_hs_update, K_alu)
    #T_0_fn_update = np.array([update_boundary_convection(np.matrix(fin_k), K_alu) for fin_k in T_k_fn_update])

    
    # update boudaries shared
    #T_0_pr_update, T_0_cs_update = calc_update_shared_bdd(T_k_pr_update[2].A1, T_k_cs_update[-3].A1, K_si, K_ceramic, T_0_pr_update, T_0_cs_update)
    #T_0_cs_update, T_0_hs_update = calc_update_shared_bdd(T_k_cs_update[2].A1, T_k_hs_update[-3].A1, K_ceramic, K_alu, T_0_cs_update, T_0_hs_update)
    #T_above_list = np.array([fin_k[-3] for fin_k in mirror_array(T_k_fn_update, True)])
    #T_below = T_k_hs_update[-3].A1
    #T_0_fn_update_2 = mirror_array(T_0_fn_update, True)
    #T_0_hs_update, T_0_fn_update = heat_sinc_boundary(T_below, T_above_list, T_0_hs_update, T_0_fn_update_2, 6, 3)
    #T_0_fn_update = np.array([T_0_fn_update[i] for i in range(n_fins)])
    
    
    # update inside temps
    T_k_pr = update_with_source(T_k_pr, T_0_pr, q, h)
    T_k_cs = update_without_source(T_k_cs, T_0_cs)
    T_k_hs = update_without_source(T_k_hs, T_0_hs)
    T_k_fn = np.array([update_without_source(np.matrix(fin_k), fin_0) for fin_0, fin_k in zip(T_0_fn, T_k_fn)])

    
    if i%1000 == 0:
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
        avT_k_1 = avT_k

    if i%500 == 0:

        diff_array_big.append(diff)
        av_array_big.append(avT_k)
        av_pr_array.append(average_T_pr)
        av_cs_array.append(average_T_cs)
        av_hs_array.append(average_fin_hs)


#%%


for i in range(20000):
    
    # update boundaries convection
    T_0_pr, T_0_cs = update_boundary_convection(T_k_pr, K_si), update_boundary_convection(T_k_cs, K_ceramic)
    T_0_hs = update_boundary_convection(T_k_hs, K_alu)
    T_0_fn = np.array([update_boundary_convection(np.matrix(fin_k), K_alu) for fin_k in T_k_fn])

    
    # update boudaries shared
    T_0_pr, T_0_cs = calc_update_shared_bdd(T_k_pr[2].A1, T_k_cs[-3].A1, K_si, K_ceramic, T_0_pr, T_0_cs)
    T_0_cs, T_0_hs = calc_update_shared_bdd(T_k_cs[2].A1, T_k_hs[-3].A1, K_ceramic, K_alu, T_0_cs, T_0_hs)
    T_above_list = np.array([fin_k[-3] for fin_k in mirror_array(T_k_fn, True)])
    T_below = T_k_hs[-3].A1
    T_0_fn_2 = mirror_array(T_0_fn, True)
    T_0_hs, T_0_fn = heat_sinc_boundary(T_below, T_above_list, T_0_hs, T_0_fn_2, 6, 3)
    T_0_fn = np.array([T_0_fn[i] for i in range(n_fins)])
    
    
    # update inside temps
    T_k_pr = update_with_source(T_k_pr, T_0_pr, q, h)
    T_k_cs = update_without_source(T_k_cs, T_0_cs)
    T_k_hs = update_without_source(T_k_hs, T_0_hs)
    T_k_fn = np.array([update_without_source(np.matrix(fin_k), fin_0) for fin_0, fin_k in zip(T_0_fn, T_k_fn)])

    
    if i%1000 == 0:
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
        avT_k_1 = avT_k

    if i%500 == 0:

        diff_array_big.append(diff)
        av_array_big.append(avT_k)
        av_pr_array.append(average_T_pr)
        av_cs_array.append(average_T_cs)
        av_hs_array.append(average_fin_hs)






