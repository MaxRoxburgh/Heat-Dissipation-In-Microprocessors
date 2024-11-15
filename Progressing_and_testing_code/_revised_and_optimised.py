from Updates import calc_update_shared_bdd, heat_sinc_boundary, initialise_boundary
from Updates import update_with_source, update_without_source, np, update_boundary_convection
from Calculators import flip_fins, mirror_array
from Plotting import temp_calc


q = 0.5*1e9 # W/m^3
T_a = 293 # K
K_si = 150 # W/mK
K_ceramic = 230 # W/mK
K_alu = 250 # W/mK
h = (1/3)*1e-3 # m         

#%% initialising the system

ambient = 4

# set up the processor

vol_pr = 14*3*3

top = np.array([ambient for i in range(14*3)])
bottom = top.copy()
left = np.array([ambient for i in range(3)])
right = left.copy()

T_0_pr, T_k_pr = initialise_boundary(top, left, right, bottom, ambient,[ambient,ambient,ambient,ambient])

# set up the case

vol_cs = 60*6

top = np.array([ambient for i in range(60)])
bottom = top.copy()
left = np.array([ambient for i in range(6)])
right = left.copy()

T_0_cs, T_k_cs = initialise_boundary(top, left, right, bottom, ambient,[ambient,ambient,ambient,ambient])


# set up the main body of the heat sinc (34x4)

sinc_length = 34*3
spacing = 2*3
vol_sinc = sinc_length*4*3

top = np.array([ambient for i in range(sinc_length)])
bottom = top.copy()
left = np.array([ambient for i in range(12)])
right = left.copy()

T_0_hs, T_k_hs = initialise_boundary(top, left, right, bottom, ambient,[ambient,ambient,ambient,ambient])

# set up fins for heat sinc

fin_length = 90

top = np.array([ambient for i in range(3)])
bottom = top.copy()
left = np.array([ambient for i in range(fin_length)])
right = left.copy()

n_fins = int(((sinc_length+spacing)/(3+spacing))/2)

vol_fins = n_fins*3*fin_length
fin_0, fin_k = initialise_boundary(top, left, right, bottom, ambient,[ambient,ambient,ambient,ambient])
T_0_fn = np.array([fin_0 for i in range(n_fins)])
T_k_fn = np.array([fin_k for i in range(n_fins)])


avT_k = ambient
#%% let it smooth out to begin with assuming neuman BC, it stops when the change is 0
for i in range(100000):
    T_0_pr, T_0_cs = calc_update_shared_bdd(T_k_pr[2].A1, T_k_cs[-3].A1, K_si, K_ceramic, T_0_pr, T_0_cs)
    T_0_cs, T_0_hs = calc_update_shared_bdd(T_k_cs[2].A1, T_k_hs[-3].A1, K_ceramic, K_alu, T_0_cs, T_0_hs)
    T_above_list = np.array([fin_k[-3] for fin_k in mirror_array(T_k_fn, True)])
    T_below = T_k_hs[-3].A1
    T_0_fn_2 = mirror_array(T_0_fn, True)
    T_0_hs, T_0_fn = heat_sinc_boundary(T_below, T_above_list, T_0_hs, T_0_fn_2, 6, 3)
    T_0_fn = np.array([T_0_fn[i] for i in range(n_fins)])
    T_k_pr = update_with_source(T_k_pr, T_0_pr, q, h)
    T_k_cs = update_without_source(T_k_cs, T_0_cs)
    T_k_hs = update_without_source(T_k_hs, T_0_hs)
    T_k_fn = np.array([update_without_source(np.matrix(fin_k), fin_0) for fin_0, fin_k in zip(T_0_fn, T_k_fn)])
    if i%2000 == 0:
        avT_k,diff = temp_calc(T_k_pr,T_k_cs,T_k_hs,T_k_fn, avT_k, vol_fins, vol_cs, vol_pr, vol_sinc, i)
        
    if diff == 0.:
        break
ambient = 1
#%%
# import matplotlib.pyplot as plt
# plt.imshow(T_k_cs, vmin=)
# plt.show()
# plt.imshow(T_0_cs, vmin=9)
# plt.show()

#%% plot to see initial system

import Plotting as plot

T_k_fn_whole = mirror_array(T_k_fn, True)
plot.plot_whole_system(T_k_cs, T_k_pr, T_k_fn_whole, T_k_hs, 8.999, 6)

#%% fist boundary update: nat convection

T_0_pr_update, T_0_cs_update = update_boundary_convection(T_k_pr, K_si), update_boundary_convection(T_k_cs, K_ceramic)
T_0_hs_update = update_boundary_convection(T_k_hs, K_alu)
T_0_fn_update = np.array([update_boundary_convection(np.matrix(fin_k), K_alu) for fin_k in T_k_fn])
    
# plt.imshow(T_0_cs_update)
# print(min(T_0_cs_update.A1))
#%% fist boundary update: shared boundies

T_0_pr_update, T_0_cs_update = calc_update_shared_bdd(T_k_pr[2].A1, T_k_cs[-3].A1, K_si, K_ceramic, T_0_pr_update, T_0_cs_update)
T_0_cs_update, T_0_hs_update = calc_update_shared_bdd(T_k_cs[2].A1, T_k_hs[-3].A1, K_ceramic, K_alu, T_0_cs_update, T_0_hs_update)

T_above_list = np.array([fin_k[-3] for fin_k in mirror_array(T_k_fn, True)])
T_below = T_k_hs[-3].A1
T_0_fn_update_2 = mirror_array(T_0_fn_update, True)
T_0_hs_update, T_0_fn_update = heat_sinc_boundary(T_below, T_above_list, T_0_hs_update, T_0_fn_update_2, 6, 3)
T_0_fn_update = np.array([T_0_fn_update[i] for i in range(n_fins)])

#%% fist inside update

T_k_pr_update = update_with_source(T_k_pr, T_0_pr_update, q, h)
T_k_cs_update = update_without_source(T_k_cs, T_0_cs_update)
T_k_hs_update = update_without_source(T_k_hs, T_0_hs_update)
T_k_fn_update = np.array([update_without_source(np.matrix(fin_k), fin_0) for fin_0, fin_k in zip(T_0_fn_update, T_k_fn)])

#%% iterations!!!

from Calculators import average_temp

# avT_k_1 = ambient
# diff_array = []
# av_array = []
#avT_k_1 = ambient
diff_array_big = []
av_array_big = []
av_pr_array = []
av_cs_array = []
av_hs_array = []


for i in range(250000):
    
    # update boundaries convection
    T_0_pr_update, T_0_cs_update = update_boundary_convection(T_k_pr_update, K_si), update_boundary_convection(T_k_cs_update, K_ceramic)
    T_0_hs_update = update_boundary_convection(T_k_hs_update, K_alu)
    T_0_fn_update = np.array([update_boundary_convection(np.matrix(fin_k), K_alu) for fin_k in T_k_fn_update])

    
    # update boudaries shared
    T_0_pr_update, T_0_cs_update = calc_update_shared_bdd(T_k_pr_update[2].A1, T_k_cs_update[-3].A1, K_si, K_ceramic, T_0_pr_update, T_0_cs_update)
    T_0_cs_update, T_0_hs_update = calc_update_shared_bdd(T_k_cs_update[2].A1, T_k_hs_update[-3].A1, K_ceramic, K_alu, T_0_cs_update, T_0_hs_update)
    T_above_list = np.array([fin_k[-3] for fin_k in mirror_array(T_k_fn_update, True)])
    T_below = T_k_hs_update[-3].A1
    T_0_fn_update_2 = mirror_array(T_0_fn_update, True)
    T_0_hs_update, T_0_fn_update = heat_sinc_boundary(T_below, T_above_list, T_0_hs_update, T_0_fn_update_2, 6, 3)
    T_0_fn_update = np.array([T_0_fn_update[i] for i in range(n_fins)])
    
    
    # update inside temps
    T_k_pr_update = update_with_source(T_k_pr_update, T_0_pr_update, q, h)
    T_k_cs_update = update_without_source(T_k_cs_update, T_0_cs_update)
    T_k_hs_update = update_without_source(T_k_hs_update, T_0_hs_update)
    T_k_fn_update = np.array([update_without_source(np.matrix(fin_k), fin_0) for fin_0, fin_k in zip(T_0_fn_update, T_k_fn_update)])

    
    if i%2000 == 0:
        average_T_pr = average_temp(T_k_pr_update)
        average_T_cs = average_temp(T_k_cs_update)
        average_T_fn = sum([average_temp(np.matrix(T_k)) for T_k in T_k_fn_update])/len(T_k_fn_update)
        average_T_hs = average_temp(T_k_cs_update)
        
        average_fin_hs = (2*vol_fins*average_T_fn + vol_sinc*average_T_hs)/(2*vol_fins+vol_sinc)
        
        print('T:     pr     |     cs     |      hs      |\n', average_T_pr, average_T_cs, average_fin_hs)
        ambient = min(T_k_fn_update[0][0])
        plot.plot_whole_system(T_k_cs_update, T_k_pr_update, mirror_array(T_k_fn_update, True), T_k_hs_update, ambient, 6)    
        
        avT_k = (2*vol_fins*average_T_fn + vol_sinc*average_T_hs + vol_cs*average_T_cs + vol_pr* average_T_pr)/(2*vol_fins+vol_sinc+vol_cs+vol_pr)
        diff = avT_k-avT_k_1
        print(f'\nAverage T, difference from 2000 runs ago: {i}\n',avT_k, diff)
        avT_k_1 = avT_k

    if i%2000 == 0:

        diff_array_big.append(diff)
        av_array_big.append(avT_k)
        av_pr_array.append(average_T_pr)
        av_cs_array.append(average_T_cs)
        av_hs_array.append(average_fin_hs)
        
#%%
import pandas as pd
import numpy as np

df = pd.DataFrame(np.matrix([diff_array, av_array]).transpose())
df.columns = ["diff_array", "av_array"]
df.to_csv('Optimised_T_k/Corrected_func_med', index=False)

#%%
from Plotting import save

title = 'corrected_low_2'

df = pd.DataFrame(np.matrix([diff_array_big, av_array_big, av_pr_array, av_cs_array, av_hs_array]).transpose())
df.columns = ["diff_array_big", "av_array_big", "av_pr_array", "av_cs_array", "av_hs_array"]
df.to_csv('Optimised_T_k/Corrected_func_high_2', index=False)


save(T_k_pr_update, T_k_cs_update, T_k_hs_update, T_k_fn_update, title)

        

#%%

import pandas as pd
import numpy as np
from Plotting import take_outside

T_k_pr_update = np.matrix(pd.read_csv('Optimised_T_k/T_k_pr_udpate_small_above_2'))
T_k_cs_update = np.matrix(pd.read_csv('Optimised_T_k/T_k_cs_udpate_small_above_2'))
T_k_hs_update = np.matrix(pd.read_csv('Optimised_T_k/T_k_hs_udpate_small_above_2'))
T_k_fn_update = np.array([np.matrix(pd.read_csv(f'Optimised_T_k/T_k_fn_{i}_small_above_2')) for i in range(n_fins)])

T_0_pr_update = take_outside(T_k_pr_update)
T_0_cs_update = take_outside(T_k_cs_update)
T_0_hs_update = take_outside(T_k_hs_update)
T_0_fn_update = np.array([take_outside(np.matrix(T_k)) for T_k in T_k_fn_update])


#%%

from Calculators import average_temp

# avT_k_1 = ambient
# diff_array_big = []
# av_array_big = []
# av_pr_array = []
# av_cs_array = []
# av_hs_array = []

# every = 5000
# change = 0


# T_k_pr_update = update_with_source(T_k_pr_update+change, take_outside(T_k_pr_update+change), q, h)
# T_k_cs_update = update_without_source(T_k_cs_update+change, take_outside(T_k_cs_update+change))
# T_k_hs_update = update_without_source(T_k_hs_update+change, take_outside(T_k_hs_update+change))
# T_k_fn_update = np.array([update_without_source(np.matrix(fin_k)+change, take_outside(np.matrix(fin_k)+change)) for fin_k in T_k_fn_update])

for i in range(500000):
    
    # update boundaries convection
    T_0_pr_update, T_0_cs_update = update_boundary_convection(T_k_pr_update, K_si), update_boundary_convection(T_k_cs_update, K_ceramic)
    T_0_hs_update = update_boundary_convection(T_k_hs_update, K_alu)
    T_0_fn_update = np.array([update_boundary_convection(np.matrix(fin_k), K_alu) for fin_k in T_k_fn_update])

    
    # update boudaries shared
    T_0_pr_update, T_0_cs_update = calc_update_shared_bdd(T_k_pr_update[2].A1, T_k_cs_update[-3].A1, K_si, K_ceramic, T_0_pr_update, T_0_cs_update)
    T_0_cs_update, T_0_hs_update = calc_update_shared_bdd(T_k_cs_update[2].A1, T_k_hs_update[-3].A1, K_ceramic, K_alu, T_0_cs_update, T_0_hs_update)
    T_above_list = np.array([fin_k[-3] for fin_k in mirror_array(T_k_fn_update, True)])
    T_below = T_k_hs_update[-3].A1
    T_0_fn_update_2 = mirror_array(T_0_fn_update, True)
    T_0_hs_update, T_0_fn_update = heat_sinc_boundary(T_below, T_above_list, T_0_hs_update, T_0_fn_update_2, 6, 3)
    T_0_fn_update = np.array([T_0_fn_update[i] for i in range(n_fins)])
    
    
    # update inside temps
    T_k_pr_update = update_with_source(T_k_pr_update, T_0_pr_update, q, h)
    T_k_cs_update = update_without_source(T_k_cs_update, T_0_cs_update)
    T_k_hs_update = update_without_source(T_k_hs_update, T_0_hs_update)
    T_k_fn_update = np.array([update_without_source(np.matrix(fin_k), fin_0) for fin_0, fin_k in zip(T_0_fn_update, T_k_fn_update)])

    
    if i%every == 0:
        average_T_pr = average_temp(T_k_pr_update)
        average_T_cs = average_temp(T_k_cs_update)
        average_T_fn = sum([average_temp(np.matrix(T_k)) for T_k in T_k_fn_update])/len(T_k_fn_update)
        average_T_hs = average_temp(T_k_cs_update)
        
        average_fin_hs = (2*vol_fins*average_T_fn + vol_sinc*average_T_hs)/(2*vol_fins+vol_sinc)
        
        print('T:     pr     |     cs     |      hs      |\n', average_T_pr, average_T_cs, average_fin_hs)
        
        ambient = min(T_k_fn_update[0][0])
        plot.plot_whole_system(T_k_cs_update, T_k_pr_update, mirror_array(T_k_fn_update, True), T_k_hs_update, ambient, 6)
        
        avT_k = (2*vol_fins*average_T_fn + vol_sinc*average_T_hs + vol_cs*average_T_cs + vol_pr* average_T_pr)/(2*vol_fins+vol_sinc+vol_cs+vol_pr)
        diff = avT_k-avT_k_1
        print(f'\nAverage T, difference from {every} runs ago: run {i}\n',avT_k, diff)
        avT_k_1 = avT_k

    if i%every == 0 and i >= 17000:

        diff_array_big.append(diff)
        av_array_big.append(avT_k)
        av_pr_array.append(average_T_pr)
        av_cs_array.append(average_T_cs)
        av_hs_array.append(average_fin_hs)


#%%

df = pd.DataFrame(T_k_pr_update)
df.to_csv("Optimised_T_k/T_k_pr_udpate_small_above_2_contiunued", index=False, header=False)
df = pd.DataFrame(T_k_cs_update)
df.to_csv("Optimised_T_k/T_k_cs_udpate_small_above_2_contiunued", index=False, header=False)
df = pd.DataFrame(T_k_hs_update)
df.to_csv("Optimised_T_k/T_k_hs_udpate_small_above_2_contiunued", index=False, header=False)
for i, T_k in enumerate(T_k_fn_update):
    df = pd.DataFrame(T_k)
    df.to_csv(f"Optimised_T_k/T_k_fn_{i}_small_above_2_contiunued", index=False, header=False)

#%%

df = pd.DataFrame(np.matrix([diff_array_big, av_array_big, av_pr_array, av_cs_array, av_hs_array]).transpose())
df.columns = ["diff_array_big", "av_array_big", "av_pr_array", "av_cs_array", "av_hs_array"]
df.to_csv('Optimised_T_k/Temp_arrays_small_above_2_contiunued', index=False)

### Note for the below, it calculated the total av temp and the av temp of the total hs wrong, this has been finxed for the "above" data










