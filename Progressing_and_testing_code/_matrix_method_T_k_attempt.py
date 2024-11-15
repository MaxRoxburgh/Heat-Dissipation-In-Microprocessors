from Updates import calc_update_shared_bdd, heat_sinc_boundary, initialise_boundary
from Updates import update_with_source, update_without_source, np, update_boundary_convection
from Calculators import flip_fins, mirror_array
from Plotting import temp_calc
import Plotting as plot
from MatrixFunctions import update_T_matrix_method    
from Calculators import average_temp
import pandas as pd
import matplotlib.pyplot as plt


q = 0.5*1e9 # W/m^3
T_a = 293 # K
K_si = 150 # W/mK
K_ceramic = 230 # W/mK
K_alu = 250 # W/mK
h = (1/3)*1e-3 # m         

#%%
n_fins = 6

M_inv_pr = np.matrix(pd.read_csv('inverses/M_inv_pr', header=None))
M_inv_cs = np.matrix(pd.read_csv('inverses/M_inv_cs', header=None))
M_inv_hs = np.matrix(pd.read_csv('inverses/M_inv_hs', header=None))
M_inv_fn = np.matrix(pd.read_csv('inverses/M_inv_fn', header=None))

#%% initialising the system

ambient = 3.9

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

#%% fist boundary update: nat convection

T_0_pr_update, T_0_cs_update = update_boundary_convection(T_k_pr, K_si), update_boundary_convection(T_k_cs, K_ceramic)
T_0_hs_update = update_boundary_convection(T_k_hs, K_alu)
T_0_fn_update = np.array([update_boundary_convection(np.matrix(fin_k), K_alu) for fin_k in T_k_fn])
    

#%% fist boundary update: shared boundies

T_0_pr_update, T_0_cs_update = calc_update_shared_bdd(T_k_pr[2].A1, T_k_cs[-3].A1, K_si, K_ceramic, T_0_pr_update, T_0_cs_update)
T_0_cs_update, T_0_hs_update = calc_update_shared_bdd(T_k_cs[2].A1, T_k_hs[-3].A1, K_ceramic, K_alu, T_0_cs_update, T_0_hs_update)

T_above_list = np.array([fin_k[-3] for fin_k in mirror_array(T_k_fn, True)])
T_below = T_k_hs[-3].A1
T_0_fn_update_2 = mirror_array(T_0_fn_update, True)
T_0_hs_update, T_0_fn_update = heat_sinc_boundary(T_below, T_above_list, T_0_hs_update, T_0_fn_update_2, 6, 3)
T_0_fn_update = np.array([T_0_fn_update[i] for i in range(n_fins)])
#%% fist T_k update

T_k_pr_update = update_T_matrix_method(T_0_pr_update, M_inv_pr, q, h)
T_k_cs_update = update_T_matrix_method(T_0_cs_update, M_inv_cs)
T_k_hs_update = update_T_matrix_method(T_0_hs_update, M_inv_hs)
T_k_fn_update = np.array([update_T_matrix_method(np.matrix(fin_0), M_inv_fn) for fin_0 in T_0_fn_update])
    

#%% iterations!!!
avT_k_1 = ambient
diff_array = []
av_array = []
avT_k_1 = ambient
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
    T_k_pr_update = update_T_matrix_method(T_0_pr_update, M_inv_pr, q, h)
    T_k_cs_update = update_T_matrix_method(T_0_cs_update, M_inv_cs)
    T_k_hs_update = update_T_matrix_method(T_0_hs_update, M_inv_hs)
    T_k_fn_update = np.array([update_T_matrix_method(np.matrix(fin_0), M_inv_fn) for fin_0 in T_0_fn_update])
    
    if i%2500 == 0:
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
        print(f'\nAverage T, difference from 2500 runs ago: {i}\n',avT_k, diff)
        avT_k_1 = avT_k
            
        
    if i%2500 == 0 and i>0:

        diff_array_big.append(diff)
        av_array_big.append(avT_k)
        av_pr_array.append(average_T_pr)
        av_cs_array.append(average_T_cs)
        av_hs_array.append(average_fin_hs)
        # plt.plot(av_array_big, diff_array_big)
        # plt.show()
        
#%%
import pandas as pd
import numpy as np
from Plotting import save



title = 'Initial_test_3p9_start'

df = pd.DataFrame(np.matrix([diff_array_big, av_array_big, av_pr_array, av_cs_array, av_hs_array]).transpose())
df.columns = ["diff_array_big", "av_array_big", "av_pr_array", "av_cs_array", "av_hs_array"]
df.to_csv(f'Optimised_T_k_matrix_method/data_{title}', index=False)


save(T_k_pr_update, T_k_cs_update, T_k_hs_update, T_k_fn_update, title, "Optimised_T_k_matrix_method")

        

#%%

