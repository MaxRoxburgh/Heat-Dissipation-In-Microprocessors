from Initializing import full_system_initialisation, heat_sinc_initialisation, np
from Updates import full_matrix_iteration_method
import pandas as pd

q = 0.5*1e9 # W/m^3
T_a = 293 # K
K_si = 150 # W/mK
K_ceramic = 230 # W/mK
K_alu = 250 # W/mK
h = (1/3)*1e-3 # m         

M_inv_pr = np.matrix(pd.read_csv('inverses/M_inv_pr', header=None))
M_inv_cs = np.matrix(pd.read_csv('inverses/M_inv_cs', header=None))
M_inv_fn = np.matrix(pd.read_csv('inverses/M_inv_fn', header=None))

# for the different heat sinc sizes
n_fins_list = np.array([12])
ambients = [1.1]

for n_fins, ambient in zip(n_fins_list, ambients):
    
    path = 'Donotsave'
    title = f'{n_fins}_fins_5p4_s3'
    
    print(f'Run for {n_fins} fins')

    T_0_hs, T_k_hs, M_inv_hs, vol_sinc = heat_sinc_initialisation(n_fins, ambient)
    T_0_pr, T_0_cs, T_0_fn, T_k_pr, T_k_cs, T_k_fn, vol_pr, vol_cs, vol_fins = full_system_initialisation(ambient, int(n_fins/2), scaling=3, get_heat_sinc=False)
    
    full_matrix_iteration_method(T_0_pr, T_0_cs, T_0_hs, T_0_fn, T_k_pr, T_k_cs, T_k_hs, T_k_fn, h, M_inv_pr
                                 , M_inv_cs, M_inv_hs, M_inv_fn, vol_pr, vol_cs, vol_sinc, vol_fins
                                 , path, title, iterations=10000, nat=False)