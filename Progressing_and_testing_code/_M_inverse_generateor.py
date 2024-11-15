from MatrixFunctions import matrix_maker, gaus_jordan
import numpy as np
import time
import pandas as pd
#%%
t = []
scalings = [3,4,5,6,7,8,9,10,11]

t.append(time.time()/60)

for scaling in scalings:
    
    print(f'\n\n\nInverses for scaling: {scaling}')

    T_0_pr, T_0_cs, T_0_fn = np.matrix(np.ones((1*scaling + 2,14*scaling + 2))),np.matrix(np.ones((2*scaling + 2,20*scaling + 2))), np.matrix(np.ones((30*scaling + 2,1*scaling + 2)))
    
    print('\nDone with: M_inv_hs_34')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    M_pr = matrix_maker(T_0_pr)
    M_cs = matrix_maker(T_0_cs)
    
    M_inv_pr = gaus_jordan(M_pr)
    print('\nDone with: M_inv_pr\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    M_inv_cs = gaus_jordan(M_cs)
    print('\nDone with: M_inv_cs\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    M_fn_30 = matrix_maker(T_0_fn)
    M_inv_fn_30 = gaus_jordan(M_fn_30)
    print('\nDone with: M_inv_fn_30\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    T_k_fn_32 = np.matrix(np.ones((32*scaling + 2,1*scaling + 2)))
    M_fn_32 = matrix_maker(T_k_fn_32)
    M_inv_fn_32 = gaus_jordan(M_fn_32)
    print('\nDone with: M_inv_fn_32\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    T_k_fn_34 = np.matrix(np.ones((34*scaling + 2,1*scaling + 2)))
    M_fn_34 = matrix_maker(T_k_fn_34)
    M_inv_fn_34 = gaus_jordan(M_fn_34)
    print('\nDone with: M_inv_fn_34\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    T_k_fn_36 = np.matrix(np.ones((36*scaling + 2,1*scaling + 2)))
    M_fn_36 = matrix_maker(T_k_fn_36)
    M_inv_fn_36 = gaus_jordan(M_fn_36)
    print('\nDone with: M_inv_fn_36\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    T_k_fn_38 = np.matrix(np.ones((38*scaling + 2,1*scaling + 2)))
    M_fn_38 = matrix_maker(T_k_fn_38)
    M_inv_fn_38 = gaus_jordan(M_fn_38)
    print('\nDone with: M_inv_fn_38\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    T_k_fn_40 = np.matrix(np.ones((40*scaling + 2,1*scaling + 2)))
    M_fn_40 = matrix_maker(T_k_fn_40)
    M_inv_fn_40 = gaus_jordan(M_fn_40)
    print('\nDone with: M_inv_fn_40\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    T_k_fn_42 = np.matrix(np.ones((42*scaling + 2,1*scaling + 2)))
    M_fn_42 = matrix_maker(T_k_fn_42)
    M_inv_fn_42 = gaus_jordan(M_fn_42)
    print('\nDone with: M_inv_fn_42\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    T_k_fn_44 = np.matrix(np.ones((44*scaling + 2,1*scaling + 2)))
    M_fn_44 = matrix_maker(T_k_fn_44)
    M_inv_fn_44 = gaus_jordan(M_fn_44)
    print('\nDone with: M_inv_fn_44\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    T_k_hs_22 = np.matrix(np.ones((4*scaling + 2, 22*scaling + 2)))
    M_hs_22 = matrix_maker(T_k_hs_22)
    M_inv_hs_22 = gaus_jordan(M_hs_22)
    print('\nDone with: M_inv_hs_22\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    T_k_hs_28 = np.matrix(np.ones((4*scaling + 2, 28*scaling + 2)))
    M_hs_28 = matrix_maker(T_k_hs_22)
    M_inv_hs_28 = gaus_jordan(M_hs_22)
    print('\nDone with: M_inv_hs_28\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    T_k_hs_34 = np.matrix(np.ones((4*scaling + 2, 34*scaling + 2)))
    M_hs_34 = matrix_maker(T_k_hs_34)
    M_inv_hs_34 = gaus_jordan(M_hs_34)
    print('\nDone with: M_inv_hs_34\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    T_k_hs_40 = np.matrix(np.ones((4*scaling + 2, 40*scaling + 2)))
    M_hs_40 = matrix_maker(T_k_hs_40)
    M_inv_hs_40 = gaus_jordan(M_hs_40)
    print('\nDone with: M_inv_hs_40\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    T_k_hs_46 = np.matrix(np.ones((4*scaling + 2, 46*scaling + 2)))
    M_hs_46 = matrix_maker(T_k_hs_46)
    M_inv_hs_46 = gaus_jordan(M_hs_46)
    print('\nDone with: M_inv_hs_46\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    T_k_hs_52 = np.matrix(np.ones((4*scaling + 2, 52*scaling + 2)))
    M_hs_52 = matrix_maker(T_k_hs_52)
    M_inv_hs_52 = gaus_jordan(M_hs_52)
    print('\nDone with: M_inv_hs_52\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    df = pd.DataFrame(M_inv_pr)
    df.to_csv(f"M_inverses/M_inv_pr_s{scaling}", index=False, header=False)
    
    df = pd.DataFrame(M_inv_cs)
    df.to_csv(f"M_inverses/M_inv_cs_s{scaling}", index=False, header=False)
    
    ls = [30,32,34,36,38,40,42,44]
    Ms = [M_inv_fn_30,M_inv_fn_32,M_inv_fn_34,M_inv_fn_36,M_inv_fn_38,M_inv_fn_40,M_inv_fn_42,M_inv_fn_44]
    for i, j in zip(ls, Ms):
        df = pd.DataFrame(j)
        df.to_csv(f"M_inverses/M_inv_fn_{i}_s{scaling}", index=False, header=False)
        
    ls = [22,28,34, 40, 46]
    Ms = [M_inv_hs_22,M_inv_hs_28,M_inv_hs_34,M_inv_hs_40,M_inv_hs_46]
    for i, j in zip(ls, Ms):
        df = pd.DataFrame(j)
        df.to_csv(f"M_inverses/M_inv_hs_{i}_s{scaling}", index=False, header=False)
#%%
from MatrixFunctions import matrix_maker_2
t = []
scalings = [3,4,5,6,7,8,9]

t.append(time.time()/60)

for scaling in scalings:
    
    print(f'\n\n\nInverses for scaling: {scaling}')

    T_0_pr, T_0_cs, T_0_fn = np.matrix(np.ones((1*scaling + 2,14*scaling + 2))),np.matrix(np.ones((2*scaling + 2,20*scaling + 2))), np.matrix(np.ones((30*scaling + 2,1*scaling + 2)))
    
    print('\nDone with: M_inv_hs_34')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    M_pr = matrix_maker_2(T_0_pr)
    M_cs = matrix_maker_2(T_0_cs)
    
    M_inv_pr = gaus_jordan(M_pr)
    print('\nDone with: M_inv_pr\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    M_inv_cs = gaus_jordan(M_cs)
    print('\nDone with: M_inv_cs\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    M_fn_30 = matrix_maker_2(T_0_fn)
    M_inv_fn_30 = gaus_jordan(M_fn_30)
    print('\nDone with: M_inv_fn_30\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    T_k_fn_32 = np.matrix(np.ones((32*scaling + 2,1*scaling + 2)))
    M_fn_32 = matrix_maker_2(T_k_fn_32)
    M_inv_fn_32 = gaus_jordan(M_fn_32)
    print('\nDone with: M_inv_fn_32\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    T_k_fn_34 = np.matrix(np.ones((34*scaling + 2,1*scaling + 2)))
    M_fn_34 = matrix_maker_2(T_k_fn_34)
    M_inv_fn_34 = gaus_jordan(M_fn_34)
    print('\nDone with: M_inv_fn_34\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    T_k_fn_36 = np.matrix(np.ones((36*scaling + 2,1*scaling + 2)))
    M_fn_36 = matrix_maker_2(T_k_fn_36)
    M_inv_fn_36 = gaus_jordan(M_fn_36)
    print('\nDone with: M_inv_fn_36\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    T_k_fn_38 = np.matrix(np.ones((38*scaling + 2,1*scaling + 2)))
    M_fn_38 = matrix_maker_2(T_k_fn_38)
    M_inv_fn_38 = gaus_jordan(M_fn_38)
    print('\nDone with: M_inv_fn_38\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    T_k_fn_40 = np.matrix(np.ones((40*scaling + 2,1*scaling + 2)))
    M_fn_40 = matrix_maker_2(T_k_fn_40)
    M_inv_fn_40 = gaus_jordan(M_fn_40)
    print('\nDone with: M_inv_fn_40\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    T_k_fn_42 = np.matrix(np.ones((42*scaling + 2,1*scaling + 2)))
    M_fn_42 = matrix_maker_2(T_k_fn_42)
    M_inv_fn_42 = gaus_jordan(M_fn_42)
    print('\nDone with: M_inv_fn_42\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    T_k_fn_44 = np.matrix(np.ones((44*scaling + 2,1*scaling + 2)))
    M_fn_44 = matrix_maker_2(T_k_fn_44)
    M_inv_fn_44 = gaus_jordan(M_fn_44)
    print('\nDone with: M_inv_fn_44\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    T_k_hs_22 = np.matrix(np.ones((4*scaling + 2, 22*scaling + 2)))
    M_hs_22 = matrix_maker_2(T_k_hs_22)
    M_inv_hs_22 = gaus_jordan(M_hs_22)
    print('\nDone with: M_inv_hs_22\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    T_k_hs_28 = np.matrix(np.ones((4*scaling + 2, 28*scaling + 2)))
    M_hs_28 = matrix_maker_2(T_k_hs_28)
    M_inv_hs_28 = gaus_jordan(M_hs_28)
    print('\nDone with: M_inv_hs_28\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    T_k_hs_34 = np.matrix(np.ones((4*scaling + 2, 34*scaling + 2)))
    M_hs_34 = matrix_maker_2(T_k_hs_34)
    M_inv_hs_34 = gaus_jordan(M_hs_34)
    print('\nDone with: M_inv_hs_34\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    T_k_hs_40 = np.matrix(np.ones((4*scaling + 2, 40*scaling + 2)))
    M_hs_40 = matrix_maker_2(T_k_hs_40)
    M_inv_hs_40 = gaus_jordan(M_hs_40)
    print('\nDone with: M_inv_hs_40\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    T_k_hs_46 = np.matrix(np.ones((4*scaling + 2, 46*scaling + 2)))
    M_hs_46 = matrix_maker_2(T_k_hs_46)
    M_inv_hs_46 = gaus_jordan(M_hs_46)
    print('\nDone with: M_inv_hs_46\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    T_k_hs_52 = np.matrix(np.ones((4*scaling + 2, 52*scaling + 2)))
    M_hs_52 = matrix_maker_2(T_k_hs_52)
    M_inv_hs_52 = gaus_jordan(M_hs_52)
    print('\nDone with: M_inv_hs_52\n\n')
    t.append(time.time()/60)
    print('time elapsed:', t[-1]-t[-2], '\n\n')
    
    df = pd.DataFrame(M_inv_pr)
    df.to_csv(f"M_inverses/M_inv_pr_s{scaling}_FDO4", index=False, header=False)
    
    df = pd.DataFrame(M_inv_cs)
    df.to_csv(f"M_inverses/M_inv_cs_s{scaling}_FDO4", index=False, header=False)
    
    ls = [30,32,34,36,38,40,42,44]
    Ms = [M_inv_fn_30,M_inv_fn_32,M_inv_fn_34,M_inv_fn_36,M_inv_fn_38,M_inv_fn_40,M_inv_fn_42,M_inv_fn_44]
    for i, j in zip(ls, Ms):
        df = pd.DataFrame(j)
        df.to_csv(f"M_inverses/M_inv_fn_{i}_s{scaling}_FDO4", index=False, header=False)
        
    ls = [22,28,34, 40, 46]
    Ms = [M_inv_hs_22,M_inv_hs_28,M_inv_hs_34,M_inv_hs_40,M_inv_hs_46]
    for i, j in zip(ls, Ms):
        df = pd.DataFrame(j)
        df.to_csv(f"M_inverses/M_inv_hs_{i}_s{scaling}_FDO4", index=False, header=False)

