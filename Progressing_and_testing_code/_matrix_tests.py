from MatrixFunctions import BC_vector_creator, gaus_jordan, matrix_maker, update_T_matrix_method    
import matplotlib.pyplot as plt   
from Initializing import initialise_boundary
import numpy as np


#%%
a = np.matrix([[4,4,4,4,4,4,4,4,4,4],
                [3.99,0,0,0,0,0,0,0,0,3.99],
                [3.99,0,0,0,0,0,0,0,0,3.99],
                [3.99,0,0,0,0,0,0,0,0,3.99],
                [3.99,0,0,0,0,0,0,0,0,3.99],
                [3.99,3.99,4,4,4,4,4,4,3.99,3.99]], dtype='float')

# ambient = 4
# top = np.array([ambient for i in range(60)])
# bottom = top.copy()
# left = np.array([ambient for i in range(6)])
# right = left.copy()

# a, T_k_cs = initialise_boundary(top, left, right, bottom, ambient,[ambient,ambient,ambient,ambient])
# M_inv_cs = np.matrix(pd.read_csv('inverses/M_inv_cs', header=None))
# a = T_0_cs_update

r, c = np.shape(a)
#%%
BC_M = BC_vector_creator(a)
print(BC_M)
#%%
plt.imshow(BC_M, 'plasma')
plt.show()
#%%
M = matrix_maker(BC_M, False)
plt.imshow(M, 'plasma')
plt.show()
inv_M = gaus_jordan(M)
plt.imshow(inv_M, 'plasma')
plt.show()
M_k = update_T_matrix_method(a, inv_M)
plt.imshow(M_k, 'plasma', vmin=3.7, vmax=4.1)
plt.show()

#%%
import pandas as pd
M_inv_cs = np.matrix(pd.read_csv('inverses/M_inv_cs', header=None))

check = inv_M - M_inv_cs
plt.imshow(check, vmin=0, vmax = 0.01)


#%%
ambient = 1

average_temp_under = []
ks = [i for i in range(20,70,3)]

for k in range(6,20):
    top = np.array([ambient for i in range(k)])
    bottom = top.copy()
    left = np.array([ambient for i in range(k)])
    right = left.copy()
    
    a, T_k_pr = initialise_boundary(top, left, right, bottom, ambient,[ambient,ambient,ambient,ambient])
    r, c = np.shape(a)
    plt.imshow(a)
    plt.show()
    BC_M = BC_vector_creator(a).transpose().reshape(r-2,c-2)
    M = matrix_maker(BC_M, False)
    inv_M = gaus_jordan(M)
    wb, wob = update_T_matrix_method(a, inv_M, 10, 1 * (20/k))
    plt.imshow(wb) 
    plt.show()
    
    average_temp_under.append(sum(wob.A1)/len(wob.A1))
    
#%%
for i in average_temp:
    average_temp_under.append(i)
#%%

ks1 = [i for i in range(6,20)]
for i in ks:
    ks1.append(i)

#%%
import pandas as pd
plt.plot(ks1, average_temp_under)
Matri = np.matrix([ks1, average_temp_under]).transpose()
df = pd.DataFrame(Matri)
df.columns = ['square size', 'average temp']
df.to_csv('inverses/consistency_analysis', index=False)
    
#%%   
row = [1 for i in range(14*3)]
M_pr = np.matrix([row for i in range(3)])
M_pr = matrix_maker(M_pr, False)

row = [1 for i in range(60)]
M_cs = np.matrix([row for i in range(6)])
M_cs = matrix_maker(M_cs, False)

row = [1 for i in range(34*3)]
M_hs = np.matrix([row for i in range(12)])
M_hs = matrix_maker(M_hs, False)

row = [1 for i in range(3)]
M_fn = np.matrix([row for i in range(90)])
M_fn = matrix_maker(M_fn, False)

M_inv_pr = gaus_jordan(M_pr)
M_inv_cs = gaus_jordan(M_cs)
M_inv_hs = gaus_jordan(M_hs)
M_inv_fn = gaus_jordan(M_fn)
#%%
import pandas as pd

df = pd.DataFrame(M_inv_cs)
df.to_csv('inverses/M_inv_cs', header=False, index=False)

df = pd.DataFrame(M_inv_pr)
df.to_csv('inverses/M_inv_pr', header=False, index=False)

df = pd.DataFrame(M_inv_hs)
df.to_csv('inverses/M_inv_hs', header=False, index=False)

df = pd.DataFrame(M_inv_fn)
df.to_csv('inverses/M_inv_fn', header=False, index=False)

#%% testing inverse matrix method

from Calculators import average_temp
from Updates import update_boundary_convection
import pandas as pd
from MatrixFunctions import matrix_maker, gaus_jordan
from Plotting import trim


ambient = 20
q = 0.5*1e9 # W/m^3
T_a = 293 # K
K_si = 150 # W/mK
K_ceramic = 230 # W/mK
K_alu = 250 # W/mK
h = (1/3)*1e-3 # m  

scaling = [3,4,5,10]
T_k_list = np.array([initialise_boundary(np.array([ambient for i in range(14*i)]), np.array([ambient for i in range(i)]), np.array([ambient for i in range(i)]), np.array([ambient for i in range(14*i)]), ambient,[ambient,ambient,ambient,ambient])[1] for i in scaling], dtype = 'object')
M_inv = np.array([gaus_jordan(matrix_maker(np.matrix(i))) for i in T_k_list], dtype='object')


#%%
for i, j in zip(scaling, M_inv):
    h = (1/i)*1e-3

    top = np.array([ambient for i in range(14*i)])
    bottom = top.copy()
    left = np.array([ambient for i in range(i)])
    right = left.copy()
    
    T_0, T_k = initialise_boundary(top, left, right, bottom, ambient,[ambient,ambient,ambient,ambient])
    # M = matrix_maker(T_k)
    M_inv_pr = j
    
    
    estimate = (7e6/(30*1.31))**(3/4)/293 + 1
    print(f'\n\nExpected estimate of T:   {estimate}')
    print(f'Average starting T:       {average_temp(T_k)}')
    
    for z in range(10000*i**2):
        T_0, T_k = update_boundary_convection(T_k, 150), update_T_matrix_method(T_0, M_inv_pr, q, h)
    
    av_T = sum(trim(T_k).A1)/len(trim(T_k).A1)
    print(f'\nScaling = {i}')    
    print(f'Average T after 10k iter: {av_T}')
    fig, ax = plt.subplots()
    im = ax.imshow(T_k, cmap='plasma')
    cbar = ax.figure.colorbar(im, orientation='horizontal')
    cbar.set_label('Temperature Scale', loc='center')
    plt.show()

#%% update T jacobi method

from MatrixFunctions import update_T_k_jacobi_matrix, find_D_Tl_Tr
from Calculators import average_temp
from Updates import update_boundary_convection
import pandas as pd
from MatrixFunctions import matrix_maker, gaus_jordan
from Plotting import trim
from Updates import update_with_source

scaling = [3,4,5,6,10]
estimate = (7e6/(30*1.31))**(3/4)/293 + 1
print(f'\nRough estimate of T:      {estimate}')
print('Average starting T:           23.18')
 
iterations = 2000
q = 0.5*1e9 # W/m^3
T_a = 293 # K
K_si = 150 # W/mK
K_ceramic = 230 # W/mK
K_alu = 250 # W/mK
# h = (1/3)*1e-3 # m  

# scaling = [3]
iterations = 50000
for i in scaling:
    h = (1/i)*1e-3
    ambient = 23.18

    top = np.array([ambient for i in range(14*i)])
    bottom = top.copy()
    left = np.array([ambient for i in range(i)])
    right = left.copy()
    
    T_0, T_k = initialise_boundary(top, left, right, bottom, ambient,[ambient,ambient,ambient,ambient])
    M = matrix_maker(T_k)
    D, L, R = find_D_Tl_Tr(M)
    M = L + R
    
    print('\n\n\n----------------------------------------\n\n\n')
    
    
    for z in range(iterations):
        T_0, T_k = update_boundary_convection(T_k, 150, 293, h), update_T_k_jacobi_matrix(T_k, T_0, M, q, h)
        if z% 10000 == 0:
            plt.imshow(T_k)
            plt.show()
    
    av_T = sum(trim(T_k).A1)/len(trim(T_k).A1)
    print(f'Scaling = {i}:\n')    
    print(f'Average T jacobi matrix method: {av_T}')
    fig, ax = plt.subplots()
    im = ax.imshow(T_k, cmap='plasma')
    cbar = ax.figure.colorbar(im, orientation='horizontal')
    cbar.set_label('Temperature Scale', loc='center')
    plt.show()
    
    
    T_0, T_k = initialise_boundary(top, left, right, bottom, ambient,[ambient,ambient,ambient,ambient])

    
    for z in range(iterations):
        T_0, T_k = update_boundary_convection(T_k, 150, 293, h), update_with_source(T_k, T_0, q, h)
        if z% 10000 == 0:
            plt.imshow(T_k)
            plt.show()
    
    av_T = sum(trim(T_k).A1)/len(trim(T_k).A1) 
    print(f'Average T pictoral operator method: {av_T}')
    fig, ax = plt.subplots()
    im = ax.imshow(T_k, cmap='plasma')
    cbar = ax.figure.colorbar(im, orientation='horizontal')
    cbar.set_label('Temperature Scale', loc='center')
    plt.show()
    
    T_0, T_k = initialise_boundary(top, left, right, bottom, ambient,[ambient,ambient,ambient,ambient])
    M = matrix_maker(T_k)
    M_inv_pr = gaus_jordan(M)
    
    for z in range(iterations):
        T_0, T_k = update_boundary_convection(T_k, 150, 293, h), update_T_matrix_method(T_0, M_inv_pr, q, h)
        if z% 10000 == 0:
            plt.imshow(T_k)
            plt.show()
    
    av_T = sum(trim(T_k).A1)/len(trim(T_k).A1)
    print(f'Average T inverse matrix method: {av_T}')
    fig, ax = plt.subplots()
    im = ax.imshow(T_k, cmap='plasma')
    cbar = ax.figure.colorbar(im, orientation='horizontal')
    cbar.set_label('Temperature Scale', loc='center')
    plt.show()

#%%

estimate = (7e6/(30*1.31))**(3/4)/293 + 1
print(f'\nRough estimate of T:      {estimate}')
print('Average starting T:           23.184')
 
q = 0.5*1e9 # W/m^3
T_a = 293 # K
K_si = 150 # W/mK
K_ceramic = 230 # W/mK
K_alu = 250 # W/mK
# h = (1/3)*1e-3 # m  

scaling = [3]
iterations = 50000
for i in scaling:
    h = (1/i)*1e-3
    ambient = 23.184
    
    top = np.array([ambient for i in range(14*i)])
    bottom = top.copy()
    left = np.array([ambient for i in range(i)])
    right = left.copy()
    
    print('\n\n\n----------------------------------------\n\n\n')
    print(f'Scaling = {i}:\n')    
    
    
    T_0, T_k = initialise_boundary(top, left, right, bottom, ambient,[ambient,ambient,ambient,ambient])
    M = matrix_maker(T_k)
    M_inv_pr = gaus_jordan(M)
    
    for z in range(iterations*4):
        T_0, T_k = update_boundary_convection(T_k, 150, 293, h), update_T_matrix_method(T_0, M_inv_pr, q, h)
        # if z% 10000 == 0:
        #     plt.imshow(T_k)
        #     plt.show()
    
    av_T = sum(trim(T_k).A1)/len(trim(T_k).A1)
    print(f'Average T inverse matrix method: {av_T}')
    fig, ax = plt.subplots()
    im = ax.imshow(T_k, cmap='plasma')
    cbar = ax.figure.colorbar(im, orientation='horizontal')
    cbar.set_label('Temperature Scale', loc='center')
    plt.show()

    


    T_0_j, T_k_j = T_0.copy(), T_k.copy()
    M = matrix_maker(T_k)
    D, L, R = find_D_Tl_Tr(M)
    M = L + R
    
    
    for z in range(iterations):
        T_0_j, T_k_j = update_boundary_convection(T_k_j, 150, 293, h), update_T_k_jacobi_matrix(T_k_j, T_0_j, M, q, h)
        # if z% 15000 == 0:
        #     plt.imshow(T_k)
        #     plt.show()
    
    av_T = sum(trim(T_k_j).A1)/len(trim(T_k_j).A1)
    print(f'Average T jacobi matrix method: {av_T}')
    fig, ax = plt.subplots()
    im = ax.imshow(T_k_j, cmap='plasma')
    cbar = ax.figure.colorbar(im, orientation='horizontal')
    cbar.set_label('Temperature Scale', loc='center')
    plt.show()
    
    
    T_0_p, T_k_p = T_0.copy(), T_k.copy()

    
    for z in range(iterations):
        T_0_p, T_k_p = update_boundary_convection(T_k_p, 150, 293, h), update_with_source(T_k_p, T_0_p, q, h)
        # if z% 10000 == 0:
        #     plt.imshow(T_k)
        #     plt.show()
    
    av_T = sum(trim(T_k_p).A1)/len(trim(T_k_p).A1) 
    print(f'Average T pictoral operator method: {av_T}')
    fig, ax = plt.subplots()
    im = ax.imshow(T_k_p, cmap='plasma')
    cbar = ax.figure.colorbar(im, orientation='horizontal')
    cbar.set_label('Temperature Scale', loc='center')
    plt.show()
    

#%%
top = np.array([ambient for i in range(14*10)])
bottom = top.copy()
left = np.array([ambient for i in range(10)])
right = left.copy()

T_0, T_k = initialise_boundary(top, left, right, bottom, ambient,[ambient,ambient,ambient,ambient])
M = matrix_maker(T_k)
M_inv_pr = gaus_jordan(M)

#%%
top = np.array([ambient for i in range(14*15)])
bottom = top.copy()
left = np.array([ambient for i in range(15)])
right = left.copy()

T_0, T_k = initialise_boundary(top, left, right, bottom, ambient,[ambient,ambient,ambient,ambient])
M = matrix_maker(T_k)
M_inv_pr = gaus_jordan(M)


#%%
ambient = 29.5
q=0.5e9
iterations = 50000

scaling = [15]

for i in scaling:
    h = (1/i)*1e-3
    top = np.array([ambient for i in range(14*i)])
    bottom = top.copy()
    left = np.array([ambient for i in range(i)])
    right = left.copy()
    
    T_0, T_k = initialise_boundary(top, left, right, bottom, ambient,[ambient,ambient,ambient,ambient])
    
    # M = matrix_maker(T_k)
    # M_inv_pr = gaus_jordan(M)
    for n in range(iterations):
        T_0, T_k = update_boundary_convection(T_k, 150, 293, h), update_T_matrix_method(T_0, M_inv_pr, q, h)
        if n%7500 == 0:
            print(average_temp(T_k), f'iteration {n}')
            fig, ax = plt.subplots()
            im = ax.imshow(T_k, cmap='plasma')
            cbar = ax.figure.colorbar(im, orientation='horizontal')
            cbar.set_label('Temperature Scale', loc='center')
            plt.show()

            
    
    # print(average_temp(T_k))
    # fig, ax = plt.subplots()
    # im = ax.imshow(T_k, cmap='plasma')
    # cbar = ax.figure.colorbar(im, orientation='horizontal')
    # cbar.set_label('Temperature Scale', loc='center')
    # plt.show()


#%%

from MatrixFunctions import matrix_maker_2, gaus_jordan
import matplotlib.pyplot as plt   

a = np.ones((9,9), dtype='float')
b = matrix_maker_2(a)
b_inv = gaus_jordan(b)

#%%
fig, (ax1, ax2) = plt.subplots(1,2)

# vmin, vmax = -1.7, 0.4
cmap = 'jet_r' 

im1 = ax1.imshow(b, cmap)#, vmin=vmin, vmax=vmax)
im2 = ax2.imshow(b_inv, cmap)#, vmin=vmin, vmax=vmax)
# im3 = ax3.imshow(np.matrix(a.flatten()).transpose())
ax1.set_title('a)')
ax2.set_title('b)')
# ax3.set_title('c)')

# cbar1 = ax1.figure.colorbar(im1, orientation='horizontal')
# cbar2 = ax2.figure.colorbar(im2, orientation='horizontal')
# cbar.set_label('Temperature Scale', loc='center')


ax1.axes.xaxis.set_visible(False)
ax1.axes.yaxis.set_visible(False)
ax2.axes.xaxis.set_visible(False)
ax2.axes.yaxis.set_visible(False)
# ax3.axes.xaxis.set_visible(False)
# ax3.axes.yaxis.set_visible(False)
#%%

fig, ax = plt.subplots(1,1)
image = ax.imshow(b_inv, 'magma_r')
cbar = ax.figure.colorbar(image, orientation='vertical')
cbar.set_label('Temp', loc='center')
ax.axes.xaxis.set_visible(False)
ax.axes.yaxis.set_visible(False)













