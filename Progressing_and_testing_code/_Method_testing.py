from Initializing import initialise_boundary
from Updates import update_boundary_convection
from Calculators import average_temp
import numpy as np
import time

K=150
q=0.5e9
h=(1/6)*1e-3

timeings = []

ambient = 30
top = np.array([ambient for i in range(14*6)])
bottom = top.copy()
left = np.array([ambient for i in range(6)])
right = left.copy()

T_0_ab, T_k_ab = initialise_boundary(top, left, right, bottom, ambient)

ambient = 20
top = np.array([ambient for i in range(14*6)])
bottom = top.copy()
left = np.array([ambient for i in range(6)])
right = left.copy()

T_0_bl, T_k_bl = initialise_boundary(top, left, right, bottom, ambient)

#         test GS FDO2 above, below, at eq
from MatrixFunctions import update_T_k_GS_matrix, matrix_maker, find_D_Tl_Tr, gaus_jordan


d_ab_GS_O2 = []
d_bl_GS_O2 = []
T_ab_GS_O2 = []
T_bl_GS_O2 = []

T_0_GS_O2_ab = T_0_ab.copy()
T_k_GS_O2_ab = T_k_ab.copy()

T_0_GS_O2_bl = T_0_bl.copy()
T_k_GS_O2_bl = T_k_bl.copy()

M = matrix_maker(T_k_GS_O2_ab)
D, Tl, Tu = find_D_Tl_Tr(M)
M = D+Tl
M_inv = gaus_jordan(M)


T_ab = 30
T_bl = 20
t_before = time.time()
for n in range(20):
    for m in range(500):
        T_0_GS_O2_ab, T_k_GS_O2_ab = update_boundary_convection(T_k_GS_O2_ab, K, h=h), update_T_k_GS_matrix(T_k_GS_O2_ab, T_0_GS_O2_ab, M_inv, Tu, q=q,h=h)
        T_0_GS_O2_bl, T_k_GS_O2_bl = update_boundary_convection(T_k_GS_O2_bl, K, h=h), update_T_k_GS_matrix(T_k_GS_O2_bl, T_0_GS_O2_bl, M_inv, Tu, q=q,h=h)
    T_ab_next = average_temp(T_k_GS_O2_ab)
    T_bl_next = average_temp(T_k_GS_O2_bl)
    
    T_ab_GS_O2.append(T_ab_next)
    T_bl_GS_O2.append(T_bl_next)
    
    d_ab_GS_O2.append(T_ab-T_ab_next)
    d_bl_GS_O2.append(T_bl-T_bl_next)
    
    T_ab = T_ab_next
    T_bl = T_bl_next
t_after = time.time()
timeings.append(t_after-t_before)

print(f'\n\nDone w/ GSO2')
print(f'time taken: {timeings[len(timeings)-2]}')
#         test GS FDO4 above, below, at eq
from MatrixFunctions import matrix_maker_2


d_ab_GS_O4 = []
d_bl_GS_O4 = []
T_ab_GS_O4 = []
T_bl_GS_O4 = []

T_0_GS_O4_ab = T_0_ab.copy()
T_k_GS_O4_ab = T_k_ab.copy()

T_0_GS_O4_bl = T_0_bl.copy()
T_k_GS_O4_bl = T_k_bl.copy()

M = matrix_maker_2(T_k_GS_O2_ab)
D, Tl, Tu = find_D_Tl_Tr(M)
M = D+Tl
M_inv = gaus_jordan(M)


T_ab = 30
T_bl = 20
t_before = time.time()
for n in range(20):
    for m in range(500):
        T_0_GS_O4_ab, T_k_GS_O4_ab = update_boundary_convection(T_k_GS_O4_ab, K, h=h), update_T_k_GS_matrix(T_k_GS_O4_ab, T_0_GS_O4_ab, M_inv, Tu, q=q,h=h)
        T_0_GS_O4_bl, T_k_GS_O4_bl = update_boundary_convection(T_k_GS_O4_bl, K, h=h), update_T_k_GS_matrix(T_k_GS_O4_bl, T_0_GS_O4_bl, M_inv, Tu, q=q,h=h)
    T_ab_next = average_temp(T_k_GS_O4_ab)
    T_bl_next = average_temp(T_k_GS_O4_bl)
    
    T_ab_GS_O4.append(T_ab_next)
    T_bl_GS_O4.append(T_bl_next)
    
    d_ab_GS_O4.append(T_ab-T_ab_next)
    d_bl_GS_O4.append(T_bl-T_bl_next)
    
    T_ab = T_ab_next
    T_bl = T_bl_next
t_after = time.time()
timeings.append(t_after-t_before)
print(f'\n\nDone w/ GSO4')
print(f'time taken: {timeings[len(timeings)-2]}')
#         test pictoral above, below, at eq
from Updates import update_with_source

d_ab_pic = []
d_bl_pic = []
T_ab_pic = []
T_bl_pic = []

T_0_pic_ab = T_0_ab.copy()
T_k_pic_ab = T_k_ab.copy()

T_0_pic_bl = T_0_bl.copy()
T_k_pic_bl = T_k_bl.copy()

T_ab = 30
T_bl = 20
t_before = time.time()
for n in range(20):
    for m in range(500):
        T_0_pic_ab, T_k_pic_ab = update_boundary_convection(T_k_pic_ab, K, h=h), update_with_source(T_k_pic_ab, T_0_pic_ab, q=q,h=h)
        T_0_pic_bl, T_k_pic_bl = update_boundary_convection(T_k_pic_bl, K, h=h), update_with_source(T_k_pic_bl, T_0_pic_bl, q=q,h=h)
    T_ab_next = average_temp(T_k_pic_ab)
    T_bl_next = average_temp(T_k_pic_bl)
    
    T_ab_pic.append(T_ab_next)
    T_bl_pic.append(T_bl_next)
    
    d_ab_pic.append(T_ab-T_ab_next)
    d_bl_pic.append(T_bl-T_bl_next)
    
    T_ab = T_ab_next
    T_bl = T_bl_next


t_after = time.time()
timeings.append(t_after-t_before)
print(f'\n\nDone w/ pictoral')
print(f'time taken: {timeings[len(timeings)-2]}')
#         test jacobi FDO2 above, below, at eq
from MatrixFunctions import update_T_k_jacobi_matrix


d_ab_JC_O2 = []
d_bl_JC_O2 = []
T_ab_JC_O2 = []
T_bl_JC_O2 = []

T_0_JC_O2_ab = T_0_ab.copy()
T_k_JC_O2_ab = T_k_ab.copy()

T_0_JC_O2_bl = T_0_bl.copy()
T_k_JC_O2_bl = T_k_bl.copy()

M = matrix_maker(T_k_GS_O2_ab)
D, Tl, Tu = find_D_Tl_Tr(M)
M = Tu+Tl


T_ab = 30
T_bl = 20
t_before = time.time()
for n in range(20):
    for m in range(500):
        T_0_JC_O2_ab, T_k_JC_O2_ab = update_boundary_convection(T_k_JC_O2_ab, K, h=h), update_T_k_jacobi_matrix(T_k_JC_O2_ab, T_0_JC_O2_ab, M ,q=q,h=h)
        T_0_JC_O2_bl, T_k_JC_O2_bl = update_boundary_convection(T_k_JC_O2_bl, K, h=h), update_T_k_jacobi_matrix(T_k_JC_O2_bl, T_0_JC_O2_bl, M, q=q,h=h)
    T_ab_next = average_temp(T_k_JC_O2_ab)
    T_bl_next = average_temp(T_k_JC_O2_bl)
    
    T_ab_JC_O2.append(T_ab_next)
    T_bl_JC_O2.append(T_bl_next)
    
    d_ab_JC_O2.append(T_ab-T_ab_next)
    d_bl_JC_O2.append(T_bl-T_bl_next)
    
    T_ab = T_ab_next
    T_bl = T_bl_next

t_after = time.time()
timeings.append(t_after-t_before)
print(f'\n\nDone w/ jacobi')
print(f'time taken: {timeings[len(timeings)-2]}')
#         test inverse matrix FDO2 above, below, at eq

from MatrixFunctions import update_T_matrix_method
import pandas as pd


d_ab_IM_02 = []
d_bl_IM_02 = []
T_ab_IM_02 = []
T_bl_IM_02 = []

T_0_IM_02_ab = T_0_ab.copy()
T_k_IM_02_ab = T_k_ab.copy()

T_0_IM_02_bl = T_0_bl.copy()
T_k_IM_02_bl = T_k_bl.copy()

M = np.matrix(pd.read_csv('M_inverses/M_inv_pr_s6', header=None))

T_ab = 30
T_bl = 20
t_before = time.time()
for n in range(20):
    for m in range(500):
        T_0_IM_02_ab, T_k_IM_02_ab = update_boundary_convection(T_k_IM_02_ab, K, h=h), update_T_matrix_method(T_0_IM_02_ab, M ,q=q,h=h)
        T_0_IM_02_bl, T_k_IM_02_bl = update_boundary_convection(T_k_IM_02_bl, K, h=h), update_T_matrix_method(T_0_IM_02_bl, M, q=q,h=h)
    T_ab_next = average_temp(T_k_IM_02_ab)
    T_bl_next = average_temp(T_k_IM_02_bl)
    
    T_ab_IM_02.append(T_ab_next)
    T_bl_IM_02.append(T_bl_next)
    
    d_ab_IM_02.append(T_ab-T_ab_next)
    d_bl_IM_02.append(T_bl-T_bl_next)
    
    T_ab = T_ab_next
    T_bl = T_bl_next


t_after = time.time()
timeings.append(t_after-t_before)
print(f'\n\nDone w/ Inverse O2')
print(f'time taken: {timeings[len(timeings)-2]}')
#         test inverse matrix FDO4 above, below, at eq

d_ab_IM_04 = []
d_bl_IM_04 = []
T_ab_IM_04 = []
T_bl_IM_04 = []

T_0_IM_04_ab = T_0_ab.copy()
T_k_IM_04_ab = T_k_ab.copy()

T_0_IM_04_bl = T_0_bl.copy()
T_k_IM_04_bl = T_k_bl.copy()

M = np.matrix(pd.read_csv('M_inverses/M_inv_pr_s6_FDO4', header=None))

T_ab = 30
T_bl = 20
t_before = time.time()
for n in range(20):
    for m in range(500):
        T_0_IM_04_ab, T_k_IM_04_ab = update_boundary_convection(T_k_IM_04_ab, K, h=h), update_T_matrix_method(T_0_IM_04_ab, M ,q=q,h=h)
        T_0_IM_04_bl, T_k_IM_04_bl = update_boundary_convection(T_k_IM_04_bl, K, h=h), update_T_matrix_method(T_0_IM_04_bl, M, q=q,h=h)
    T_ab_next = average_temp(T_k_IM_04_ab)
    T_bl_next = average_temp(T_k_IM_04_bl)
    
    T_ab_IM_04.append(T_ab_next)
    T_bl_IM_04.append(T_bl_next)
    
    d_ab_IM_04.append(T_ab-T_ab_next)
    d_bl_IM_04.append(T_bl-T_bl_next)
    
    T_ab = T_ab_next
    T_bl = T_bl_next


t_after = time.time()
timeings.append(t_after-t_before)
print(f'\n\nDone w/ Inverse O4')
print(f'time taken: {timeings[len(timeings)-2]}')


import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
plt.style.use('bmh')

def lin(x,a,b):
    return a*x + b


def plot_the_temp(x1, y1, x2, y2, label, colour, alpha, ax, labeling = True, start=0):
    if label == 'IM $O(h^4)$' or label == 'IM $O(h^2)$':
        lw = 0.8
        # return
    else:
        lw=3
    if labeling:
        ax.plot(x1, y1, color=f'{colour}', alpha=alpha, label=label, lw=lw)
    else:
        ax.plot(x1, y1, color=f'{colour}', alpha=alpha, lw=lw)
    ax.plot(x2, y2, color=f'{colour}', alpha=alpha, lw=lw)
    
    (a,b), cov = curve_fit(lin, np.append(x1,x2)[start:], np.append(y1,y2)[start:])
    ax.plot(np.append(x1,x2), [lin(i,a,b) for i in np.append(x1,x2)], alpha=0.2, lw=lw*1.25)
    ax.plot([-b/a],[0], 'o', color=colour)
    print(f'Converge point for {label} is {-b/a}')

# plt.plot(T_bl_IM_04, d_bl_IM_04, color='b', alpha=0.7)
# plt.plot(T_ab_IM_04, d_ab_IM_04, color='b', alpha=0.7)

# plt.plot(T_bl_IM_02, d_bl_IM_02, color='r', alpha=0.7)
# plt.plot(T_ab_IM_02, d_ab_IM_02, color='r', alpha=0.7)

# plt.plot(T_bl_GS_O2, d_bl_GS_O2, color='y', alpha=0.7)
# plt.plot(T_ab_GS_O2, d_ab_GS_O2, color='y', alpha=0.7)

# plt.plot(T_bl_pic, d_bl_pic, color='g', alpha=0.7)
# plt.plot(T_ab_pic, d_ab_pic, color='g', alpha=0.7)

# plt.plot(T_bl_JC_O2, d_bl_JC_O2, color='m', alpha=0.7)
# plt.plot(T_ab_JC_O2, d_ab_JC_O2, color='m', alpha=0.7)
bl_list = [T_bl_IM_04, T_bl_IM_02, T_bl_GS_O2, T_bl_pic, T_bl_JC_O2]
ab_list = [T_ab_IM_04, T_ab_IM_02, T_ab_GS_O2, T_ab_pic, T_ab_JC_O2]
d_bl_list = [d_bl_IM_04, d_bl_IM_02, d_bl_GS_O2, d_bl_pic, d_bl_JC_O2]
d_ab_list = [d_ab_IM_04, d_ab_IM_02, d_ab_GS_O2, d_ab_pic, d_ab_JC_O2]
label_list = ['IM $O(h^4)$', 'IM $O(h^2)$', 'GS $O(h^2)$', 'pictoral', 'JC $O(h^2)$']
colour_list = ['b','r','y','g','m']

fig = plt.figure(1, [7,5])
(ax1, ax2) = fig.subplots(1,2)

ax1.set_xlabel('Temperature (293K)')
ax1.set_ylabel('Difference from 500 updates ago')

def save_data(bl_list, d_bl_list, ab_list, d_ab_list, label, file_name, labels=['T_bl','d_bl','T_ab', 'd_ab']):
    M = [bl_list, d_bl_list, ab_list, d_ab_list]
    # print(M)
    M = np.matrix(np.array([bl_list, d_bl_list, ab_list, d_ab_list])).transpose()
    M = pd.DataFrame(M, index=None)
    M.columns = labels
    M.to_csv(f'Method_testing/{file_name}_{label}')

file_labels = ['IMO4', 'IMO2', 'GS', 'pictoral', 'JC']
for x1, y1, x2, y2, label, colour,file_label in zip(bl_list, d_bl_list, ab_list, d_ab_list, label_list, colour_list,file_labels):
    plot_the_temp(x1, y1, x2, y2, label, colour, 0.7, ax1)
    save_data(x1, y1, x2, y2, file_label, 'starting_20_30')

# save_data(bl_list, d_bl_list, ab_list, d_ab_list, label_list, 'starting_20_30')
# plt.legend()
# plt.show()

print('saved')

#        

ambient = 27.1
top = np.array([ambient for i in range(14*6)])
bottom = top.copy()
left = np.array([ambient for i in range(6)])
right = left.copy()

T_0_ab, T_k_ab = initialise_boundary(top, left, right, bottom, ambient)

ambient = 26.9
top = np.array([ambient for i in range(14*6)])
bottom = top.copy()
left = np.array([ambient for i in range(6)])
right = left.copy()

T_0_bl, T_k_bl = initialise_boundary(top, left, right, bottom, ambient)

#         test GS FDO2 above, below, at eq
from MatrixFunctions import update_T_k_GS_matrix, matrix_maker, find_D_Tl_Tr, gaus_jordan


d_ab_GS_O2 = []
d_bl_GS_O2 = []
T_ab_GS_O2 = []
T_bl_GS_O2 = []

T_0_GS_O2_ab = T_0_ab.copy()
T_k_GS_O2_ab = T_k_ab.copy()

T_0_GS_O2_bl = T_0_bl.copy()
T_k_GS_O2_bl = T_k_bl.copy()

M = matrix_maker(T_k_GS_O2_ab)
D, Tl, Tu = find_D_Tl_Tr(M)
M = D+Tl
M_inv = gaus_jordan(M)


T_ab = 27.1
T_bl = 26.9
t_before = time.time()
for n in range(20):
    for m in range(500):
        T_0_GS_O2_ab, T_k_GS_O2_ab = update_boundary_convection(T_k_GS_O2_ab, K, h=h), update_T_k_GS_matrix(T_k_GS_O2_ab, T_0_GS_O2_ab, M_inv, Tu, q=q,h=h)
        T_0_GS_O2_bl, T_k_GS_O2_bl = update_boundary_convection(T_k_GS_O2_bl, K, h=h), update_T_k_GS_matrix(T_k_GS_O2_bl, T_0_GS_O2_bl, M_inv, Tu, q=q,h=h)
    T_ab_next = average_temp(T_k_GS_O2_ab)
    T_bl_next = average_temp(T_k_GS_O2_bl)
    
    T_ab_GS_O2.append(T_ab_next)
    T_bl_GS_O2.append(T_bl_next)
    
    d_ab_GS_O2.append(T_ab-T_ab_next)
    d_bl_GS_O2.append(T_bl-T_bl_next)
    
    T_ab = T_ab_next
    T_bl = T_bl_next
t_after = time.time()
timeings.append(t_after-t_before)

print(f'\n\nDone w/ GSO2')
print(f'time taken: {timeings[len(timeings)-2]}')
#         test GS FDO4 above, below, at eq
from MatrixFunctions import matrix_maker_2


d_ab_GS_O4 = []
d_bl_GS_O4 = []
T_ab_GS_O4 = []
T_bl_GS_O4 = []

T_0_GS_O4_ab = T_0_ab.copy()
T_k_GS_O4_ab = T_k_ab.copy()

T_0_GS_O4_bl = T_0_bl.copy()
T_k_GS_O4_bl = T_k_bl.copy()

M = matrix_maker_2(T_k_GS_O2_ab)
D, Tl, Tu = find_D_Tl_Tr(M)
M = D+Tl
M_inv = gaus_jordan(M)


T_ab = 27.1
T_bl = 26.9
t_before = time.time()
for n in range(20):
    for m in range(500):
        T_0_GS_O4_ab, T_k_GS_O4_ab = update_boundary_convection(T_k_GS_O4_ab, K, h=h), update_T_k_GS_matrix(T_k_GS_O4_ab, T_0_GS_O4_ab, M_inv, Tu, q=q,h=h)
        T_0_GS_O4_bl, T_k_GS_O4_bl = update_boundary_convection(T_k_GS_O4_bl, K, h=h), update_T_k_GS_matrix(T_k_GS_O4_bl, T_0_GS_O4_bl, M_inv, Tu, q=q,h=h)
    T_ab_next = average_temp(T_k_GS_O4_ab)
    T_bl_next = average_temp(T_k_GS_O4_bl)
    
    T_ab_GS_O4.append(T_ab_next)
    T_bl_GS_O4.append(T_bl_next)
    
    d_ab_GS_O4.append(T_ab-T_ab_next)
    d_bl_GS_O4.append(T_bl-T_bl_next)
    
    T_ab = T_ab_next
    T_bl = T_bl_next
t_after = time.time()
timeings.append(t_after-t_before)
print(f'\n\nDone w/ GSO4')
print(f'time taken: {timeings[len(timeings)-2]}')
#         test pictoral above, below, at eq
from Updates import update_with_source

d_ab_pic = []
d_bl_pic = []
T_ab_pic = []
T_bl_pic = []

T_0_pic_ab = T_0_ab.copy()
T_k_pic_ab = T_k_ab.copy()

T_0_pic_bl = T_0_bl.copy()
T_k_pic_bl = T_k_bl.copy()

T_ab = 27.1
T_bl = 26.9
t_before = time.time()
for n in range(20):
    for m in range(500):
        T_0_pic_ab, T_k_pic_ab = update_boundary_convection(T_k_pic_ab, K, h=h), update_with_source(T_k_pic_ab, T_0_pic_ab, q=q,h=h)
        T_0_pic_bl, T_k_pic_bl = update_boundary_convection(T_k_pic_bl, K, h=h), update_with_source(T_k_pic_bl, T_0_pic_bl, q=q,h=h)
    T_ab_next = average_temp(T_k_pic_ab)
    T_bl_next = average_temp(T_k_pic_bl)
    
    T_ab_pic.append(T_ab_next)
    T_bl_pic.append(T_bl_next)
    
    d_ab_pic.append(T_ab-T_ab_next)
    d_bl_pic.append(T_bl-T_bl_next)
    
    T_ab = T_ab_next
    T_bl = T_bl_next


t_after = time.time()
timeings.append(t_after-t_before)
print(f'\n\nDone w/ pictoral')
print(f'time taken: {timeings[len(timeings)-2]}')
#         test jacobi FDO2 above, below, at eq
from MatrixFunctions import update_T_k_jacobi_matrix


d_ab_JC_O2 = []
d_bl_JC_O2 = []
T_ab_JC_O2 = []
T_bl_JC_O2 = []

T_0_JC_O2_ab = T_0_ab.copy()
T_k_JC_O2_ab = T_k_ab.copy()

T_0_JC_O2_bl = T_0_bl.copy()
T_k_JC_O2_bl = T_k_bl.copy()

M = matrix_maker(T_k_GS_O2_ab)
D, Tl, Tu = find_D_Tl_Tr(M)
M = Tu+Tl


T_ab = 27.1
T_bl = 26.9
t_before = time.time()
for n in range(20):
    for m in range(500):
        T_0_JC_O2_ab, T_k_JC_O2_ab = update_boundary_convection(T_k_JC_O2_ab, K, h=h), update_T_k_jacobi_matrix(T_k_JC_O2_ab, T_0_JC_O2_ab, M ,q=q,h=h)
        T_0_JC_O2_bl, T_k_JC_O2_bl = update_boundary_convection(T_k_JC_O2_bl, K, h=h), update_T_k_jacobi_matrix(T_k_JC_O2_bl, T_0_JC_O2_bl, M, q=q,h=h)
    T_ab_next = average_temp(T_k_JC_O2_ab)
    T_bl_next = average_temp(T_k_JC_O2_bl)
    
    T_ab_JC_O2.append(T_ab_next)
    T_bl_JC_O2.append(T_bl_next)
    
    d_ab_JC_O2.append(T_ab-T_ab_next)
    d_bl_JC_O2.append(T_bl-T_bl_next)
    
    T_ab = T_ab_next
    T_bl = T_bl_next

t_after = time.time()
timeings.append(t_after-t_before)
print(f'\n\nDone w/ jacobi')
print(f'time taken: {timeings[len(timeings)-2]}')
#         test inverse matrix FDO2 above, below, at eq

from MatrixFunctions import update_T_matrix_method
import pandas as pd


d_ab_IM_02 = []
d_bl_IM_02 = []
T_ab_IM_02 = []
T_bl_IM_02 = []

T_0_IM_02_ab = T_0_ab.copy()
T_k_IM_02_ab = T_k_ab.copy()

T_0_IM_02_bl = T_0_bl.copy()
T_k_IM_02_bl = T_k_bl.copy()

M = np.matrix(pd.read_csv('M_inverses/M_inv_pr_s6', header=None))

T_ab = 27.1
T_bl = 26.9
t_before = time.time()
for n in range(20):
    for m in range(500):
        T_0_IM_02_ab, T_k_IM_02_ab = update_boundary_convection(T_k_IM_02_ab, K, h=h), update_T_matrix_method(T_0_IM_02_ab, M ,q=q,h=h)
        T_0_IM_02_bl, T_k_IM_02_bl = update_boundary_convection(T_k_IM_02_bl, K, h=h), update_T_matrix_method(T_0_IM_02_bl, M, q=q,h=h)
    T_ab_next = average_temp(T_k_IM_02_ab)
    T_bl_next = average_temp(T_k_IM_02_bl)
    
    T_ab_IM_02.append(T_ab_next)
    T_bl_IM_02.append(T_bl_next)
    
    d_ab_IM_02.append(T_ab-T_ab_next)
    d_bl_IM_02.append(T_bl-T_bl_next)
    
    T_ab = T_ab_next
    T_bl = T_bl_next


t_after = time.time()
timeings.append(t_after-t_before)
print(f'\n\nDone w/ Inverse O2')
print(f'time taken: {timeings[len(timeings)-2]}')
#         test inverse matrix FDO4 above, below, at eq

d_ab_IM_04 = []
d_bl_IM_04 = []
T_ab_IM_04 = []
T_bl_IM_04 = []

T_0_IM_04_ab = T_0_ab.copy()
T_k_IM_04_ab = T_k_ab.copy()

T_0_IM_04_bl = T_0_bl.copy()
T_k_IM_04_bl = T_k_bl.copy()

M = np.matrix(pd.read_csv('M_inverses/M_inv_pr_s6_FDO4', header=None))

T_ab = 27.1
T_bl = 26.9
t_before = time.time()
for n in range(20):
    for m in range(500):
        T_0_IM_04_ab, T_k_IM_04_ab = update_boundary_convection(T_k_IM_04_ab, K, h=h), update_T_matrix_method(T_0_IM_04_ab, M ,q=q,h=h)
        T_0_IM_04_bl, T_k_IM_04_bl = update_boundary_convection(T_k_IM_04_bl, K, h=h), update_T_matrix_method(T_0_IM_04_bl, M, q=q,h=h)
    T_ab_next = average_temp(T_k_IM_04_ab)
    T_bl_next = average_temp(T_k_IM_04_bl)
    
    T_ab_IM_04.append(T_ab_next)
    T_bl_IM_04.append(T_bl_next)
    
    d_ab_IM_04.append(T_ab-T_ab_next)
    d_bl_IM_04.append(T_bl-T_bl_next)
    
    T_ab = T_ab_next
    T_bl = T_bl_next


t_after = time.time()
timeings.append(t_after-t_before)
print(f'\n\nDone w/ Inverse O4')
print(f'time taken: {timeings[len(timeings)-2]}')

#        
bl_list = [T_bl_IM_04, T_bl_IM_02, T_bl_GS_O2, T_bl_pic, T_bl_JC_O2]
ab_list = [T_ab_IM_04, T_ab_IM_02, T_ab_GS_O2, T_ab_pic, T_ab_JC_O2]
d_bl_list = [d_bl_IM_04, d_bl_IM_02, d_bl_GS_O2, d_bl_pic, d_bl_JC_O2]
d_ab_list = [d_ab_IM_04, d_ab_IM_02, d_ab_GS_O2, d_ab_pic, d_ab_JC_O2]
label_list = ['IM $O(h^4)$', 'IM $O(h^2)$', 'GS $O(h^2)$', 'pictoral', 'JC $O(h^2)$']
colour_list = ['b','r','y','g','m']


ax2.set_xlabel('Temperature (293K)')
# ax2.set_ylabel('Difference from 500 updates ago')

for x1, y1, x2, y2, label, colour,file_label in zip(bl_list, d_bl_list, ab_list, d_ab_list, label_list, colour_list,file_labels):
    plot_the_temp(x1, y1, x2, y2, label, colour, 0.7, ax2, labeling=False, start=0)
    save_data(x1, y1, x2, y2, file_label, 'starting_26p91_26p92')
ax1.set_title('a)')
ax2.set_title('b)')
ax1.legend()
plt.show()
# save_data(bl_list, d_bl_list, ab_list, d_ab_list, label_list, 'starting_26p91_26p92')

#%%