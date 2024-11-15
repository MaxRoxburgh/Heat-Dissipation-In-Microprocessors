import matplotlib.pyplot as plt
import numpy as np


#%%  plt.show()  Fig 1
from MatrixFunctions import matrix_maker_2, gaus_jordan

fig, (ax1, ax2) = plt.subplots(1,2)

a = np.ones((9,9), dtype='float')
b = matrix_maker_2(a)
b_inv = gaus_jordan(b)

cmap = 'jet_r' 

im1 = ax1.imshow(b, cmap)
im2 = ax2.imshow(b_inv, cmap)
ax1.set_title('a)')
ax2.set_title('b)')

ax1.axes.xaxis.set_visible(False)
ax1.axes.yaxis.set_visible(False)
ax2.axes.xaxis.set_visible(False)
ax2.axes.yaxis.set_visible(False)

#%%  
plt.show() 
import pandas as pd
from Plotting import trim, combine, plot_3d_matrix, add_boundary
import matplotlib.pyplot as plt
import numpy as np
plt.style.use('seaborn-poster')

T_cs = np.matrix(pd.read_csv('pr_cs_data_FDO4/scaling_8/T_k_cs_above_s8', header=None))*293
T_pr = np.matrix(pd.read_csv('pr_cs_data_FDO4/scaling_8/T_k_pr_above_s8', header=None))*293
combined = np.matrix(combine(T_cs, T_pr, 0))

combined[combined==0] = np.nan
cmap = plt.cm.get_cmap('plasma')
cmap.set_bad('#ffffff')

fig, (ax1, ax2) = plt.subplots(1,2)
im = plt.imshow(combined, 'plasma')
ax1.axes.xaxis.set_visible(False)
ax1.axes.yaxis.set_visible(False)
ax2.axes.xaxis.set_visible(False)
ax2.axes.yaxis.set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.spines['left'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['left'].set_visible(False)

ax1.set_title('a)')
ax2.set_title('b)', pad=219)

cbar = ax1.figure.colorbar(im, orientation='horizontal', pad=0.3)
cbar.set_label('Temperature (K)')
print(f'The minimum temperature of the system:  {round(min(trim(T_cs).A1))}K')
print(f'The maximum temperature of the system:  {round(max(T_pr.A1))}K')

# plt.show()
# plt.style.use('bmh')
combined = np.matrix(combine(T_cs, T_pr, 6080))
combined = add_boundary(combined, 6080, 4)
plot_3d_matrix(add_boundary(combined,6080, 4), 50, 69, ax2, 121, fig)
cmap = plt.cm.get_cmap('plasma')
cmap.set_bad('#ffffff')
#%%  
plt.show() 
import pandas as pd
from Plotting import trim, combine, plot_3d_matrix, add_boundary, plot_whole_system
from Calculators import mirror_array
import matplotlib.pyplot as plt
import numpy as np
plt.style.use('seaborn-poster')

T_cs = np.matrix(pd.read_csv('pr_cs_data_FDO4/scaling_8/T_k_cs_above_s8', header=None))*293
T_pr = np.matrix(pd.read_csv('pr_cs_data_FDO4/scaling_8/T_k_pr_above_s8', header=None))*293
combined = np.matrix(combine(T_cs, T_pr, 0))

combined[combined==0] = np.nan
cmap = plt.cm.get_cmap('plasma')
cmap.set_bad('#ffffff')

fig, (ax1, ax2) = plt.subplots(1,2)
im = ax1.imshow(combined, 'plasma')
ax1.axes.xaxis.set_visible(False)
ax1.axes.yaxis.set_visible(False)
ax2.axes.xaxis.set_visible(False)
ax2.axes.yaxis.set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.spines['left'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['left'].set_visible(False)

ax1.set_title('a)', pad=215)
ax2.set_title('b)')

cbar = ax1.figure.colorbar(im, orientation='horizontal', ax=ax1, pad=0.3, ticks = [6097, 6100, 6103, 6106, 6109])
cbar.set_label('Temperature (K)')
print(f'The minimum temperature of the system:  {round(min(trim(T_cs).A1))}K')
print(f'The maximum temperature of the system:  {round(max(T_pr.A1))}K')


T_cs_1 = np.matrix(pd.read_csv('whole_system_FDO4/fin_length_test_forced/T_k_cs_above_s6_fin_length_30', header=None))*293
T_pr_1 = np.matrix(pd.read_csv('whole_system_FDO4/fin_length_test_forced/T_k_pr_above_s6_fin_length_30', header=None))*293
T_hs_1 = np.matrix(pd.read_csv('whole_system_FDO4/fin_length_test_forced/T_k_hs_above_s6_fin_length_30', header=None))*293
T_fn_1 = np.array([np.matrix(pd.read_csv(f'whole_system_FDO4/fin_length_test_forced/T_k_fn_above_s6_fin_length_30_{i}', header=None))*293 for i in range(6)])

comined_system = plot_whole_system(T_cs_1, T_pr_1, mirror_array(T_fn_1, True), T_hs_1, 0, 6, False)

comined_system[comined_system==0] = np.nan
cmap = plt.cm.get_cmap('plasma')
cmap.set_bad('#ffffff')

im2 = plt.imshow(comined_system, 'plasma')
cbar2 = ax2.figure.colorbar(im2, orientation='horizontal', pad=0.1, ax=ax2)
cbar2.set_label('Temperature (K)')
# plt.show()
# plt.style.use('bmh')
# combined = np.matrix(combine(T_cs, T_pr, 6080))
# combined = add_boundary(combined, 6080, 4)
# plot_3d_matrix(add_boundary(combined,6080, 4), 50, 69, ax2, 121, fig)
# cmap = plt.cm.get_cmap('plasma')
# cmap.set_bad('#ffffff')
#%%  
plt.show() 
import pandas as pd
from Plotting import trim, combine, plot_3d_matrix, add_boundary
import matplotlib.pyplot as plt
import numpy as np
plt.style.use('bmh')
scalings = [3,4,5,6,7,8,9,10,11]
T_above = []
T_below = []
T_above_FDO4 = []
T_below_FDO4 = []

for scaling in scalings:
    try:
        data = pd.read_csv(f'pr_data/scaling_{scaling}/data_s{scaling}')
        T_below_array, T_above_array = np.array(data['T_bellow']), np.array(data['T_above'])
        T_below.append(T_below_array[-1])
        T_above.append(T_above_array[-1])
        data = pd.read_csv(f'pr_data_FDO4/scaling_{scaling}/data_s{scaling}')
        T_below_array, T_above_array = np.array(data['T_bellow']), np.array(data['T_above'])
        T_below_FDO4.append(T_below_array[-1])
        T_above_FDO4.append(T_above_array[-1])
    except:
        next


T_above = np.array(T_above)
T_below = np.array(T_below)
T_above_FDO4 = np.array(T_above_FDO4)
T_below_FDO4 = np.array(T_below_FDO4)

T_FDO2 = (T_above + T_below)*293/2
err_FDO2 = np.array([0.5*293 for i in T_FDO2])
T_FDO4 = (T_above_FDO4 + T_below_FDO4)*293/2
err_FDO4 = np.array([0.5*293 for i in T_FDO4])

plt.figure(1, [7,5])
plt.errorbar(scalings, T_FDO2, err_FDO2, capsize=8, fmt='o--', alpha=1, label='O($h^2$)')
plt.errorbar([scalings[i] for i in range(len(T_FDO4))], T_FDO4, err_FDO4, capsize=8, fmt='o--', alpha=0.5, label='O($h^4$)')
plt.legend()
plt.xlabel('Points per mm')
plt.ylabel('Temperature K')
plt.show()
#%%  
plt.show() 

plt.style.use('bmh')
scalings = [3,4,5,6,7,8,9,10,11]
T_above_2 = []
T_below_2 = []
T_above_FDO4_2 = []
T_below_FDO4_2 = []

for scaling in scalings:
    try:
        data = pd.read_csv(f'pr_cs_data/scaling_{scaling}/data_s{scaling}')
        T_below_array_2, T_above_array_2 = np.array(data['T_bellow']), np.array(data['T_above'])
        T_below_2.append(T_below_array_2[-1])
        T_above_2.append(T_above_array_2[-1])
        data = pd.read_csv(f'pr_cs_data_FDO4/scaling_{scaling}/data_s{scaling}')
        T_below_array_2, T_above_array_2 = np.array(data['T_bellow']), np.array(data['T_above'])
        T_below_FDO4_2.append(T_below_array_2[-1])
        T_above_FDO4_2.append(T_above_array_2[-1])
    except:
        next


T_above_2 = np.array(T_above_2)
T_below_2 = np.array(T_below_2)
T_above_FDO4_2 = np.array(T_above_FDO4_2)
T_below_FDO4_2 = np.array(T_below_FDO4_2)

T_FDO2_2 = (T_above_2 + T_below_2)*293/2
err_FDO2_2 = np.array([0.5*293 for i in T_FDO2_2])
T_FDO4_2 = (T_above_FDO4_2 + T_below_FDO4_2)*293/2
err_FDO4_2 = np.array([0.5*293 for i in T_FDO4_2])

plt.figure(1, [7,5])
scalings_2 = [scalings[i] for i in range(len(T_above_2))]
plt.errorbar(scalings_2, T_FDO2_2, err_FDO2_2, capsize=8, fmt='o--', alpha=1, label='O($h^2$)')
plt.errorbar([scalings[i] for i in range(len(T_FDO4_2))], T_FDO4_2, err_FDO4_2, capsize=8, fmt='o--', alpha=0.5, label='O($h^4$)')
plt.legend()
plt.xlabel('Points per mm')
plt.ylabel('Temperature K')
plt.show()

#%% 
plt.show() 

fig = plt.figure(1, [9,5])
ax1, ax2 = fig.subplots(1,2)

ax1.set_title('a)')
ax2.set_title('b)')


ax1.errorbar(scalings_2, T_FDO2_2, err_FDO2_2, capsize=8, fmt='o--', alpha=1, label='IMO($h^2$)')
ax1.errorbar([scalings[i] for i in range(len(T_FDO4_2))], T_FDO4_2, err_FDO4_2, capsize=8, fmt='o--', alpha=0.5, label='IMO($h^4$)')
ax1.legend()
ax1.set_xlabel('Points per mm')
ax1.set_ylabel('Temperature K')

ax2.errorbar(scalings, T_FDO2, err_FDO2, capsize=8, fmt='o--', alpha=1, label='IMO($h^2$)')
ax2.errorbar([scalings[i] for i in range(len(T_FDO4))], T_FDO4, err_FDO4, capsize=8, fmt='o--', alpha=0.5, label='IMO($h^4$)')
ax2.legend()
ax2.set_xlabel('Points per mm')
# ax2.set_ylabel('Temperature K')
plt.show()



fig = plt.figure(1, [7,5])
fig, ax = plt.subplots(1,1)
ax.errorbar(scalings, T_FDO2*(T_FDO2_2[0]/T_FDO2[0]), err_FDO2*(T_FDO2_2[0]/T_FDO2[0]), capsize=8, fmt='o--', alpha=1, label='scaled down prcessor with ceramic')
ax.errorbar(scalings_2, T_FDO2_2, err_FDO2_2, capsize=8, fmt='o--', alpha=1, label='just processor')
# ax.errorbar(scalings,  np.array([T_FDO2[i]/T_FDO2[0]*T_FDO2_2[0] for i in range(len(T_FDO2))]))
# ax.axes.yaxis.set_visible(False)
ax.set_ylabel('Temperature [K]')
ax.legend()
# ax.grid()
ax.set_xlabel('Points per mm')
plt.show()

print(f'There will be at least {round(100*abs(1-T_FDO2[-1]/T_FDO2[3]))}% difference if using 6 ppm')


#%%  initial boundary interface 
plt.show()  
plt.style.use('bmh')
def new_imaginary_point_2(points_in_order, k1, k2):
    # the points in order will be the 3rd to last in the left array and the 3rd of the right array
    return (k1*points_in_order[0] + k2*points_in_order[1])/(k1+k2)

def update(T_0, T_k):
    n = len(T_0)
    T_kPlus1 = T_0.copy()
    for index in range(1,n-1):
        T_kPlus1[index] += 1/2*(T_k[index-1] + T_k[index+1])
    
    return T_kPlus1

def lin(x, a, b):
    return a*x+b

from scipy.optimize import curve_fit

k1 = 50
k2 = 100

def example(n):
    # x = np.linspace(0+1/((n-1)*2),1,(n-1)*2)
    # h = x[1]-x[0]


    T_k_Material1 = np.array([0.5 for i in range(n+2)])
    T_k_Material2 = np.array([0.5 for i in range(n+2)])


    # we may initialise the starting values of the imaginary points like this as the entire system is in equilibrium

    new_point = new_imaginary_point_2([T_k_Material1[n-3], T_k_Material2[2]], k1, k2)


    T_0_Material1 = np.append(np.append([0.5], [0 for  i in range(n)]), [new_point])
    T_0_Material2 = np.append(np.append([new_point], [0 for  i in range(n)]), [1])

    # where the 0's in the seperate array are is where the new values will be put
    # notice that the array is one longer as we are generating an imaginary point for each of the materials


    for i in range(500*(n-1)):
        # start by updating material 1
        T_k_Material1_update = update(T_0_Material1, T_k_Material1)

        # update material 2
        T_k_Material2_update = update(T_0_Material2, T_k_Material2)

        # with the updated arrays, update the T_0 with the new BC
        new_point = new_point = new_imaginary_point_2([T_k_Material1_update[n-3], T_k_Material2_update[2]], k1, k2)
        T_0_Material1 = np.append(np.append([0.5], [0 for  i in range(n-2)]), [new_point])
        T_0_Material2 = np.append(np.append([new_point], [0 for  i in range(n-2)]), [1])

        # update the origonal T_k's to the new values of the updated arrrays
        T_k_Material1 = np.array([i for i in T_k_Material1_update])
        T_k_Material2 = np.array([i for i in T_k_Material2_update])

    T = np.append(T_k_Material1[1:-1],T_k_Material2[1:-1])
    print(len(T))
    x = [(i+1)/(len(T)+1) for i in range(len(T))]
    plt.plot(x, T*100, 'o--', ms=1, label=f'{len(T)} points per mm')
    print(T_k_Material1[-2])


plt.figure(9, [7,4])
# example(4)
example(5)
# example(6)
example(9)
# example(17)
example(72)

plt.legend()
plt.xlabel('Distance [mm]')
plt.ylabel('Temp [K]')


#%%  timings for different methods plot
plt.show()

timings = np.array([29.43551898, 29.37344545, 48.1245181 , 27.00512168, 26.28961517, 26.09180913])
label_list = np.array(['GS $O(h^2)$', 'GS $O(h^4)$','pictoral', 'JC $O(h^2)$','IM $O(h^2)$', 'IM $O(h^4)$'])
# label_list = np.array(np.fliplr(np.matrix(label_list)).A1)
sort = True
while sort:
    sort = False
    for i in range(len(timings)-1):
        if timings[i+1] > timings[i]:
            sort = True
            label_list[i], label_list[i+1] = label_list[i+1], label_list[i]
            timings[i], timings[i+1] = timings[i+1], timings[i]

# print(timings, label_list)
plt.bar(label_list, timings)
plt.ylabel('time taken for 10k itereations (s)')


#%%  
plt.show() 
from scipy.optimize import curve_fit
from Calculators import average_temp
import pandas as pd
plt.style.use('bmh')

def lin(x,a,b):
    return a*x + b
def one_over_x_3_4(x, a, b):
    return (a/(x))**(3/4)/293 + b

def one_over_(x, a, b):
    return (a/(x))**(7/4)/293 + b
def err_finder(l, above, path='whole_system_FDO4/fin_length_test', name='data_s6_fin_length_'):
    df = pd.read_csv(f'{path}/{name}{l}')
    if above:
        diff = np.array(df['diff_above'])
        return abs(diff[-1])
    else:
        diff = np.array(df['diff_bellow'])
        return abs(diff[-1])

T_points = []

def plot_the_temp(path, n, starting,ax ,col= 'r'):
    df = pd.read_csv(f'{path}/data_s6_nfins_{n}_2')
    
    T_a, T_b = df['T_above'][starting:], df['T_bellow'][starting:]
    D_a, D_b = df['diff_above'][starting:], df['diff_bellow'][starting:]

    ax.plot(T_a*293, D_a, color=col, label=f'{n} fins')
    ax.plot(T_b*293, D_b, color=col)
    
    (a,b), cov = curve_fit(lin, np.append(T_a,T_b), np.append(D_a,D_b))
    ax.plot(np.append(T_a,T_b)*293, [lin(i,a,b) for i in np.append(T_a,T_b)], alpha=0.2, color=col)
    ax.plot([-b*293/a],[0], 'o', color=col)
    print(f'Converge point for {n} number of fins is {-b/a}')
    n=len(T_a)
    print(f'difference from below to mid = {(-b/a - T_b[n-1])*293}')
    print(f'difference from above to mid = {(-b/a - T_a[n-1])*293}\n')
    return -b/a
    
x = np.array([8,10,12,14,16]) 
col = ['r','g','b','y','m']
fig = plt.figure(99,[9,5])
ax1, ax2 = fig.subplots(1,2)
for i, c in zip(x, col):
    T_points.append(plot_the_temp('whole_system_FDO4/n_fins_test', i, 8, ax1,c))
ax1.set_ylim([-0.001,0.0016])
ax1.legend()
ax1.set_xlabel('Temperature [K]')
ax1.set_ylabel('difference per 1k runs')
ax2.set_xlabel('no. of fins')
ax2.set_ylabel('Temperature [K]')
ax1.set_title('a)')
ax2.set_title('b)')

# for xx, y in zip(x, T_points):
ax2.plot(x, np.array(T_points)*293, 'o', label='system average Temp', color='gray')

fit, cov = curve_fit(one_over_x_3_4, x,np.array(T_points)*293)
xx = np.linspace(min(x),max(x), 1000)
# ax2.plot(xx, [one_over_x_3_4(i, fit[0], fit[1]) for i in xx], alpha=0.8, color='r', label='fitted f($x^{-3/4}$)')

def plot_pr_temp(path, n, ax ,col= 'r', above=True):
    if above:
        df = np.matrix(pd.read_csv(f'{path}/T_k_pr_above_s6_nfin_{n}_final', header=None))
    else:
        df = np.matrix(pd.read_csv(f'{path}/T_k_pr_bellow_s6_nfin_{n}_final', header=None))
    T = average_temp(df)
    # ax.plot(n,T*293,'o',color=col)
    return T


T_pr_above= []    
T_pr_below=[]
above_or_below=[False,False,False,False,False]
for n in x:
    T_pr_above.append(plot_pr_temp('whole_system_FDO4/n_fins_test', n, ax2, 'y', False))
for n in x:
    T_pr_below.append(plot_pr_temp('whole_system_FDO4/n_fins_test', n, ax2, 'y', True))
    
T_above_err = []
T_bellow_err = []

for l in x:
    T_above_err.append(err_finder(l,True, path='whole_system_FDO4/n_fins_test', name='data_s6_nfins_'))
for l in x:
    T_bellow_err.append(err_finder(l,False, path='whole_system_FDO4/n_fins_test', name='data_s6_nfins_'))

# (a,b), cov = curve_fit(one_over_x_3_4, x, T_pr)
# (a,b), cov = curve_fit(one_over_, x, T_pr)
# ax2.plot(xx, [one_over_(i, a, b)*293 for i in xx], alpha=0.8, color='r', label='fitted f($x^{-3/4}$)')

T_pr = np.array( [(T_pr_above[i]*1/T_above_err[i] + T_pr_below[i]*1/T_bellow_err[i])/((1/T_above_err[i] + 1/T_bellow_err[i])) for i in range(len(x))])
T_pr_err = [abs(T_pr_above[i] - T_pr_below[i])*293/2 for i in range(len(x))]

plt.errorbar(x, np.array(T_pr)*293, yerr=T_pr_err, fmt='o', capsize=5, color='black', label='Processor average temp')

(a,b), cov = curve_fit(one_over_x_3_4, x, T_pr, sigma=T_pr_err)
ax2.plot(xx, [one_over_x_3_4(i, a, b)*293 for i in xx], alpha=0.8, color='r', label='fitted f($x^{-3/4}$)')
ax2.legend()

print(f'The distribution apears to plataue and extrapolating the fit, the data sugests that with an infinite array of processors, the temperature could at best get to {round(b*293-273)} degrees C')

#%%  
plt.show() 
from scipy.optimize import curve_fit

T_points_fin = []

def plot_the_temp_2(path, l, starting,ax ,col= 'r'):
    df = pd.read_csv(f'{path}/data_s6_fin_length_{l}')
    if l == 42 or l ==30:
        starting = 18
    T_a, T_b = df['T_above'][starting:], df['T_bellow'][starting:]
    D_a, D_b = df['diff_above'][starting:], df['diff_bellow'][starting:]
    ax.plot(T_a*293, D_a, color=col, label=f'{l} length fins')
    ax.plot(T_b*293, D_b, color=col)
    
    (a,b), cov = curve_fit(lin, np.append(T_a,T_b), np.append(D_a,D_b))
    ax.plot(np.append(T_a,T_b)*293, [lin(i,a,b) for i in np.append(T_a,T_b)], alpha=0.2, color=col)
    ax.plot([-b*293/a],[0], 'o', color=col)
    print(f'Converge point for {l} length fin is {-b/a}')
    # n=len(T_b)
    print(f'difference from below to mid = {(-b/a - np.array(T_b)[-1])*293}')
    print(f'difference from above to mid = {(-b/a - np.array(T_a)[-1])*293}\n')
    return -b/a
    
x_fin = np.array([30,32, 34, 36, 38, 40, 42, 44]) 
col = ['r','g','b','y','m', 'r', 'g', 'b', 'y']
fig = plt.figure(2,[9,5])
ax = fig.subplots(1,1)
# T_points.append(plot_the_temp('whole_system_FDO4/n_fins_test', 8, 8, ax, 'm'))
for i, c in zip(x_fin, col):
    T_points_fin.append(plot_the_temp_2('whole_system_FDO4/fin_length_test', i, 2, ax,c))
xx = np.append([30], x)
#%%  
plt.show() 
# ax1.set_ylim([-0.001,0.0016])
# ax1.legend()
# ax1.set_xlabel('Temperature [K]')
# ax1.set_ylabel('difference per 1k runs')
# ax2.set_xlabel('no. of fins')
# ax2.set_ylabel('Temperature [K]')
# ax1.set_title('a)')
# ax2.set_title('b)')

# # for xx, y in zip(x, T_points):
# ax2.plot(x, np.array(T_points)*293, 'o', label='system average Temp', color='gray')

# fit, cov = curve_fit(one_over_x_3_4, x,np.array(T_points)*293)
# xx = np.linspace(min(x),max(x), 1000)
# ax2.plot(xx, [one_over_x_3_4(i, fit[0], fit[1]) for i in xx], alpha=0.8, color='r', label='fitted f($x^{-3/4}$)')

def plot_pr_temp(path, l, ax ,col= 'r', above=True):
    if above:
        df = np.matrix(pd.read_csv(f'{path}/T_k_pr_above_s6_fin_length_{l}', header=None))
    else:
        df = np.matrix(pd.read_csv(f'{path}/T_k_pr_bellow_s6_fin_length_{l}', header=None))
    T = average_temp(df)
    # ax.plot(n,T*293,'o',color=col)
    return T

def err_finder(l, above):
    path = 'whole_system_FDO4/fin_length_test'
    df = pd.read_csv(f'{path}/data_s6_fin_length_{l}')
    if above:
        diff = np.array(df['diff_above'])
        return abs(diff[-1])
    else:
        diff = np.array(df['diff_bellow'])
        return abs(diff[-1])
        
        

T_pr_above= []    
T_pr_below=[]
# above_or_below=[False,False,False,False,False]
for l in x_fin:
    T_pr_above.append(plot_pr_temp('whole_system_FDO4/fin_length_test', l, ax2, 'y', False))
for l in x_fin:
    T_pr_below.append(plot_pr_temp('whole_system_FDO4/fin_length_test', l, ax2, 'y', True))

T_above_err = []
T_bellow_err = []

for l in x_fin:
    T_above_err.append(err_finder(l,True))
for l in x_fin:
    T_bellow_err.append(err_finder(l,False))

T_pr_err_fin = [abs(T_pr_above[i] - T_pr_below[i])*293/2 for i in range(len(x_fin))]
T_pr_fin = np.array( [(T_pr_above[i]*1/T_above_err[i] + T_pr_below[i]*1/T_bellow_err[i])/((1/T_above_err[i] + 1/T_bellow_err[i])) for i in range(len(x_fin))])*293
plt.errorbar(x_fin, T_pr_fin, yerr=T_pr_err_fin, capsize=5, fmt='o',  markeredgewidth=2)
# plt.show()

# (a,b), cov = curve_fit(one_over_x_3_4, x, T_pr)
# (a,b), cov = curve_fit(one_over_, x, T_pr)
# ax2.plot(xx, [one_over_(i, a, b)*293 for i in xx], alpha=0.8, color='r', label='fitted f($x^{-3/4}$)')

# T_pr = [(T_pr_above[i] + T_pr_below[i])/2 for i in range(len(x))]
# T_pr_err = [abs(T_pr_above[i] - T_pr_below[i])*293/2 for i in range(len(x))]

# plt.errorbar(x, np.array(T_pr)*293, yerr=T_pr_err, fmt='o', capsize=10, color='black', label='Processor average temp')
xx_fin = np.linspace(min(x_fin),max(x_fin), 1000)
(a,b), cov = curve_fit(one_over_x_3_4, x_fin, T_pr_fin, sigma=T_pr_err_fin)
plt.plot(xx_fin, [one_over_x_3_4(i, a, b) for i in xx_fin], alpha=0.8, color='r', label='fitted f($x^{-3/4}$)')
# ax2.legend()

print(f'The distribution apears to plataue and extrapolating the fit, the data sugests that with an infinite array of processors, the temperature could at best get to {round(b*293-273)} degrees C')

#%%  
plt.show() 

fig = plt.figure(99, [9,5])

ax1, ax2 = fig.subplots(1,2, sharey=True)
(a,b), cov = curve_fit(one_over_x_3_4, x_fin, T_pr_fin, sigma=T_pr_err_fin)
ax2.errorbar(x_fin, np.array(T_pr_fin), yerr=T_pr_err_fin, capsize=5, fmt='o',  markeredgewidth=2, color='black', ms=5)
ax2.plot(xx_fin, [one_over_x_3_4(i, a, b) for i in xx_fin], alpha=0.8, color='r', label='fitted f($x^{-3/4}$)')
ax2.set_xlabel('Length of fins [mm]')
ax2.set_xticks(x_fin) 
ax1.set_xticks(x)
ax1.set_xlabel('No. of fins')
ax1.errorbar(x, np.array(T_pr)*293, yerr=T_pr_err, fmt='o', capsize=5, color='black', label='Processor average temp')
xx = np.linspace(min(x),max(x), 1000)
(a,b), cov = curve_fit(one_over_x_3_4, x, T_pr, sigma=T_pr_err)
ax1.plot(xx, [one_over_x_3_4(i, a, b)*293 for i in xx], alpha=0.8, color='r', label='fitted f($x^{-3/4}$)')
ax1.legend()

#%%  
plt.show() 

mm_added_fins = [0,132,132*2,132*3,132*4]
mm_added_length = [32,64,96,128,160,32*6,32*7]

fig = plt.figure(64, [9,5])
ax2, ax = fig.subplots(1,2, sharey=True)

ax.errorbar(mm_added_length, np.array(T_pr_fin)[1:], yerr=T_pr_err_fin[1:], capsize=5, fmt='o',  markeredgewidth=2, color='blue', ms=3, label='length of fins study')
ax.errorbar(mm_added_fins, np.array(T_pr)*293, yerr=np.array(T_pr_err), fmt='o', capsize=5, color='black', label='Number of fins study', ms=3,  markeredgewidth=2)
# plt.errorbar([0], )
mm_added_fins = [132,132*2,132*3,132*4]
mm_added_length = [0,32,64,96,128,160,32*6,32*7]
xxx = np.append(mm_added_fins, mm_added_length)
yyy = np.append(np.array(T_pr)[1:]*293, np.array(T_pr_fin))
yyerr = np.append(np.array(T_pr_err)[1:], T_pr_err_fin)
# plt.errorbar(xxx, yyy, yerr=yyerr, capsize=5, fmt='o',  markeredgewidth=2, color='black', ms=3)

Fit, cov = curve_fit(one_over_x_3_4, xxx+538, yyy, sigma=yyerr)

x_range = np.linspace(min(xxx), max(xxx), 1000)
ax.plot(x_range, [one_over_x_3_4(i+538, Fit[0], Fit[1]) for i in x_range], color='r')
ax.set_xlabel('Difference in mm')
ax2.set_ylabel('Temperature of processor [K]')
ax2.set_xlabel('length of fins [mm]')
ax2.errorbar(x_fin, np.array(T_pr_fin), yerr=T_pr_err_fin, capsize=5, fmt='o',  markeredgewidth=2, color='black', ms=5)
ax.legend()
ax.set_title('b)')
ax2.set_title('a)')

#%%  
plt.show() 


def plot_the_temp_3(path, l, starting, ax ,col= 'r'):
    df = pd.read_csv(f'{path}/data_s6_fin_length_{l}')
    # pr = pd.read_csv(f'{path}/T_k_pr')
    
    T_b, T_a = df['T_above'][starting:], df['T_bellow'][starting:]
    D_b, D_a = df['diff_above'][starting:], df['diff_bellow'][starting:]

    ax.plot(T_a*293, D_a, color=col, label=f'{l} length fins')
    ax.plot(T_b*293, D_b, color=col)
    
    (a,b), cov = curve_fit(lin, np.append(T_a,T_b), np.append(D_a,D_b))
    ax.plot(np.append(T_a,T_b)*293, [lin(i,a,b) for i in np.append(T_a,T_b)], alpha=0.2, color=col)
    ax.plot([-b*293/a],[0], 'o', color=col)
    print(f'Converge point for {l} fin length is {-b/a}')
    n=len(T_a)
    print(f'difference from below to mid = {(-b/a - np.array(T_b)[-1])*293}')
    print(np.array(T_b)[-1])
    print(f'difference from above to mid = {(-b/a - np.array(T_a)[-1])*293}\n')
    return -b/a

T_forced = []
lengths = np.array([30,32,34,36,38,40,42,44])
for l in lengths:
    T_forced.append(plot_the_temp_3('whole_system_FDO4/fin_length_test_forced', l, 10, plt))
plt.xlim([334,365])
plt.show()
#%%

mms = np.array([60+24*i for i in lengths])
plt.plot(mms, 293*np.array(T_forced), 'o--')
T_pr_forced = []
path = 'whole_system_FDO4/fin_length_test_forced'
lengths = np.array([30,32,34,36,38,40,42,44])
for l in lengths:
    T_pr_forced.append(average_temp(np.matrix(pd.read_csv(f'{path}/T_k_pr_bellow_s6_fin_length_{l}', header=None)))*293)

plt.plot(mms, T_pr_forced, 'o--')

def inv_lin(x, a, b):
    return a/x + b

(a1, b1), cov = curve_fit(inv_lin, mms, np.array(T_pr_forced))
plt.plot(mms, [inv_lin(i, a1, b1) for i in mms])

#%%
def plot_the_temp_5(path, l, starting, ax ,col= 'r'):
    df = pd.read_csv(f'{path}/data_s6_fin_length_{l}')
    # pr = pd.read_csv(f'{path}/T_k_pr')
    
    T_b, T_a = df['T_pr_above'][starting:]*293, df['T_pr_below'][starting:]*293
    D_b, D_a = df['diff_above'][starting:], df['diff_bellow'][starting:]

    ax.plot(T_a, D_a, color=col, label=f'{l} length fins')
    ax.plot(T_b, D_b, color=col)
    
    (a,b), cov = curve_fit(lin, np.append(T_a,T_b), np.append(D_a,D_b))
    ax.plot(np.append(T_a,T_b), [lin(i,a,b) for i in np.append(T_a,T_b)], alpha=0.2, color=col)
    ax.plot([-b/a],[0], 'o', color=col)
    print(f'Converge point for {l} fin length is {-b/a}')
    n=len(T_a)
    print(f'difference from below to mid = {(-b/a - np.array(T_b)[-1])}')
    print(np.array(T_b)[-1])
    print(f'difference from above to mid = {(-b/a - np.array(T_a)[-1])}\n')
    return -b/a

T_forced = []
lengths = np.array([22,24,30,34,38,42,44])
for l in lengths:
    T_forced.append(plot_the_temp_5('whole_system_FDO4/fin_length_test_forced', l, 10, plt))
# plt.xlim([334,365])
plt.show()
#%%
# def one_over_x_3_4(x, a, b):
#     return (a/(x))**(3/4) + b
mms = np.array([60+24*i for i in lengths])
(a1, b1), cov = curve_fit(inv_lin, mms, np.array(T_forced))
plt.plot(mms, [inv_lin(i, a1, b1) for i in mms], '-')
plt.plot(mms, T_forced, 'o--')
#
scaling_factor = 0.985
operational_temp = (80 + 273)*scaling_factor
print(round(1/((operational_temp-b1)/a1)), ' is the number of mm area for the system to be at 80 degrees')
n = np.array([i+25 for i in range(8)])
l = np.array([(round(1/((operational_temp-b1)/a1)) - 13)/(2*i) - 1 for i in n])
area = np.array([(l[i] + 7)*(2*n[i]-1) for i in range(len(n))])
print(n)
print(l)

#%%  plt.show() 
x = 540
n = [i+8 for i in range(13)]
l = [(x-13)/(2*i) - 1 for i in n]

#%%
starting = 8
path = 'whole_system_FDO4/fin_length_test_forced'
df = pd.read_csv(f'{path}/data_s10_fin_length_30')
T_b, T_a = df['T_pr_above'][starting:]*293, df['T_pr_below'][starting:]*293
D_b, D_a = df['diff_above'][starting:], df['diff_bellow'][starting:]

plt.plot(T_b, D_b, color='r')
plt.plot(T_a, D_a, color='r')

(A,B), cov = curve_fit(lin, np.append(T_b,T_a), np.append(D_b,D_a))
xxxx = np.linspace(370,400,400)
plt.plot(xxxx, [lin(i, A, B) for i in xxxx], alpha=0.1)
print(-B/A)

starting = 6
path = 'whole_system_FDO4/fin_length_test_forced'
df = pd.read_csv(f'{path}/data_s6_fin_length_30')
T_b, T_a = df['T_pr_above'][starting:]*293, df['T_pr_below'][starting:]*293
D_b, D_a = df['diff_above'][starting:], df['diff_bellow'][starting:]

plt.plot(T_b, D_b, color='r')
plt.plot(T_a, D_a, color='r')

(A1,B1), cov = curve_fit(lin, np.append(T_b,T_a), np.append(D_b,D_a))
xxxx = np.linspace(370,400,400)
plt.plot(xxxx, [lin(i, A1, B1) for i in xxxx], alpha=0.1)
print(-B1/A1)

print('therefore the percentage difference in scaling is', round((B*A1/(B1*A)-1)*100, 1), '%')

scaling_factor = 1/(B*A1/(B1*A))
target_temp_C = 80
operational_temp = (target_temp_C + 273)*scaling_factor
print(round(1/((operational_temp-b1)/a1)), f' is the number of mm area for the system to be at {target_temp_C} degrees')
n = np.array([i+25 for i in range(50)])
l = np.array([(round(1/((operational_temp-b1)/a1)) - 13)/(2*i) - 1 for i in n])
area = np.array([(l[i] + 7)*(2*n[i]-1) for i in range(len(n))])

target_temp_C = 85
operational_temp = (target_temp_C + 273)*scaling_factor
print(round(1/((operational_temp-b1)/a1)), f' is the number of mm area for the system to be at {target_temp_C} degrees')
n = np.array([i+25 for i in range(50)])
l = np.array([(round(1/((operational_temp-b1)/a1)) - 13)/(2*i) - 1 for i in n])
area = np.array([(l[i] + 7)*(2*n[i]-1) for i in range(len(n))])
# print(n)
# print(l)
# n = np.array([i for i in range(25)])
# l = np.array([(round(1/((operational_temp-b1)/a1)) - 13)/(2*i) - 1 for i in n])
# area = np.array([(l[i] + 7)*(2*n[i]-1) for i in range(len(n))])

#%%  
plt.show() 

def plot_the_temp_4(path, l, starting, ax ,col= 'r'):
    if l ==2:
        return
    df = pd.read_csv(f'{path}/data_s6_spacing_{l}')
    
    T_a, T_b = df['T_above'][starting:], df['T_bellow'][starting:]
    D_a, D_b = df['diff_above'][starting:], df['diff_bellow'][starting:]

    ax.plot(T_a*293, D_a, color=col, label=f'{l} length fins')
    ax.plot(T_b*293, D_b, color=col)
    
    (a,b), cov = curve_fit(lin, np.append(T_a,T_b), np.append(D_a,D_b))
    ax.plot(np.append(T_a,T_b)*293, [lin(i,a,b) for i in np.append(T_a,T_b)], alpha=0.2, color=col)
    ax.plot([-b*293/a],[0], 'o', color=col)
    print(f'Converge point for {l} fin length is {-b/a}')
    n=len(T_a)
    print(f'difference from below to mid = {(-b/a - np.array(T_b)[-1])*293}')
    print(np.array(T_b)[-1])
    print(f'difference from above to mid = {(-b/a - np.array(T_a)[-1])*293}\n')
    return -b/a


spacings = np.array([1,2,3,4])
for l in spacings:
    plot_the_temp_4('whole_system_FDO4/spacing_test', l, 8, plt)
# plt.xlim([334,355])
plt.show()



#%%
def plot_the_temp_5(path, n, starting, ax ,col= 'r'):
    df = pd.read_csv(f'{path}/data_s6_nfins_{n}')
    # pr = pd.read_csv(f'{path}/T_k_pr')
    
    T_b, T_a = df['T_pr_above'][starting:]*293, df['T_pr_below'][starting:]*293
    D_b, D_a = df['diff_above'][starting:], df['diff_bellow'][starting:]

    ax.plot(T_a, D_a, color=col, label=f'{l} length fins')
    ax.plot(T_b, D_b, color=col)
    
    (a,b), cov = curve_fit(lin, np.append(T_a,T_b), np.append(D_a,D_b))
    ax.plot(np.append(T_a,T_b), [lin(i,a,b) for i in np.append(T_a,T_b)], alpha=0.2, color=col)
    ax.plot([-b/a],[0], 'o', color=col)
    print(f'Converge point for {l} fin length is {-b/a}')
    n=len(T_a)
    print(f'difference from below to mid = {(-b/a - np.array(T_b)[-1])}')
    print(np.array(T_b)[-1])
    print(f'difference from above to mid = {(-b/a - np.array(T_a)[-1])}\n')
    return -b/a
def inv_lin(x,a,b):
    return a/x + b
# def inv_lin(x,a,b):
    # return a*x**(-3/4) + b

T_forced = []
n_fins = np.array([8,10,12,14,16])
# n_fins = np.array([])
for n in n_fins:
    T_forced.append(plot_the_temp_5('whole_system_FDO4/n_fins_test_forced', n, 10, plt))
# plt.xlim([334,365])
plt.show()

(A,B), cov = curve_fit(lin, np.append(T_b,T_a), np.append(D_b,D_a))
# xxxx = np.linspace(370,400,400)
plt.plot(xxxx, [lin(i, A, B) for i in xxxx], alpha=0.1)

mms = np.array([66*i+42 for i in n_fins])
(a1, b1), cov = curve_fit(inv_lin, mms, np.array(T_forced))
plt.plot(mms, [inv_lin(i, a1, b1) for i in mms])
plt.plot(mms, T_forced)

starting = 8
path = 'whole_system_FDO4/fin_length_test_forced'
df = pd.read_csv(f'{path}/data_s10_fin_length_30')
T_b, T_a = df['T_pr_above'][starting:]*293, df['T_pr_below'][starting:]*293
D_b, D_a = df['diff_above'][starting:], df['diff_bellow'][starting:]

plt.plot(T_b, D_b, color='r')
plt.plot(T_a, D_a, color='r')

(A,B), cov = curve_fit(lin, np.append(T_b,T_a), np.append(D_b,D_a))
xxxx = np.linspace(370,400,400)
plt.plot(xxxx, [lin(i, A, B) for i in xxxx], alpha=0.1)
print(-B/A)

starting = 6
path = 'whole_system_FDO4/fin_length_test_forced'
df = pd.read_csv(f'{path}/data_s6_fin_length_30')
T_b, T_a = df['T_pr_above'][starting:]*293, df['T_pr_below'][starting:]*293
D_b, D_a = df['diff_above'][starting:], df['diff_bellow'][starting:]

plt.plot(T_b, D_b, color='r')
plt.plot(T_a, D_a, color='r')

(A1,B1), cov = curve_fit(lin, np.append(T_b,T_a), np.append(D_b,D_a))
xxxx = np.linspace(370,400,400)
plt.plot(xxxx, [lin(i, A1, B1) for i in xxxx], alpha=0.1)
print(-B1/A1)

print('therefore the percentage difference in scaling is', round((B*A1/(B1*A)-1)*100, 1), '%')

scaling_factor = 1/(B*A1/(B1*A))
target_temp_C = 80
operational_temp = (target_temp_C + 273)*scaling_factor
print(round(1/((operational_temp-b1)/a1)), f' is the number of mm area for the system to be at {target_temp_C} degrees')
n = np.array([i+25 for i in range(50)])
l = np.array([(round(1/((operational_temp-b1)/a1)) - 13)/(2*i) - 1 for i in n])
area = np.array([(l[i] + 7)*(2*n[i]-1) for i in range(len(n))])

target_temp_C = 85
operational_temp = (target_temp_C + 273)*scaling_factor
print(round(1/((operational_temp-b1)/a1)), f' is the number of mm area for the system to be at {target_temp_C} degrees')
n = np.array([i+25 for i in range(50)])
l = np.array([(round(1/((operational_temp-b1)/a1)) - 13)/(2*i) - 1 for i in n])
area = np.array([(l[i] + 7)*(2*n[i]-1) for i in range(len(n))])