import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Plotting import trim, plot
from Calculators import average_temp


plt.style.use('bmh')

#%%
def pr_T_analysis(T1, T2, title):
    avg_T1 = average_temp(T1)
    avg_T2 = average_temp(T2)
    
    
    
    if avg_T2 > avg_T1:
        _T1, _T2 = T2.copy(), T1.copy()
    else:
        _T1, _T2 = T1.copy(), T2.copy()
    
    
    
    fig, (ax1,ax2) = plt.subplots(2,1)
    # fig.suptitle(f"{title}".title())
    plt.text(x=0.5, y=0.84, s=f"{title}".title(), fontsize=18, ha="center", transform=fig.transFigure)
    ax1.set_title('High starting Temperature')
    ax2.set_title('Low starting Temperature')
    
    im1 = ax1.imshow(_T1, cmap='plasma')#, vmin=vmin, vmax=vmax)
    im2 = ax2.imshow(_T2, cmap='plasma')#, vmin=vmin, vmax=vmax)
    
    cbar1 = ax1.figure.colorbar(im1, orientation='horizontal')
    cbar1.ax.set_xlabel('T emperature  Scale [K]')
    cbar2 = ax2.figure.colorbar(im2, orientation='horizontal')
    cbar2.ax.set_xlabel('T emperature  Scale [K]')
    
    ax1.axes.xaxis.set_visible(False)
    ax1.axes.yaxis.set_visible(False)
    ax2.axes.xaxis.set_visible(False)
    ax2.axes.yaxis.set_visible(False)
    
    plt.show()
    
    
    avg_T1 = average_temp(_T1)
    avg_T2 = average_temp(_T2)
    print("\n\n\n-------------------------------------------------------\n")
    print(f"{title}")
    print('Average T (K) of processor from high run: ', avg_T1)
    print('Average T (K) of processor from low run:  ', avg_T2)
    
    Temp = (avg_T1+avg_T2)/2
    Temp_error = (avg_T1-avg_T2)/2
    
    return Temp, Temp_error

def T_over_time(x1,x2,y1,y2,fist_issue=False):
    
    n1 = len(y1)
    if fist_issue:
        y1_copy = y1.copy()
        y1 = [y1[0]]
        for i in range(1,n1):
            y1.append(y1_copy[i] - y1_copy[i-1])

    plt.title('whole system'.title())
    plt.plot(x1,y1)
    plt.plot(x2,y2)
    plt.yscale('log')
    plt.xlabel('Temperature [293K]')
    plt.ylabel('change from previous 2000 updates (log scale)')
    
    plt.show()
    
    
    plt.title('processor'.title())
    plt.plot(x1,y1)
    plt.plot(x2,y2)
    plt.xlabel('Temperature [293K]')
    plt.ylabel('change from previous 2000 updates')
    # plt.yscale('log')
    plt.show()


#%%
n_fins_list=[8,10,12,14,16]
T_list = []
T_err_list = []

for n_fins in n_fins_list:
    if n_fins < 11:
        first_issue = True
    else:
        first_issue = False
        
    pr_fin_4p1 = trim(np.matrix(pd.read_csv(f'4_Fins_study_data/T_k_pr_{n_fins}_fins_4p1_s3', header=None)))*293
    pr_fin_5p4 = trim(np.matrix(pd.read_csv(f'4_Fins_study_data/T_k_pr_{n_fins}_fins_5p4_s3', header=None)))*293
    data_4p1 = pd.read_csv(f'4_Fins_study_data/data_{n_fins}_fins_4p1_s3')
    data_5p4 = pd.read_csv(f'4_Fins_study_data/data_{n_fins}_fins_5p4_s3')
    
    T, T_err = pr_T_analysis(pr_fin_4p1, pr_fin_5p4, f'{n_fins} fin run')
    y1, x1 = abs(data_4p1["diff_array_big"]), data_4p1["av_array_big"]
    y2, x2 = data_5p4["diff_array_big"], data_5p4["av_array_big"]
    T_over_time(x1, x2, y1, y2, first_issue)
    
    T_list.append(T)
    T_err_list.append(T_err)
    
#%%
# x, y = n_fins_list, np.log(np.array(T_list))
x,y = n_fins_list, np.array(T_list)
# x1 = np.array([0,2,4,6])
# y1 = np.array([8700-900, 4000, 2000, 1000])
# plt.errorbar(x,y, fmt='o-')
# plt.plot(np.append(x1,x), np.append(y1,y))
# plt.yscale('log')


# def lin(x, a, b):
#     return x*a + b

def exp(x, a, b, c):
    return a*np.exp(b*x) + c
from scipy.optimize import curve_fit

#%%
plt.plot(x, y)
(a,b,c), cov = curve_fit(exp, x, y)
plt.plot(x, [exp(i, a,b,c) for i in x])


#%%
# plt.plot(x,y)
def exp_2(x, a, b):
    return a*np.exp(b*x)
dy_dx = np.array([y[i]-y[i+1] for i in range(len(x)-1)])/2
x_diff = [(x[i+1]+x[i])/2 for i in range(len(x)-1)]
plt.plot(x_diff, dy_dx)
# plt.yscale('log')

(a,b), cov = curve_fit(exp_2, x_diff, dy_dx)

#%%

plt.plot(x_diff, [exp_2(i,a,b) for i in x_diff])
plt.plot(x_diff, dy_dx, label='differential')
plt.legend()

#%%
x,y = n_fins_list, np.array(T_list)
d = a/b
plt.plot(x,y)
plt.plot(x, [exp(i, -d, b, 838) for i in x])


#%%

def difference(value):
    x,y1 = n_fins_list, np.array(T_list)
    y2 = np.array([exp(i, -d, b, value) for i in x])
    diff_array = y1-y2
    diff = sum([i for i in diff_array])
    return diff

xx = np.linspace(820,850,1000)
plt.plot(xx, [difference(i) for i in xx])



#%%
from scipy.optimize import curve_fit
def lin(x,a,b):
    return a*x + b
def one_over_x_3_4(x, a, b):
    return (a/(x))**(3/4)/293 + b

def one_over_(x, a, b):
    return (a/(x))**(7/4)/293 + b

T_points = []

def plot_the_temp(path, n, starting,ax ,col= 'r'):
    df = pd.read_csv(f'{path}/data_s6_nfins_{n}_2')
    
    T_a, T_b = df['T_above'][starting:], df['T_bellow'][starting:]
    D_a, D_b = df['diff_above'][starting:], df['diff_bellow'][starting:]
    n=len(T_a)
    ax.plot(T_a*293, D_a, color=col, label=f'{n} fins')
    ax.plot(T_b*293, D_b, color=col)
    
    (a,b), cov = curve_fit(lin, np.append(T_a,T_b), np.append(D_a,D_b))
    ax.plot(np.append(T_a,T_b)*293, [lin(i,a,b) for i in np.append(T_a,T_b)], alpha=0.2, color=col)
    ax.plot([-b*293/a],[0], 'o', color=col)
    print(f'Converge point for {n} number of fins is {-b/a}')
    print(f'difference from below to mid = {(-b/a - T_b[n-1])*293}')
    print(f'difference from above to mid = {(-b/a - T_a[n-1])*293}\n')
    return -b/a
    
x = np.array([8,10,12,14,16]) 
col = ['r','g','b','y','m']
fig = plt.figure(2,[9,5])
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

# (a,b), cov = curve_fit(one_over_x_3_4, x, T_pr)
# (a,b), cov = curve_fit(one_over_, x, T_pr)
# ax2.plot(xx, [one_over_(i, a, b)*293 for i in xx], alpha=0.8, color='r', label='fitted f($x^{-3/4}$)')

T_pr = [(T_pr_above[i] + T_pr_below[i])/2 for i in range(len(x))]
T_pr_err = [abs(T_pr_above[i] - T_pr_below[i])*293/2 for i in range(len(x))]

plt.errorbar(x, np.array(T_pr)*293, yerr=T_pr_err, fmt='o', capsize=5, color='black', label='Processor average temp')

(a,b), cov = curve_fit(one_over_x_3_4, x, T_pr, sigma=T_pr_err)
ax2.plot(xx, [one_over_x_3_4(i, a, b)*293 for i in xx], alpha=0.8, color='r', label='fitted f($x^{-3/4}$)')
ax2.legend()

print(f'The distribution apears to plataue and extrapolating the fit, the data sugests that with an infinite array of processors, the temperature could at best get to {round(b*293-273)} degrees C')
#%%
x = np.array([8,10,12,14,16])
y = np.array([9, 7.78, 6.92,6.28,5.78])
plt.plot(x,y)


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
#     above_T = 2.7
#     bellow_T = 2.8

#%%
def one_over_x_3_4(x, a, c, d):
    return (a/(x))**(3/4)/293 + c

# def one_over_x(x, a, c):
#     return (a/x) + c
# # plt.plot(x, T_points, lw=0.8)
# fit1, cov = curve_fit(one_over_x, x, T_points)

# plt.plot(x, T_points, 'x-', ms=3, lw=1)
# plt.plot(x, [one_over_x(i, fit1[0], 1) for i in x],'--', lw=0.8)

fit2, cov = curve_fit(one_over_x_3_4, x, T_points)

plt.plot(x, T_points)
plt.plot(x, [one_over_x_3_4(i, fit2[0], 1, fit2[2]) for i in x],'--', lw=0.8)

# plt.plot(x, T_points, lw=0.8)

#expected distribution 
def expectex_dist(n):
    return (7000000/(1.31*((n-1)*60+(n-1)*6+2 + 74)))**(3/4)/293 + 0.65
    
plt.plot(x, [expectex_dist(i) for i in x])
# plt.xlim([334,355])
    
    
    
    
    
    