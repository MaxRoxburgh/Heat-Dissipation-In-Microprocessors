import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv('Optimised_T_k/Corrected_func_med_2')



#%%

y = df['diff_array_big'][5:]
x = df['av_array_big'][5:]
#x = df['av_hs_array'][5:]
print(np.array(x)[-1])
plt.plot(x,y)

c,d = np.array(df['av_pr_array'][5:]), np.array(df['av_hs_array'][5:])
print(c[-1] - d[-1])

#%%

def func(x, a, b):
    return a*x + b

from scipy.optimize import curve_fit

(a, b), cov = curve_fit(func, x, y)

fy = [func(i, a, b) for i in x]

#%%
plt.plot(x,y)
plt.plot(x, fy, '--', lw=1)

x_intercept = -b/a
print(x_intercept)
print(x_intercept - np.array(x)[-1])

#%%

df = pd.read_csv('Optimised_T_k/Temp_arrays_Big')



#%%

y = df['diff_array_big'][-10:]
x = df['av_array_big'][-10:]
#x = df['av_cs_array'][5:]
print(np.array(x)[-1])
print(len(x))
plt.plot(x,y)

#%%

def func(x, a, b):
    return a*x + b

from scipy.optimize import curve_fit

(a, b), cov = curve_fit(func, x, y)

fy = [func(i, a, b) for i in x]

#%%
plt.plot(x,y)
plt.plot(x, fy, '--', lw=1)

x_intercept = -b/a
print(x_intercept)
print(x_intercept - np.array(x)[-1])

