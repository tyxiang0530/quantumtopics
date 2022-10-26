# %%
import numpy as np
import matplotlib.pyplot as plt
import math
import cmath

# %%
r_vals_i = [0, 1, 2]
alpha_i = 0

r_vals_ii = [0, 1, 2]
alpha_ii = 5
theta_ii = [0, math.pi]

delta = 1
omega_k = 1
theta_i = 1
radius = 1
k = 1
vol = 1
t_arr = np.linspace(0, 5 * math.pi, 100)
hbar = 1
# %%
# part i
'''
alpha = 0, r = 0, 1, 2
'''
fsize = 15
lsize = 10
tdir = 'in'
major = 5.0
minor = 3.0
style = 'default'

plt.style.use(style)
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = fsize
plt.rcParams['legend.fontsize'] = lsize
plt.rcParams['xtick.direction'] = tdir
plt.rcParams['ytick.direction'] = tdir
plt.rcParams['xtick.major.size'] = major
plt.rcParams['xtick.minor.size'] = minor
plt.rcParams['ytick.major.size'] = major
plt.rcParams['ytick.minor.size'] = minor

for r in r_vals_i:
    eqn_E_arr = np.zeros(len(t_arr))
    uncertainty_arr_imag = np.zeros(len(t_arr))
    uncertainty_arr_real = np.zeros(len(t_arr))
    for i in range(len(t_arr)):
        t = t_arr[i]    
        eqn_E = 1j * math.sqrt(2*math.pi * hbar * omega_k / vol) * \
                (alpha_i * cmath.exp(1j * (k * radius - omega_k * t)) - \
                        np.conj(alpha_i) * cmath.exp(-1j * (k * radius - omega_k * t)))
        eqn_E_arr[i] = eqn_E
        
        eqn_delta_E = math.sqrt((2*math.pi * hbar * omega_k / vol) * (math.cosh(2 * r) + math.sinh(2 * r) \
                * math.cos(2 * (k * radius - omega_k * t + theta_i / 2))))
        uncertainty_arr_real[i] = (eqn_delta_E.real / 2)
        
    plt.plot(t_arr, eqn_E_arr)
    plt.fill_between(t_arr, eqn_E_arr - uncertainty_arr_real, eqn_E_arr + uncertainty_arr_real, alpha = 0.5)
    plt.ylabel("$<\hat{\overrightarrow{E}}>$")
    plt.xlabel("Time")
    plt.ylim(-30, 30)
    txt="Expectation value $<\hat{\overrightarrow{E}}>$ with its uncertainty plotted over time with r = " \
            + str(r) + " and $\\alpha$ = " + str(alpha_i)
    plt.figtext(0.5, -0.03, txt, wrap=True, horizontalalignment='center', fontsize=12)
    plt.savefig("hw_6_computational_r=" + str(r) + "_alpha=" + str(alpha_i), bbox_inches = "tight")
    plt.show()

# %%
'''
alpha = 5, r = 0, 1, 2, theta = 0, pi
'''
for r in r_vals_ii:
        for theta in theta_ii:
                eqn_E_arr = np.zeros(len(t_arr))
                uncertainty_arr_imag = np.zeros(len(t_arr))
                uncertainty_arr_real = np.zeros(len(t_arr))
                for i in range(len(t_arr)):
                        t = t_arr[i]    
                        eqn_E = -2 * math.sqrt((2*math.pi * hbar * omega_k / vol) * alpha_ii * math.sin(k * r - omega_k * t + delta))
                        eqn_E_arr[i] = eqn_E
                        
                        eqn_delta_E = math.sqrt((2*math.pi * hbar * omega_k / vol) * (math.cosh(2 * r) + math.sinh(2 * r) \
                                * math.cos(2 * (k * radius - omega_k * t + theta_i / 2))))
                        uncertainty_arr_real[i] = (eqn_delta_E.real / 2)
                        uncertainty_arr_imag[i] = (eqn_delta_E.imag / 2)
                        
                plt.plot(t_arr, eqn_E_arr)
                plt.fill_between(t_arr, eqn_E_arr - uncertainty_arr_real, eqn_E_arr + uncertainty_arr_real, alpha = 0.5)
                plt.ylabel("$<\hat{\overrightarrow{E}}>$")
                plt.xlabel("Time")
                txt="Expectation value $<\hat{\overrightarrow{E}}>$ with its uncertainty plotted over time with r = " \
                        + str(r) + ", $\\theta$ = " + str(theta) + " and $\\alpha$ = " + str(alpha_ii)
                plt.figtext(0.5, -0.03, txt, wrap=True, horizontalalignment='center', fontsize=12)
                plt.ylim(-50, 50)
                plt.savefig("hw_6_computational_r=" + str(r) + "_alpha=" + str(alpha_ii) + "_theta=" 
                            + str(theta) + ".png", bbox_inches = "tight")
                plt.show()
# %%
