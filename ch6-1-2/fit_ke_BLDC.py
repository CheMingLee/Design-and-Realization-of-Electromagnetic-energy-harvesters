import numpy as np
from numpy import pi, sin, cos
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit

rg = 5
Ri = 2.7

Re_list = [1.5, 3.3, 6.8, 10.0, 15.0, 30.0, 51.0]
T_list = np.array([1.2359, 0.9838, 0.6998, 0.6014, 0.5442, 0.4734, 0.4412])
w_list = 2*pi/T_list

for Re, w in zip(Re_list, w_list):
    path_Re = '{0}ohm'.format(Re)
    data_V = os.path.join(path_Re, 'scope_2.csv')

    df = pd.read_csv(data_V, header=1)
    t_exp = df['second']
    V_exp = df['Volt']
    P_exp = V_exp**2/Re

    R = Re+Ri
    W_m = w*rg
    W_e = W_m*3*8

    def P_sim(t, ke):
        i_phase = (ke*W_m*cos(W_e*t))/R
        return i_phase**2*Re

    popt, pcov = curve_fit(P_sim, t_exp, P_exp)
    ke = popt[0]

    time = np.linspace(0, 5, 1e4)

    plt.figure()
    plt.title('Re={0} (ohm)'.format(Re))
    plt.xlabel('Time (s)')
    plt.ylabel('Power (W)')
    plt.grid()
    plt.plot(t_exp, P_exp)
    plt.plot(time, P_sim(time, ke))
    plt.xlim([0, 2])
    plt.legend(['Exp.', 'fit: ke={:.4f}'.format(ke)])
    plt.savefig('./fit/{0}ohm.png'.format(Re))
