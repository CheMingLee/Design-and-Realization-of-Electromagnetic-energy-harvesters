import pandas as pd
import numpy as np
from numpy import pi, sin, cos
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os

Re_list = [6.8, 10.0, 15.0, 30.0]
feq_list = [1.0, 1.5, 2.0, 2.5, 3.0]
amp_list = [5, 7, 10]

Jfw = 932e-7
Jbs = 60e-7
Jcp = 50e-7
m_move = (360+597)*1e-3
p = 10e-3
dm = 12.3e-3

fc = 12
f = 0.1

rg = 5
Jm = 80e-7
Bm = 5e-4

b = (2/dm)*((pi*dm+f*p)/(p-pi*f*dm))*(2*pi/p)*(Jfw+Jbs+Jcp+Jm*rg)
m_e = m_move + b

Ri = 2.7
ke = 0.028
kt = ke

with open('./fig/efficiency/eff.csv', 'w') as f_csv:
    f_csv.write('Ohm,Hz,mm,Eff_Elec,Eff_Mech,Eff_EH,P_in_avg_exp,P_out_avg_exp,P_out_max_exp\n')

    for Re in Re_list:
        path_Re = '{0}ohm'.format(Re)
        R = Re+Ri
        c_e = (2/dm)*((pi*dm+f*p)/(p-pi*f*dm))*(2*pi/p)*(3*ke*kt/(2*R)+Bm)*rg

        for feq in feq_list:
            for amp in amp_list:
                path_X = '{0}hz{1}mm'.format(feq, amp)
                data_MTS = os.path.join(path_Re, path_X, 'specimen.dat')
                data_GW = os.path.join(path_Re, path_X, 'DS0001.CSV')

                MTS = []
                with open(data_MTS, 'r') as data:
                    for line in data:
                        ine = line.strip('\n')
                        line = line.split('\t')
                        MTS.append(line)

                MTS = np.array(MTS[5:], dtype=np.float)
                F_exp = MTS[:, 0]
                t_exp = MTS[:, 1]
                X_exp = MTS[:, 2]

                A = amp*1e-3
                w = feq*2*pi
                T = 1/feq

                def X_sim(t, phi):
                    return amp*sin(w*t+phi)

                popt, pcov = curve_fit(X_sim, t_exp, X_exp)
                phi = popt[0]

                dxdt = A*w*cos(w*t_exp+phi)
                dx2dt2 = -A*w**2*sin(w*t_exp+phi)
                fric = fc*np.sign(dxdt)
                F_sim = m_e*dx2dt2+c_e*dxdt+fric

                GW_exp = []
                with open(data_GW, 'r') as data:
                    for line in data:
                        line = line.strip('\n')
                        line = line.split(',')
                        line.pop()
                        GW_exp.append(line)

                GW_exp = np.array(GW_exp[28:], dtype=np.float)
                tGW_exp = GW_exp[:, 0]+2.5
                V_exp = GW_exp[:, 1]
                P_exp = V_exp**2/Re

                def P_Envelope(t, phi_a):
                    dxdt = A*w*cos(w*t+phi_a*2*pi)
                    W_m = (2*pi/p)*dxdt*rg
                    i_phase = ke*W_m/R
                    return i_phase**2*Re

                popt, pcov = curve_fit(P_Envelope, tGW_exp, P_exp, bounds=(0., 1.))
                phi_a = popt[0]
                dxdt = A*w*cos(w*tGW_exp+phi_a*2*pi)
                W_m = (2*pi/p)*dxdt*rg
                W_e = W_m*3*8
                i_phase = ke*W_m*cos(W_e*tGW_exp)/R
                P_sim = i_phase**2*Re

                y_F_exp = F_exp[200:round(201+104/feq)]
                x_X_exp = X_exp[200:round(201+104/feq)]
                P_in_exp = np.trapz(y_F_exp, x=x_X_exp)*1e-3*feq

                y_P_exp = P_exp[:len(tGW_exp[tGW_exp <= T])]
                x_tGW_exp = tGW_exp[tGW_exp <= T]
                P_out_exp = 3*np.trapz(y_P_exp, x=x_tGW_exp)/T

                eff_Elec = Re/R
                eff_EH = P_out_exp/P_in_exp
                eff_Mech = eff_EH/eff_Elec

                f_csv.write('{0},{1},{2},{3},{4},{5},{6},{7},{8}\n'.format(Re, feq, amp,
                            eff_Elec, eff_Mech, eff_EH, P_in_exp, P_out_exp, P_exp.max()))

df = pd.read_csv('./fig/efficiency/eff.csv')

for feq in feq_list:
    plt.figure()
    for amp in amp_list:
        data = df[(df['Hz'] == feq) & (df['mm'] == amp)]

        plt.plot(data['Ohm'], data['Eff_Mech'], 'o-',
                 label='{0}mm'.format(amp))

    plt.title('{0}Hz'.format(feq))
    plt.xlabel('External electrical load (ohm)')
    plt.ylabel('Mechanical Efficiency')
    plt.xlim([6.8, 30])
    plt.legend()
    plt.grid()
    plt.savefig('./fig/efficiency/eff_{0}hz.png'.format(feq))

for amp in amp_list:
    plt.figure()
    for Re in Re_list:
        data = df[(df['Ohm'] == Re) & (df['mm'] == amp)]

        plt.plot(data['Hz'], data['Eff_Mech'], 'o-',
                 label='{0}ohm'.format(Re))

    plt.title('{0}mm'.format(amp))
    plt.xlabel('Excitation frequency (Hz)')
    plt.ylabel('Mechanical Efficiency')
    plt.xlim([1, 3])
    plt.ylim([0, 1])
    plt.legend()
    plt.grid()
    plt.savefig('./fig/efficiency/eff_{0}mm.png'.format(amp))
