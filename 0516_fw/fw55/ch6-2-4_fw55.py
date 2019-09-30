import numpy as np
from numpy import pi, sin, cos
import matplotlib.pyplot as plt
import os
import pandas as pd

Ce_csv = 'Ce_fw55.csv'

Re_list = [10.0, 30.0, 'open']
feq_list = [1.0, 2.0, 3.0]
amp_list = [5, 7, 10]

Jfw = 1425e-7
Jbs = 60e-7
Jcp = 50e-7
m_move = 1051e-3
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

x_list = Re_list[:-1]
x_list.append(50)
x = np.array(x_list)

rate = 1024
num = 2*rate

with open('./ch6-2-4/{0}'.format(Ce_csv), 'w') as f_csv:
    f_csv.write('Ohm,Hz,mm,Ce_Sim,Ceq_Sim,Ceq_Exp\n')

    for Re in Re_list:
        if Re == 'open':
            path_Re = '{0}'.format(Re)
            c_e = (2/dm)*((pi*dm+f*p)/(p-pi*f*dm))*(2*pi/p)*Bm*rg
        else:
            path_Re = '{0}ohm'.format(Re)
            R = Re+Ri
            c_e = (2/dm)*((pi*dm+f*p)/(p-pi*f*dm))*(2*pi/p)*(3*ke*kt/(2*R)+Bm)*rg

        for feq in feq_list:
            for amp in amp_list:
                path_X = '{0}hz{1}mm'.format(feq, amp)
                data_MTS = os.path.join(path_Re, path_X, 'specimen.dat')

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

                time = np.arange(0, T+0.01, 0.01)
                X_sim = A*sin(w*time)
                dxdt = A*w*cos(w*time)
                dx2dt2 = -A*w**2*sin(w*time)
                fric = fc*np.sign(dxdt)
                F_sim = m_e*dx2dt2+c_e*dxdt+fric

                Ceq_sim = (np.trapz(F_sim, x=X_sim)-4*fc*A)/(pi*w*A**2)

                y_F_exp = F_exp[num:round(num+1+rate/feq)]
                x_X_exp = X_exp[num:round(num+1+rate/feq)]*1e-3
                Ceq_exp = (np.trapz(y_F_exp, x=x_X_exp)-4*fc*A)/(pi*w*A**2)

                f_csv.write('{0},{1},{2},{3},{4},{5}\n'.format(Re, feq, amp, c_e, Ceq_sim, Ceq_exp))

df = pd.read_csv('./ch6-2-4/{0}'.format(Ce_csv))

for feq in feq_list:
    for amp in amp_list:
        data = df[(df['Hz'] == feq) & (df['mm'] == amp)]

        plt.figure()
        plt.plot(x, data['Ceq_Sim'], 'C1')
        plt.plot(x, data['Ceq_Exp'], 'o-')
        plt.title('{0}Hz, {1}mm'.format(feq, amp))
        plt.xlabel('External electrical load (ohm)')
        plt.ylabel('Energy harvester damping (Ns/m)')
        plt.xticks(x, Re_list)
        # plt.ylim([800, 3000])
        plt.legend(['Sim.', 'Exp.'])
        plt.grid()
        plt.savefig('./ch6-2-4/Ceq_{0}hz{1}mm.png'.format(feq, amp))
