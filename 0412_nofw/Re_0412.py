import numpy as np
from numpy import pi, sin, cos
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os

Re_list = [6.8, 10.0, 15.0, 30.0]
feq_list = [1.0, 1.5, 2.0, 2.5, 3.0]
amp_list = [5, 7, 10]

Jfw = 0
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

            t_sim = np.arange(0, 5+1e-4, 1e-4)
            popt, pcov = curve_fit(P_Envelope, tGW_exp, P_exp, bounds=(0., 1.))
            phi_a = popt[0]
            dxdt = A*w*cos(w*t_sim+phi_a*2*pi)
            W_m = (2*pi/p)*dxdt*rg
            W_e = W_m*3*8
            i_phase = ke*W_m*cos(W_e*t_sim)/R
            P_sim = i_phase**2*Re

            plt.figure()
            plt.plot(t_exp, F_exp)
            plt.plot(t_exp, F_sim)
            plt.title('{0}(ohm),{1}(Hz),{2}(mm)'.format(Re, feq, amp))
            plt.xlabel('time (s)')
            plt.ylabel('Force (N)')
            plt.xlim([T, 3*T])
            plt.ylim([-420, 420])
            plt.legend(['Exp.', 'Sim.'])
            plt.grid()
            plt.savefig('./fig/{0}ohm/{1}hz{2}mm_F_t.png'.format(Re, feq, amp))

            plt.figure()
            plt.plot(tGW_exp, P_exp)
            plt.plot(t_sim, P_Envelope(t_sim, phi_a))
            plt.title('{0}(ohm),{1}(Hz),{2}(mm)'.format(Re, feq, amp))
            plt.xlabel('time (s)')
            plt.ylabel('Power (W)')
            plt.xlim([T, 3*T])
            plt.legend(['Exp.', 'Sim.'])
            plt.grid()
            plt.savefig('./fig/{0}ohm/{1}hz{2}mm_Pe_t.png'.format(Re, feq, amp))

            plt.figure()
            plt.plot(t_sim, P_sim, 'C1')
            plt.plot(tGW_exp, P_exp, 'C0')
            plt.title('{0}(ohm),{1}(Hz),{2}(mm)'.format(Re, feq, amp))
            plt.xlabel('time (s)')
            plt.ylabel('Power (W)')
            plt.xlim([T, 3*T])
            plt.legend(['Sim.', 'Exp.'])
            plt.grid()
            plt.savefig('./fig/{0}ohm/{1}hz{2}mm_P_t.png'.format(Re, feq, amp))

            plt.figure()
            plt.plot(X_exp[200:round(201+104/feq)], F_exp[200:round(201+104/feq)])
            plt.plot(X_sim(t_exp, phi)[200:round(201+104/feq)], F_sim[200:round(201+104/feq)])
            plt.title('{0}(ohm),{1}(Hz),{2}(mm)'.format(Re, feq, amp))
            plt.xlabel('Disp. (mm)')
            plt.ylabel('Force (N)')
            plt.ylim([-420, 420])
            plt.legend(['Exp.', 'Sim.'])
            plt.grid()
            plt.savefig('./fig/{0}ohm/{1}hz{2}mm_F_X.png'.format(Re, feq, amp))
