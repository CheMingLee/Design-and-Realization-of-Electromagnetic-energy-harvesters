import numpy as np
from numpy import pi, sin, cos
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os

Re_list = [6.8, 10.0, 15.0, 30.0]
feq_list = [1.0, 2.0, 3.0]
amp_list = [5, 7, 10]

Jfw = 0
Jbs = 60e-7
Jcp = 50e-7
m_move = 1051e-3
p = 10e-3
dm = 12.3e-3

Fs = 18
alpha = 0.8
vs = 2.3e-3
Fv = 2.5e-3
f = 0.1

rg = 5
Jm = 80e-7
Bm = 5e-4

b = (2/dm)*((pi*dm+f*p)/(p-pi*f*dm))*(2*pi/p)*(Jfw+Jbs+Jcp+Jm*rg)
m_e = m_move + b

Ri = 2.7
ke = 0.028
kt = ke

# y_lim = 800
rate = 1024
num = 2*rate

for Re in Re_list:
    path_Re = '{0}ohm'.format(Re)
    R = Re+Ri
    c_e = (2/dm)*((pi*dm+f*p)/(p-pi*f*dm))*(2*pi/p)*(3*ke*kt/(2*R)+Bm)*rg

    for amp in amp_list:
        plt.figure('{0}_{1}mm_Exp'.format(Re, amp))
        plt.title(r'{0}$\Omega$, {1}mm, Exp.'.format(Re, amp), fontsize=14)
        plt.xlabel('Disp. (mm)', fontsize=14)
        plt.ylabel('Force (N)', fontsize=14)
        plt.grid()

        plt.figure('{0}_{1}mm_Sim'.format(Re, amp))
        plt.title(r'{0}$\Omega$, {1}mm, Sim.'.format(Re, amp), fontsize=14)
        plt.xlabel('Disp. (mm)', fontsize=14)
        plt.ylabel('Force (N)', fontsize=-14)
        plt.grid()
        for feq in feq_list:
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
            fric = Fs*(alpha+(1-alpha)*np.exp(-np.abs(dxdt)/vs))*np.sign(dxdt)+Fv*dxdt
            F_sim = m_e*dx2dt2+c_e*dxdt+fric

            plt.figure('{0}_{1}mm_Exp'.format(Re, amp))
            plt.plot(X_exp[num:round(num+1+rate/feq)], F_exp[num:round(num+1+rate/feq)], '.-',
                     label='{0}(Hz)'.format(feq))

            plt.figure('{0}_{1}mm_Sim'.format(Re, amp))
            plt.plot(X_sim(t_exp, phi)[num:round(num+1+rate/feq)], F_sim[num:round(num+1+rate/feq)],
                     label='{0}(Hz)'.format(feq))

        plt.figure('{0}_{1}mm_Exp'.format(Re, amp))
        plt.legend(loc=1, fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        # plt.ylim([-y_lim, y_lim])
        plt.savefig('./ch6-2-3/amp/{0}ohm/{1}mm_Exp_F_X.png'.format(Re, amp), bbox_inches='tight')

        plt.figure('{0}_{1}mm_Sim'.format(Re, amp))
        plt.legend(loc=1, fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        # plt.ylim([-y_lim, y_lim])
        plt.savefig('./ch6-2-3/amp/{0}ohm/{1}mm_Sim_F_X.png'.format(Re, amp), bbox_inches='tight')

    for feq in feq_list:
        plt.figure('{0}_{1}Hz_Exp'.format(Re, feq))
        plt.title(r'{0}$\Omega$, {1}Hz, Exp.'.format(Re, feq), fontsize=14)
        plt.xlabel('Disp. (mm)', fontsize=14)
        plt.ylabel('Force (N)', fontsize=14)
        plt.grid()

        plt.figure('{0}_{1}Hz_Sim'.format(Re, feq))
        plt.title(r'{0}$\Omega$, {1}Hz, Sim.'.format(Re, feq), fontsize=14)
        plt.xlabel('Disp. (mm)', fontsize=14)
        plt.ylabel('Force (N)', fontsize=14)
        plt.grid()
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

            A = amp * 1e-3
            w = feq * 2 * pi
            T = 1 / feq

            def X_sim(t, phi):
                return amp * sin(w * t + phi)

            popt, pcov = curve_fit(X_sim, t_exp, X_exp)
            phi = popt[0]

            dxdt = A * w * cos(w * t_exp + phi)
            dx2dt2 = -A * w ** 2 * sin(w * t_exp + phi)
            fric = Fs*(alpha+(1-alpha)*np.exp(-np.abs(dxdt)/vs))*np.sign(dxdt)+Fv*dxdt
            F_sim = m_e * dx2dt2 + c_e * dxdt + fric

            plt.figure('{0}_{1}Hz_Exp'.format(Re, feq))
            plt.plot(X_exp[num:round(num+1 + rate / feq)], F_exp[num:round(num+1 + rate / feq)], '.-',
                     label='{0}(mm)'.format(amp))

            plt.figure('{0}_{1}Hz_Sim'.format(Re, feq))
            plt.plot(X_sim(t_exp, phi)[num:round(num+1 + rate / feq)], F_sim[num:round(num+1 + rate / feq)],
                     label='{0}(mm)'.format(amp))

        plt.figure('{0}_{1}Hz_Exp'.format(Re, feq))
        plt.legend(loc=1, fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        # plt.ylim([-y_lim, y_lim])
        plt.savefig('./ch6-2-3/feq/{0}ohm/{1}hz_Exp_F_X.png'.format(Re, feq), bbox_inches='tight')

        plt.figure('{0}_{1}Hz_Sim'.format(Re, feq))
        plt.legend(loc=1, fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        # plt.ylim([-y_lim, y_lim])
        plt.savefig('./ch6-2-3/feq/{0}ohm/{1}hz_Sim_F_X.png'.format(Re, feq), bbox_inches='tight')

for feq in feq_list:
    for amp in amp_list:
        plt.figure('{0}Hz_{1}mm_Exp'.format(feq, amp))
        plt.title('{0}Hz, {1}mm, Exp.'.format(feq, amp), fontsize=14)
        plt.xlabel('Disp. (mm)', fontsize=14)
        plt.ylabel('Force (N)', fontsize=14)
        plt.grid()

        plt.figure('{0}Hz_{1}mm_Sim'.format(feq, amp))
        plt.title('{0}Hz, {1}mm, Sim.'.format(feq, amp), fontsize=14)
        plt.xlabel('Disp. (mm)', fontsize=14)
        plt.ylabel('Force (N)', fontsize=14)
        plt.grid()
        for Re in Re_list:
            path_Re = '{0}ohm'.format(Re)
            R = Re+Ri
            c_e = (2/dm)*((pi*dm+f*p)/(p-pi*f*dm))*(2*pi/p)*(3*ke*kt/(2*R)+Bm)*rg

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
            fric = Fs*(alpha+(1-alpha)*np.exp(-np.abs(dxdt)/vs))*np.sign(dxdt)+Fv*dxdt
            F_sim = m_e*dx2dt2+c_e*dxdt+fric

            plt.figure('{0}Hz_{1}mm_Exp'.format(feq, amp))
            plt.plot(X_exp[num:round(num+1+rate/feq)], F_exp[num:round(num+1+rate/feq)], '.-',
                     label=r'{0}($\Omega$)'.format(Re))

            plt.figure('{0}Hz_{1}mm_Sim'.format(feq, amp))
            plt.plot(X_sim(t_exp, phi)[num:round(num+1+rate/feq)], F_sim[num:round(num+1+rate/feq)],
                     label=r'{0}($\Omega$)'.format(Re))

        plt.figure('{0}Hz_{1}mm_Exp'.format(feq, amp))
        plt.legend(loc=1, fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        # plt.ylim([-y_lim, y_lim])
        plt.savefig('./ch6-2-3/Re/{0}hz{1}mm_Exp_F_X.png'.format(feq, amp), bbox_inches='tight')

        plt.figure('{0}Hz_{1}mm_Sim'.format(feq, amp))
        plt.legend(loc=1, fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        # plt.ylim([-y_lim, y_lim])
        plt.savefig('./ch6-2-3/Re/{0}hz{1}mm_Sim_F_X.png'.format(feq, amp), bbox_inches='tight')
