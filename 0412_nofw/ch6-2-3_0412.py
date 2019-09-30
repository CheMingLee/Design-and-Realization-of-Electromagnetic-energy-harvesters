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

    for amp in amp_list:
        plt.figure('{0}_{1}mm_Exp'.format(Re, amp))
        plt.title('{0}ohm, {1}mm, Exp.'.format(Re, amp))
        plt.xlabel('Disp. (mm)')
        plt.ylabel('Force (N)')
        plt.grid()

        plt.figure('{0}_{1}mm_Sim'.format(Re, amp))
        plt.title('{0}ohm, {1}mm, Sim.'.format(Re, amp))
        plt.xlabel('Disp. (mm)')
        plt.ylabel('Force (N)')
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
            fric = fc*np.sign(dxdt)
            F_sim = m_e*dx2dt2+c_e*dxdt+fric

            plt.figure('{0}_{1}mm_Exp'.format(Re, amp))
            plt.plot(X_exp[500:round(501+104/feq)], F_exp[500:round(501+104/feq)], '.-',
                     label='{0}(Hz)'.format(feq))

            plt.figure('{0}_{1}mm_Sim'.format(Re, amp))
            plt.plot(X_sim(t_exp, phi)[500:round(501+104/feq)], F_sim[500:round(501+104/feq)],
                     label='{0}(Hz)'.format(feq))

        plt.figure('{0}_{1}mm_Exp'.format(Re, amp))
        plt.legend()
        plt.ylim([-420, 420])
        plt.savefig('./ch6-2-3/amp/{0}ohm/{1}mm_Exp_F_X.png'.format(Re, amp))

        plt.figure('{0}_{1}mm_Sim'.format(Re, amp))
        plt.legend()
        plt.ylim([-420, 420])
        plt.savefig('./ch6-2-3/amp/{0}ohm/{1}mm_Sim_F_X.png'.format(Re, amp))

    for feq in feq_list:
        plt.figure('{0}_{1}Hz_Exp'.format(Re, feq))
        plt.title('{0}ohm, {1}Hz, Exp.'.format(Re, feq))
        plt.xlabel('Disp. (mm)')
        plt.ylabel('Force (N)')
        plt.grid()

        plt.figure('{0}_{1}Hz_Sim'.format(Re, feq))
        plt.title('{0}ohm, {1}Hz, Sim.'.format(Re, feq))
        plt.xlabel('Disp. (mm)')
        plt.ylabel('Force (N)')
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
            fric = fc * np.sign(dxdt)
            F_sim = m_e * dx2dt2 + c_e * dxdt + fric

            plt.figure('{0}_{1}Hz_Exp'.format(Re, feq))
            plt.plot(X_exp[500:round(501 + 104 / feq)], F_exp[500:round(501 + 104 / feq)], '.-',
                     label='{0}(mm)'.format(amp))

            plt.figure('{0}_{1}Hz_Sim'.format(Re, feq))
            plt.plot(X_sim(t_exp, phi)[500:round(501 + 104 / feq)], F_sim[500:round(501 + 104 / feq)],
                     label='{0}(mm)'.format(amp))

        plt.figure('{0}_{1}Hz_Exp'.format(Re, feq))
        plt.legend()
        plt.ylim([-420, 420])
        plt.savefig('./ch6-2-3/feq/{0}ohm/{1}hz_Exp_F_X.png'.format(Re, feq))

        plt.figure('{0}_{1}Hz_Sim'.format(Re, feq))
        plt.legend()
        plt.ylim([-420, 420])
        plt.savefig('./ch6-2-3/feq/{0}ohm/{1}hz_Sim_F_X.png'.format(Re, feq))

for amp in amp_list:
    for feq in feq_list:
        plt.figure('{0}mm_{1}Hz_Exp'.format(amp, feq))
        plt.title('{0}mm, {1}Hz, Exp.'.format(amp, feq))
        plt.xlabel('Disp. (mm)')
        plt.ylabel('Force (N)')
        plt.grid()

        plt.figure('{0}mm_{1}Hz_Sim'.format(amp, feq))
        plt.title('{0}mm, {1}Hz, Sim.'.format(amp, feq))
        plt.xlabel('Disp. (mm)')
        plt.ylabel('Force (N)')
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
            fric = fc*np.sign(dxdt)
            F_sim = m_e*dx2dt2+c_e*dxdt+fric

            plt.figure('{0}mm_{1}Hz_Exp'.format(amp, feq))
            plt.plot(X_exp[500:round(501+104/feq)], F_exp[500:round(501+104/feq)], '.-',
                     label='{0}(ohm)'.format(Re))

            plt.figure('{0}mm_{1}Hz_Sim'.format(amp, feq))
            plt.plot(X_sim(t_exp, phi)[500:round(501+104/feq)], F_sim[500:round(501+104/feq)],
                     label='{0}(ohm)'.format(Re))

        plt.figure('{0}mm_{1}Hz_Exp'.format(amp, feq))
        plt.legend()
        plt.ylim([-420, 420])
        plt.savefig('./ch6-2-3/Re/amp/{0}mm{1}hz_Exp_F_X.png'.format(amp, feq))

        plt.figure('{0}mm_{1}Hz_Sim'.format(amp, feq))
        plt.legend()
        plt.ylim([-420, 420])
        plt.savefig('./ch6-2-3/Re/amp/{0}mm{1}hz_Sim_F_X.png'.format(amp, feq))

for feq in feq_list:
    for amp in amp_list:
        plt.figure('{0}Hz_{1}mm_Exp'.format(feq, amp))
        plt.title('{0}Hz, {1}mm, Exp.'.format(feq, amp))
        plt.xlabel('Disp. (mm)')
        plt.ylabel('Force (N)')
        plt.grid()

        plt.figure('{0}Hz_{1}mm_Sim'.format(feq, amp))
        plt.title('{0}Hz, {1}mm, Sim.'.format(feq, amp))
        plt.xlabel('Disp. (mm)')
        plt.ylabel('Force (N)')
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
            fric = fc*np.sign(dxdt)
            F_sim = m_e*dx2dt2+c_e*dxdt+fric

            plt.figure('{0}Hz_{1}mm_Exp'.format(feq, amp))
            plt.plot(X_exp[500:round(501+104/feq)], F_exp[500:round(501+104/feq)], '.-',
                     label='{0}(ohm)'.format(Re))

            plt.figure('{0}Hz_{1}mm_Sim'.format(feq, amp))
            plt.plot(X_sim(t_exp, phi)[500:round(501+104/feq)], F_sim[500:round(501+104/feq)],
                     label='{0}(ohm)'.format(Re))

        plt.figure('{0}Hz_{1}mm_Exp'.format(feq, amp))
        plt.legend()
        plt.ylim([-420, 420])
        plt.savefig('./ch6-2-3/Re/feq/{0}hz{1}mm_Exp_F_X.png'.format(feq, amp))

        plt.figure('{0}Hz_{1}mm_Sim'.format(feq, amp))
        plt.legend()
        plt.ylim([-420, 420])
        plt.savefig('./ch6-2-3/Re/feq/{0}hz{1}mm_Sim_F_X.png'.format(feq, amp))
