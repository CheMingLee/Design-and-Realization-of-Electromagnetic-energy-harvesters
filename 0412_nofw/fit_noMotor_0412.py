import numpy as np
from numpy import pi, sin, cos
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

feq_list = [0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
amp_list = [2, 5, 7, 10]

Jfw = 0
Jbs = 60e-7
Jcp = 50e-7
m_move = (360+597)*1e-3
p = 10e-3
dm = 12.3e-3

for feq in feq_list:
    for amp in amp_list:
        path_X = 'nofw_noMotor{0}hz{1}mm'.format(feq, amp)
        data_MTS = os.path.join('noMotor', path_X, 'specimen.dat')

        MTS = []
        with open(data_MTS, 'r') as data:
            for line in data:
                line = line.strip('\n')
                line = line.split('\t')
                MTS.append(line)

        MTS = np.array(MTS[5:], dtype=np.float)
        F_exp = MTS[:,0]
        t_exp = MTS[:,1]
        X_exp = MTS[:,2]

        A = amp*1e-3
        w = feq*2*pi
        T = 1/feq

        def X_sim(t, phi):
            return amp*sin(w*t+phi)

        popt, pcov = curve_fit(X_sim, t_exp, X_exp)
        phi = popt[0]

        def F_sim(t, f, Fs, alpha, vs, Fv):
            b = (2/dm)*((pi*dm+f*p)/(p-pi*f*dm))*(2*pi/p)*(Jfw+Jbs+Jcp)
            m_e = m_move + b
            dxdt = A*w*cos(w*t_exp+phi)
            dx2dt2 = -A*w**2*sin(w*t_exp+phi)
            fric = Fs*(alpha+(1-alpha)*np.exp(-np.abs(dxdt)/vs))*np.sign(dxdt)+Fv*dxdt
            return m_e*dx2dt2+fric

        popt, pcov = curve_fit(F_sim, t_exp, F_exp, bounds=(0, [1., 30., 1., 1e-2, 1e-2]))

        f = popt[0]
        # fc = popt[1]
        Fs = popt[1]
        alpha = popt[2]
        vs = popt[3]
        Fv = popt[4]

        plt.figure()
        plt.plot(t_exp, F_exp)
        plt.plot(t_exp, F_sim(t_exp, f, Fs, alpha, vs, Fv))
        plt.title('{0}(hz),{1}(mm)'.format(feq, amp))
        plt.xlabel('time (s)')
        plt.ylabel('Force (N)')
        plt.xlim([T, 3*T])
        plt.ylim([-60, 60])
        # plt.legend(['Exp.', 'fit: f={:.3f}, fc={:.3f}'.format(f, fc)])
        plt.legend(['Exp.', 'fit: f={:.3f}, Fs={:.3f}, alpha={:.3f}, vs={:.5f}, Fv={:.5f}'.format(f, Fs, alpha, vs, Fv)])
        plt.grid()
        plt.savefig('./fit/noMotor/{0}hz{1}mm_F_t.png'.format(feq, amp))
