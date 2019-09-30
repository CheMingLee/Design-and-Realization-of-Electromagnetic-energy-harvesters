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

fc = 12
f = 0.1

b = (2/dm)*((pi*dm+f*p)/(p-pi*f*dm))*(2*pi/p)*(Jfw+Jbs+Jcp)
m_e = m_move + b

for i, amp in enumerate(amp_list):
    plt.figure(i)
    plt.title('{0}mm, Exp.'.format(amp), fontsize=14)
    plt.xlabel('Disp. (mm)', fontsize=14)
    plt.ylabel('Force (N)', fontsize=14)
    plt.grid()

    plt.figure(i+10)
    plt.title('{0}mm, Sim.'.format(amp), fontsize=14)
    plt.xlabel('Disp. (mm)', fontsize=14)
    plt.ylabel('Force (N)', fontsize=14)
    plt.grid()
    for j, feq in enumerate(feq_list):
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

        def X_sim(t, phi):
            return amp*sin(w*t+phi)

        popt, pcov = curve_fit(X_sim, t_exp, X_exp)
        phi = popt[0]

        dxdt = A*w*cos(w*t_exp+phi)
        dx2dt2 = -A*w**2*sin(w*t_exp+phi)
        fric = fc*np.sign(dxdt)
        F_sim = m_e*dx2dt2+fric
        
        plt.figure(i)
        plt.plot(X_exp[500:round(501+104/feq)], F_exp[500:round(501+104/feq)], '.-', label='{0}(Hz)'.format(feq))

        plt.figure(i+10)
        plt.plot(X_sim(t_exp, phi)[500:round(501+104/feq)], F_sim[500:round(501+104/feq)], label='{0}(Hz)'.format(feq))
        
    plt.figure(i)
    plt.legend(loc=1, fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylim([-60, 60])
    plt.savefig('./ch6-2-1/amp/{0}mm_Exp_F_X.png'.format(amp), bbox_inches='tight')

    plt.figure(i+10)
    plt.legend(loc=1, fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylim([-60, 60])
    plt.savefig('./ch6-2-1/amp/{0}mm_Sim_F_X.png'.format(amp), bbox_inches='tight')

#######################################################################

for i, feq in enumerate(feq_list):
    plt.figure(i+20)
    plt.title('{0}Hz, Exp.'.format(feq), fontsize=14)
    plt.xlabel('Disp. (mm)', fontsize=14)
    plt.ylabel('Force (N)', fontsize=14)
    plt.grid()

    plt.figure(i+30)
    plt.title('{0}Hz, Sim.'.format(feq), fontsize=14)
    plt.xlabel('Disp. (mm)', fontsize=14)
    plt.ylabel('Force (N)', fontsize=14)
    plt.grid()
    for j, amp in enumerate(amp_list):
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

        def X_sim(t, phi):
            return amp*sin(w*t+phi)

        popt, pcov = curve_fit(X_sim, t_exp, X_exp)
        phi = popt[0]

        dxdt = A*w*cos(w*t_exp+phi)
        dx2dt2 = -A*w**2*sin(w*t_exp+phi)
        fric = fc*np.sign(dxdt)
        F_sim = m_e*dx2dt2+fric
        
        plt.figure(i+20)
        plt.plot(X_exp[500:round(501+104/feq)], F_exp[500:round(501+104/feq)], '.-', label='{0}(mm)'.format(amp))

        plt.figure(i+30)
        plt.plot(X_sim(t_exp, phi)[500:round(501+104/feq)], F_sim[500:round(501+104/feq)], label='{0}(mm)'.format(amp))
        
    plt.figure(i+20)
    plt.legend(loc=1, fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylim([-60, 60])
    plt.savefig('./ch6-2-1/feq/{0}hz_Exp_F_X.png'.format(feq), bbox_inches='tight')

    plt.figure(i+30)
    plt.legend(loc=1, fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylim([-60, 60])
    plt.savefig('./ch6-2-1/feq/{0}hz_Sim_F_X.png'.format(feq), bbox_inches='tight')
