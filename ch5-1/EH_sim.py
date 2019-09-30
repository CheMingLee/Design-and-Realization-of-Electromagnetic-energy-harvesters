import numpy as np
from numpy import pi, sin, cos
import matplotlib.pyplot as plt

Re = 10
feq = 3
amp = 10

A = amp*1e-3
w = feq*2*pi

Jfw_list = [0, 932e-7, 1425e-7]

for Jfw in Jfw_list:
    Jbs = 60e-7
    Jcp = 50e-7
    m_move = (360+597)*1e-3
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

    R = Re+Ri
    c_e = (2/dm)*((pi*dm+f*p)/(p-pi*f*dm))*(2*pi/p)*(3*ke*kt/(2*R)+Bm)*rg

    time = np.arange(0, 1+1e-4, 1e-4)
    X_sim = amp*np.sin(w*time)
    dxdt = A*w*cos(w*time)
    dx2dt2 = -A*w**2*sin(w*time)
    fric = Fs*(alpha+(1-alpha)*np.exp(-np.abs(dxdt)/vs))*np.sign(dxdt)+Fv*dxdt
    F_sim = m_e*dx2dt2+c_e*dxdt+fric
    W_m = (2*pi/p)*dxdt*rg
    W_e = W_m*3*8
    i_phase = ke*W_m*cos(W_e*time)/R
    V_sim = i_phase*Re
    P_sim = i_phase**2*Re

    plt.figure(1)
    plt.plot(time, F_sim)

    plt.figure(2)
    plt.plot(X_sim, F_sim)

    plt.figure(3)
    plt.plot(time, P_sim)

plt.figure(1)
plt.xlabel('time (s)', fontsize=14)
plt.ylabel('Force (N)', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(['without flywheel', 'with flywheel 1', 'with flywheel 2'], loc=1, fontsize=14)
plt.grid()
plt.savefig('Sim_F_t.png', bbox_inches='tight')

plt.figure(2)
plt.xlabel('Disp. (mm)', fontsize=14)
plt.ylabel('Force (N)', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(['without flywheel', 'with flywheel 1', 'with flywheel 2'], loc=1, fontsize=14)
plt.grid()
plt.savefig('Sim_F_X.png', bbox_inches='tight')

plt.figure(3)
plt.xlabel('time (s)', fontsize=14)
plt.ylabel('Power (W)', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(['without flywheel', 'with flywheel 1', 'with flywheel 2'], loc=1, fontsize=14)
plt.grid()
plt.savefig('Sim_P_t.png', bbox_inches='tight')
