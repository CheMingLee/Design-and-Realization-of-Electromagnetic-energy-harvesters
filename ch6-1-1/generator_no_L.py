import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Ri (Ohm)
Ri = 0.108
# ke (V*s/rad)
kv = 612
ke = 1/(kv*(2*np.pi/60))
# Re (Ohm)
R = np.array([10, 50, 100, 500])
Re = R*1e6/(R+1e6)
# w (rad/s)
w = np.array([4.228253907927043, 8.496531855550488, 11.127246706929610, 8.267349088394193,])
# solve
for i in range(len(R)):
    # simulation
    t = np.arange(0, 2.1, 0.001)
    ia = ke*w[i]*26/(Ri+Re[i])
    # experiment
    fileName = '{0}ohm'.format(R[i])
    df = pd.read_csv('{0}.csv'.format(fileName), header=1)
    # plot
    plt.figure()
    # plt.title('Re={0} (ohm)'.format(R[i]))
    plt.xlabel('Time (s)', fontsize=14)
    plt.ylabel('Power (W)', fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(True)
    plt.plot(df['second'], df['Volt']**2/Re[i], label='Exp.')
    plt.plot(t, ia**2*Re[i]*np.ones(len(t)), label='Sim.')
    plt.axis([0, 2, 0, 0.35])
    plt.legend(loc=3, fontsize=14)
    plt.savefig('Re{0}.png'.format(fileName))
