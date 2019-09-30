import numpy as np
from numpy import pi, sin, cos

Re_list = [10.0, 30.0]
feq_list = [1.0, 2.0, 3.0]
amp_list = [5, 7, 10]

Jfw_list = [932e-7, 1425e-7]
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

for Jfw in Jfw_list:
    b = (2/dm)*((pi*dm+f*p)/(p-pi*f*dm))*(2*pi/p)*(Jfw+Jbs+Jcp+Jm*rg)
    print(b)
