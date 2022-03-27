import numpy as np
from msu_esd import cross_flow_unmixed
from scipy.interpolate import interp1d

# Hot side (air)
Tin_h, Tout_h = 468, 220  # F
mdot_h = 4.7  # lbm/sec
cp_h = 0.244  # btu/(lbm F)
q = mdot_h*cp_h*(Tin_h - Tout_h)
rho_h = 0.05  # lbm/ft^3
mu_h = 1.672e-5  # lbm/(ft sec)
Pr_h = 0.701

# Cold side (water)
Tin_c = 70
cp_c = 0.997
rho_c = 62.15
mu_c = 5.46e-4
Pr_c = 5.51
mdot_c = 9.499  # lbm/sec

# Surface
Re = np.array([600, 800, 1000, 1500, 2000, 3000, 4000, 6000, 8000, 10_000])
hts = np.array([0.014, 0.012, 0.01, 0.0088, 0.008, 0.0071, 0.0067, 0.006, 0.0054, 0.0051])
hts_lamb = interp1d(Re, hts, kind='linear')

# Heat Exchanger
Cc = mdot_c*cp_c
Ch = mdot_h*cp_h
C_min, C_max = min([Cc, Ch]), max([Cc, Ch])
C = C_min/C_max
q_max = C_min*(Tin_h - Tin_c)
eff = q/q_max
ntu = cross_flow_unmixed(eff, C)
UA = ntu*C_min

# Design iteration
L1, L2, L3 = 0.45, 1.5, 1.6
vol = L1*L2*L3
Dh = 0.01352  # ft
alpha = 228  # ft^2/ft^3
A_ratio = 0.814
sigma = 0.788
Vh = mdot_h/(L2*L3*sigma*rho_h)
G = Vh*rho_h
Re_h = Vh*Dh*rho_h/mu_h

# Find heat transfer
hts_value = hts_lamb(Re_h)
h = hts_value*G*cp_h/(Pr_h**(2/3))*3600  # Btu/(ft^2 hr F)

# Get the resistance
A = alpha*vol  # ft^2
L_fin = 0.01865  # ft
k_fin = 221  # Btu/(hr ft F)
t_fin = 0.0004  # ft
m = np.sqrt(h*2/(k_fin*t_fin))
eta_fin = np.tanh(m*L_fin)/(m*L_fin)
eta_t = 1 - A_ratio*(1 - eta_fin)
R_fin = 1/(eta_t*A*h)
