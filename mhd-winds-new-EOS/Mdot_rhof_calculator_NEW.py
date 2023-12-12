import numpy as np
import math
import sys
pythonpath=str(sys.argv[4])
sys.path.insert(0, pythonpath)
import athena_read

filename = 'rho_f_'+str(sys.argv[1])+'.txt'
Mdot     = float(sys.argv[1]) # In Msun/s units
vr       = float(sys.argv[2])

Mdot_cgs = Mdot * (2e33)     # Mass accretion rate in g/s
R_out    = float(sys.argv[3]) # 1.0e8             # Outer boundary in cm
G        = 6.6743e-8         # cm^3 g^-1 s^-2
M_NS     = 1.4               # Neutron star mass in solar masses
mu       = G * M_NS * (2e33) # cm^3 s^-2

V_out   = vr #(-1.0) * np.sqrt( mu / ( 2.0 * R_out) )               # free-fall velocity at outer boundary, in cm/s
Rho_out = Mdot_cgs / (4.0 * math.pi * R_out**2 * np.abs(V_out)) # density at outer boundary, in g/cm^3

print('Rho_out = ' + '{:.4e}'.format(Rho_out) + ' g/cm^3')
c=[Rho_out]
with open(filename,"w") as file:
	for x in zip(c):
		file.write(str(*x))
