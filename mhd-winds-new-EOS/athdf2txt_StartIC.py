#import h5py
#import matplotlib.pyplot as plt
import numpy as np
import math
import sys
pythonpath=str(sys.argv[4])
sys.path.insert(0, pythonpath)
import athena_read

FILE_name_f = "StartICFiles_DIR/accretion.prim."+sys.argv[2]+".athdf"
FILE_name_g = "StartICFiles_DIR/accretion.uov."+sys.argv[2]+".athdf"
IC_name     = 'Lnu_'+sys.argv[3]+'e51_Mdot_'+sys.argv[1]+'_RPNS_30km_StartIC.txt'
f_data = FILE_name_f
g_data = FILE_name_g

data_prim = athena_read.athdf(f_data)
data_uov  = athena_read.athdf(g_data)

r1       = data_prim['x1v']
vr1      = (data_prim['vel1'])[0]
P1       = (data_prim['press'])[0]
rho1     = (data_prim['rho'])[0]
T        = (data_uov['dt1'])[0]
elecfrac = (data_prim['r0'])[0]
testrho  = rho1[0,:]
testvr1  = vr1[0,:]
testT    = T[0,:]
testX    = r1[:]
testP    = P1[0,:]
Ye       = elecfrac[0,:]

Mdot = np.zeros(np.size(testrho))
for i in range(0,np.size(testrho)):
	Mdot[i]=4.0*math.pi*testrho[i]*testX[i]**2*np.abs(testvr1[i])/(2.0e33)
#print(' ')
#print('INPUT FILE: ' + IC_name)
#print('Mdot(r_out)=' + str(Mdot[np.size(testrho)-1]) + ' Msun/s')
#
#print('rho_i = ' + '{:.12e}'.format(testrho[0])                  + ' g cm^-3')
#print('rho_f = ' + '{:.12e}'.format(testrho[np.size(testrho)-1]) + ' g cm^-3')
#print('v_f   = ' + '{:.12e}'.format(testvr1[np.size(testvr1)-1]) + ' cm s^-1')
#print('p_f   = ' + '{:.12e}'.format(testP[np.size(testP)-1])     + ' erg cm^-3')
#print('T_f   = ' + '{:.12e}'.format(testT[np.size(testT)-1])     + ' K')
#print('Ye_0  = ' + '{:.12e}'.format(Ye[0]))
#print('Ye_f  = ' + '{:.12e}'.format(Ye[np.size(Ye)-1]))

c=[testX,testrho,testvr1,testT,Ye,testP]
#c=[testX,testrho,testvr1,testT]
with open(IC_name, "w") as file:
    for x in zip(*c):
        file.write("{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\n".format(*x))
        #file.write("{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\n".format(*x))


