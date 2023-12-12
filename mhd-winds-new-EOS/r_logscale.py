import numpy as np
import math
import sys
pythonpath=str(sys.argv[4])
sys.path.insert(0, pythonpath)
import athena_read

filename = 'r_scale.txt'
R_min = float(sys.argv[1])
R_max = float(sys.argv[2])
nR    = float(sys.argv[3])
scale = pow(R_max / R_min,1.0 / nR)

print('R_max = ' + '{:.4e}'.format(R_max) + ' g/cm^3')
c=[scale]
with open(filename,"w") as file:
	for x in zip(c):
		file.write(str(*x))
