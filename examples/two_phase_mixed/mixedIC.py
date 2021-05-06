import numpy as np
import json
import sys
import ctypes

"input must be of the form: lx, ly, cell_Radius, number of cell A, total num cells,  output directory, BC"
"REMEMBER TO CHANGE THE FOLLOWING IN THE config.dat FILE: initial-cell-configuration-file = directoryname/init_config.json"

print(sys.argv[0])

# inputs
#lx=ctypes.c_uint(int(sys.argv[1]))
#ly=ctypes.c_uint(int(sys.argv[2]))
#cell_radius=float(sys.argv[3])
#num_cells=ctypes.c_double(float(sys.argv[4]))
#alpha=ctypes.c_double(float(sys.argv[5]))
#dir=sys.argv[6]
#bc=ctypes.c_uint(int(sys.argv[7]))

lx=int(sys.argv[1])
ly=int(sys.argv[2])
cell_radius=int(sys.argv[3])
num_cells_A=int(int(sys.argv[4]))
num_cells_B=int(int(sys.argv[5])-num_cells_A)
dir=sys.argv[6]
bc=int(sys.argv[7])

# filename
fname = dir + "/init_config.json"
print(fname)

d = 2*cell_radius
cell_centers = []
cell_alphas=[]

## Circle
#for ii in range(0, int(np.sqrt(num_cells))):
#    cell_alphas.append(alpha)
#for ii in range(0, int(np.sqrt(num_cells))):
#    for jj in range(0, int(np.sqrt(num_cells))):
#        cell_centers.append([int(lx/2 - num_cells/2 + 2*ii*cell_radius), int(ly/2 - num_cells/2 + 2*jj*cell_radius) ])
        #print([ii*cell_radius, jj*cell_radius])
## Mixed
for ii in range(0, num_cells_A):
    cell_centers.append([int( int(d*ii/ly)*2*d + d*(ii%2) ), int( ii*d %ly ) ])
    #print([int( int(d*ii/ly)*2*d + d*(ii%2) ), int( ii*d %ly ) ])
    #print([ii*cell_radius, jj*cell_radius])
print(num_cells_B)
for ii in range(0, num_cells_B):
    cell_centers.append([int( int(d*ii/ly)*2*d + d*((ii+1)%2)) , int( ii*d%ly ) ])
    #print([ii*cell_radius, jj*cell_radius])

#cell_centers.append([5,5])
#cell_centers.append([5,15])
#cell_centers.append([5,25])
#cell_centers.append([5,35])
#cell_centers.append([5,45])
#cell_centers.append([15,5])
#cell_centers.append([15,15])
#cell_centers.append([15,25])
#cell_centers.append([15,35])
#cell_centers.append([15,45])

#cell_centers.append([25,5])
#cell_centers.append([25,15])
#cell_centers.append([25,25])
#cell_centers.append([25,35])
#cell_centers.append([25,45])
#cell_centers.append([35,5])
#cell_centers.append([35,15])
#cell_centers.append([35,25])
#cell_centers.append([35,35])
#cell_centers.append([35,45])

jdata = {}
jdata['id']="initial_configuration"
jdata['BC']={"type":"unsigned","value":bc}
jdata['Size']={"type":"array(unsigned)","value":[lx,ly]}
jdata['nphases']={"type":"unsigned","value":(num_cells_A+num_cells_B)}
jdata['init_centers']={"type":"array(array(double))","value":cell_centers}
jdata['theta_pol']={"type":"array(double)","value":[0.0]*num_cells_A + [1.0]*num_cells_B}
jdata['theta_nem']={"type":"array(double)","value":[0.0]*num_cells_A + [0.0]*num_cells_B}
jdata['cell_ID']={"type":"array(double)","value":[1]*num_cells_A + [2]*num_cells_B}


#jdata['init_centers']={"type":"array(array(double))","value":cell_centers}

#jdata = json.loads(jdata)
#print(jdata)

with open(fname, 'w') as f:
    json.dump(jdata, f,separators=(',', ':'))

#print(jdata)


