import matplotlib.pyplot as plt  # Import library for direct plotting functions
import numpy as np, sys               # Import Numerical Python
import netCDF4 as nc
from netCDF4 import Dataset
from scipy import interpolate

pi = np.pi
mu0 = 4*pi*1e-7						#H/m
mu0_r = mu0/4/pi					#reduced vacuum permeability
in_to_cm = 2.54						#cm/inch

data = Dataset('Bfields_same.nc4', mode='r')

x = data['x-domain'][:]
xlen = len(x)
y = data['y-domain'][:]
ylen = len(y)
z = data['z-domain'][:]
zlen = len(z)
bxl = data['Bx_coil_l'][:]
byl = data['By_coil_l'][:]
bzl = data['Bz_coil_l'][:]
bxu = data['Bx_coil_u'][:]
byu = data['By_coil_u'][:]
bzu = data['Bz_coil_u'][:]

bx = bxl+bxu*1
by = byl+byu*1
bz = bzl+bzu*1

b = np.sqrt(bx*bx + by*by + bz*bz)
btemp = b[:,int(ylen/2),:]
btemp = np.transpose(btemp)
fig, ax = plt.subplots()
pcm=ax.contourf(x, z, btemp, vmin=btemp.min(), vmax=btemp.max())
fig.colorbar(pcm)
plt.show()


bx1 = 0.5*30*data['Bx_rect1'][:]
by1 = 0.5*30*data['By_rect1'][:]
bz1 = 0.5*30*data['Bz_rect1'][:]
bx3 = 0.5*30*data['Bx_rect3'][:]
by3 = 0.5*30*data['By_rect3'][:]
bz3 = 0.5*30*data['Bz_rect3'][:]

bx = 1*bx + bx1 + bx3
by = 1*by + by1 + by3
bz = 1*bz + bz1 + bz3

b = np.sqrt(bx*bx + by*by + bz*bz)
btemp = b[:,int(ylen/2),:]
btemp = np.transpose(btemp)
fig, ax = plt.subplots()
pcm=ax.contourf(x, z, btemp, vmin=btemp.min(), vmax=btemp.max())
fig.colorbar(pcm)
plt.show()

sys.exit()



btemp = b[int(xlen/2),int(ylen/2),:]*1e4
print(x[int(xlen/2)],y[int(ylen/2)])
z0 = 100*data['z_coil_l'][:]
ztemp = z*100 - z0
fig, ax = plt.subplots()
pcm=ax.plot(ztemp, btemp)
plt.show()

dbdz = np.gradient(btemp,ztemp)
dbdz_interp = interpolate.interp1d(ztemp, dbdz)
print(dbdz_interp(7.4))
sys.exit()