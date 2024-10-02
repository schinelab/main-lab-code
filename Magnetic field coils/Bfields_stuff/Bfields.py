import matplotlib.pyplot as plt  # Import library for direct plotting functions
import numpy as np, sys               # Import Numerical Python
import netCDF4 as nc

pi = np.pi
mu0 = 4*pi*1e-7						#H/m
mu0_r = mu0/4/pi					#reduced vacuum permeability
in_to_cm = 2.54						#cm/inch


#defining region over which B-field is calculated
xstep = 1e-3						#m
xmax = 2*1e-2			#m
xlen = int(2*xmax/xstep)+1
xregion = np.linspace(-xmax,xmax,xlen)

ystep = xstep						#m
ymax = xmax							#m
ylen = int(2*ymax/ystep)+1
yregion = np.linspace(-ymax,ymax,ylen)

zstep = xstep						#m
zmax = xmax						#m
zlen = int(2*zmax/zstep)+1
zregion = np.linspace(-zmax,zmax,zlen)
#################################################

#open output netCDF file, define sizes of arrays, write spatial domains
def WriteIntoNetCDF(file, var_name, var_format, var_dim, var_units, var_longname, array):
	write = file.createVariable(var_name, var_format, (var_dim))
	write.units = var_units
	write.long_name = var_longname
	write[:] = array
	
newdata = nc.Dataset('Bfields.nc4', mode='w')
x_dim = newdata.createDimension('x', xlen)     # x axis
y_dim = newdata.createDimension('y', ylen)     # y axis
z_dim = newdata.createDimension('z', zlen)     # z axis
single_value = newdata.createDimension('single', 1)     # when want to record a single value
WriteIntoNetCDF(newdata, 'x-domain', 'f8', 'x', 'm', 'spatial_domain_along_x-axis', xregion)
WriteIntoNetCDF(newdata, 'y-domain', 'f8', 'y', 'm', 'spatial_domain_along_y-axis', yregion)
WriteIntoNetCDF(newdata, 'z-domain', 'f8', 'z', 'm', 'spatial_domain_along_z-axis', zregion)
#########################

#(x0,y0,z0) - coordinates of the coil center, rcoil - radius of the coil, ic - current in the coil, nturns- #of turns, ccw = 1 (0) - counter-clockwise (clockwise) flow of current
def CircularCoil(x0,y0,z0, rcoil, ic, nturns, ccw):						
	bx = np.zeros((xlen,ylen,zlen))	#to store magnetic field x-components
	by = np.zeros((xlen,ylen,zlen))	#to store magnetic field y-components
	bz = np.zeros((xlen,ylen,zlen))	#to store magnetic field z-components

	philen = 720
	phi = np.linspace(0,2*pi,philen)
	dphi = phi[1]-phi[0]				#angular element of the loop with current ic

	for i in range(xlen):
		xt = xregion[i]
		for j in range(ylen):
			yt = yregion[j]
			for k in range(zlen):
				zt = zregion[k]
				rt = np.array([xt,yt,zt])
				integrand = 0
				for p in range(philen):
					phit = phi[p]
					xl = rcoil*np.cos(phit)-x0
					yl = rcoil*np.sin(phit)-y0
					zl = z0
					l = np.array([xl,yl,zl])
					dl = rcoil*dphi/np.sqrt(xl**2 + yl**2)*np.array([-yl,xl,0])
					if ccw:
						dl = -dl
					rprime = rt - l
					rprime_abs = np.sqrt(rprime[0]**2 + rprime[1]**2 + rprime[2]**2)
					rprime3 = rprime_abs**3
					integrand = integrand + np.cross(dl, rprime)/rprime3
				#print(integrand)
				bx[i,j,k] = nturns*mu0_r*ic*integrand[0]
				by[i,j,k] = nturns*mu0_r*ic*integrand[1]
				bz[i,j,k] = nturns*mu0_r*ic*integrand[2]
		print('Done: ', (i+1)/xlen*100)

	# bzaxis = bz[int(xlen/2),int(ylen/2),:]
	# plt.plot(zregion, bzaxis*1e4)

	# bztheory = nturns*mu0*ic*rcoil*rcoil/2/((rcoil*rcoil+(zregion-z0)**2)**1.5)
	# plt.plot(zregion, bztheory*1e4)

	# plt.show()
	return bx,by,bz
#####################################################################################################################

#(x0,y0,z0) - coordinates of the beginning, (xf,yf,zf) - of the end, ic - current, nturns - #of turns, ic_dir = 1 (0) if current flows from beginning (end) to end (beginning)
def FiniteWire(x0,y0,z0, xf, yf, zf, ic, nturns, ic_dir):						
	bx = np.zeros((xlen,ylen,zlen))	#to store magnetic field x-components
	by = np.zeros((xlen,ylen,zlen))	#to store magnetic field y-components
	bz = np.zeros((xlen,ylen,zlen))	#to store magnetic field z-components
	dl = 0.5e-3						#will be used as a length step during integration over coil perimeter, m
	wire_length = np.sqrt((xf-x0)**2 + (yf-y0)**2 + (zf-z0)**2)
	nsteps = int(wire_length/dl)
	print('nsteps = ', nsteps)
	
	if ic_dir:
		dl = dl*np.array([xf-x0,yf-y0,zf-z0])/wire_length
		l0 = np.array([x0, y0, z0]) + dl/2
	else:	
		dl = -dl*np.array([xf-x0,yf-y0,zf-z0])/wire_length
		l0 = np.array([xf, yf, zf]) + dl/2

	for i in range(xlen):
		xt = xregion[i]
		for j in range(ylen):
			yt = yregion[j]
			for k in range(zlen):
				zt = zregion[k]
				rt = np.array([xt,yt,zt])
				integrand = 0
				ni = 1
				l = l0
				while ni < nsteps:
					rprime = rt - l
					rprime_abs = np.sqrt(rprime[0]**2 + rprime[1]**2 + rprime[2]**2)
					rprime3 = rprime_abs**3
					integrand = integrand + np.cross(dl, rprime)/rprime3
					l = l + dl
					ni = ni + 1
				bx[i,j,k] = nturns*mu0_r*ic*integrand[0]
				by[i,j,k] = nturns*mu0_r*ic*integrand[1]
				bz[i,j,k] = nturns*mu0_r*ic*integrand[2]
		print('Done: ', (i+1)/xlen*100)

	return bx,by,bz
#####################################################################################################################

#LOWER CIRCULAR COIL
print('Starting lower big coil')
#location of the coil
x0 = 0
y0 = 0
z0 = -8.0645*1e-2					#m		
ic = 10								#current, A
rcoil = 8.255*1e-2					#coil radius, m
nturns = 150						#number of turns
ccw = 1								#counter-clockwise current flow if looking from top
[bx,by,bz] = CircularCoil(x0,y0,z0, rcoil, ic, nturns, ccw)
WriteIntoNetCDF(newdata, 'Bx_coil_l', 'f8', ('x','y','z'), 'T', 'x-component_of_Bfield_due_to_lower_circular_coil', bx)
WriteIntoNetCDF(newdata, 'By_coil_l', 'f8', ('x','y','z'), 'T', 'y-component_of_Bfield_due_to_lower_circular_coil', by)
WriteIntoNetCDF(newdata, 'Bz_coil_l', 'f8', ('x','y','z'), 'T', 'z-component_of_Bfield_due_to_lower_circular_coil', bz)
WriteIntoNetCDF(newdata, 'z_coil_l', 'f8', 'single', 'm', 'z-coordinate_of_lower_circular_coil', z0)
print('Finished lower big coil\n')

#UPPER CIRCULAR COIL
print('Starting upper big coil')
#location of the coil
x0 = 0
y0 = 0
z0 = -z0							#m		
ic = 10								#current, A
rcoil = 8.255*1e-2					#coil radius, m
nturns = 150						#number of turns
ccw = 0								#clockwise current flow if looking from top
[bx,by,bz] = CircularCoil(x0,y0,z0, rcoil, ic, nturns, ccw)
WriteIntoNetCDF(newdata, 'Bx_coil_u', 'f4', ('x','y','z'), 'T', 'x-component_of_Bfield_due_to_upper_circular_coil', bx)
WriteIntoNetCDF(newdata, 'By_coil_u', 'f4', ('x','y','z'), 'T', 'y-component_of_Bfield_due_to_upper_circular_coil', by)
WriteIntoNetCDF(newdata, 'Bz_coil_u', 'f4', ('x','y','z'), 'T', 'z-component_of_Bfield_due_to_upper_circular_coil', bz)
WriteIntoNetCDF(newdata, 'z_coil_u', 'f8', 'single', 'm', 'z-coordinate_of_upper_circular_coil', z0)
print('Finished upper big coil\n')

xx = 9.779*1e-2; yy = 9.779*1e-2; zz = 4.826*1e-2
#RECTANGULAR COIL || yz, at x>0
print('Starting rectangular coil at x > 0')
#4 points defining the loop
x1 = xx; y1 = -yy; z1 = zz; 
x2 = xx; y2 = yy; z2 = zz
x3 = xx; y3 = yy; z3 = -zz
x4 = xx; y4 = -yy; z4 = -zz
#current
ic = 10							#current, A
nturns = 1						#number of turns
ic_dir = 1						#current flows from beginning to the end of the wire (not reverse)

[bx1,by1,bz1] = FiniteWire(x1,y1,z1, x2, y2, z2, ic, nturns, ic_dir)
[bx2,by2,bz2] = FiniteWire(x2,y2,z2, x3, y3, z3, ic, nturns, ic_dir)
[bx3,by3,bz3] = FiniteWire(x3,y3,z3, x4, y4, z4, ic, nturns, ic_dir)
[bx4,by4,bz4] = FiniteWire(x4,y4,z4, x1, y1, z1, ic, nturns, ic_dir)

bx = bx1 + bx2 + bx3 + bx4
by = by1 + by2 + by3 + by4
bz = bz1 + bz2 + bz3 + bz4

WriteIntoNetCDF(newdata, 'Bx_rect1', 'f8', ('x','y','z'), 'T', 'x-component_of_Bfield_due_to_rectangular_loop_at_positive_x0', bx)
WriteIntoNetCDF(newdata, 'By_rect1', 'f8', ('x','y','z'), 'T', 'y-component_of_Bfield_due_to_rectangular_loop_at_positive_x0', by)
WriteIntoNetCDF(newdata, 'Bz_rect1', 'f8', ('x','y','z'), 'T', 'z-component_of_Bfield_due_to_rectangular_loop_at_positive_x0', bz)
WriteIntoNetCDF(newdata, 'x_rect1', 'f8', 'single', 'm', 'x-coordinate_of_rectangular_loop_at_positive_x0', x1)
print('Finished\n')

#RECTANGULAR COIL || yz, at x<0
print('Starting rectangular coil at x < 0')
#4 points defining the loop
x1 = -xx; y1 = yy; z1 = zz; 
x2 = -xx; y2 = -yy; z2 = zz
x3 = -xx; y3 = -yy; z3 = -zz
x4 = -xx; y4 = yy; z4 = -zz
#current
ic = 10							#current, A
nturns = 1						#number of turns
ic_dir = 1						#current flows from beginning to the end of the wire (not reverse)

[bx1,by1,bz1] = FiniteWire(x1,y1,z1, x2, y2, z2, ic, nturns, ic_dir)
[bx2,by2,bz2] = FiniteWire(x2,y2,z2, x3, y3, z3, ic, nturns, ic_dir)
[bx3,by3,bz3] = FiniteWire(x3,y3,z3, x4, y4, z4, ic, nturns, ic_dir)
[bx4,by4,bz4] = FiniteWire(x4,y4,z4, x1, y1, z1, ic, nturns, ic_dir)

bx = bx1 + bx2 + bx3 + bx4
by = by1 + by2 + by3 + by4
bz = bz1 + bz2 + bz3 + bz4

WriteIntoNetCDF(newdata, 'Bx_rect3', 'f8', ('x','y','z'), 'T', 'x-component_of_Bfield_due_to_rectangular_loop_at_negative_x0', bx)
WriteIntoNetCDF(newdata, 'By_rect3', 'f8', ('x','y','z'), 'T', 'y-component_of_Bfield_due_to_rectangular_loop_at_negative_x0', by)
WriteIntoNetCDF(newdata, 'Bz_rect3', 'f8', ('x','y','z'), 'T', 'z-component_of_Bfield_due_to_rectangular_loop_at_negative_x0', bz)
WriteIntoNetCDF(newdata, 'x_rect3', 'f8', 'single', 'm', 'x-coordinate_of_rectangular_loop_at_negative_x0', x1)
print('Finished\n')

#RECTANGULAR COIL || yz, at y>0
print('Starting rectangular coil at y > 0')
#4 points defining the loop
x1 = -xx; y1 = yy; z1 = zz; 
x2 = xx; y2 = yy; z2 = zz
x3 = xx; y3 = yy; z3 = -zz
x4 = -xx; y4 = yy; z4 = -zz
#current
ic = 10							#current, A
nturns = 1						#number of turns
ic_dir = 1						#current flows from beginning to the end of the wire (not reverse)

[bx1,by1,bz1] = FiniteWire(x1,y1,z1, x2, y2, z2, ic, nturns, ic_dir)
[bx2,by2,bz2] = FiniteWire(x2,y2,z2, x3, y3, z3, ic, nturns, ic_dir)
[bx3,by3,bz3] = FiniteWire(x3,y3,z3, x4, y4, z4, ic, nturns, ic_dir)
[bx4,by4,bz4] = FiniteWire(x4,y4,z4, x1, y1, z1, ic, nturns, ic_dir)

bx = bx1 + bx2 + bx3 + bx4
by = by1 + by2 + by3 + by4
bz = bz1 + bz2 + bz3 + bz4

WriteIntoNetCDF(newdata, 'Bx_rect2', 'f8', ('x','y','z'), 'T', 'x-component_of_Bfield_due_to_rectangular_loop_at_positive_y0', bx)
WriteIntoNetCDF(newdata, 'By_rect2', 'f8', ('x','y','z'), 'T', 'y-component_of_Bfield_due_to_rectangular_loop_at_positive_y0', by)
WriteIntoNetCDF(newdata, 'Bz_rect2', 'f8', ('x','y','z'), 'T', 'z-component_of_Bfield_due_to_rectangular_loop_at_positive_y0', bz)
WriteIntoNetCDF(newdata, 'y_rect2', 'f8', 'single', 'm', 'y-coordinate_of_rectangular_loop_at_positive_y0', y1)
print('Finished\n')

#RECTANGULAR COIL || yz, at y<0
print('Starting rectangular coil at y < 0')
#4 points defining the loop
x1 = xx; y1 = -yy; z1 = zz; 
x2 = -xx; y2 = -yy; z2 = zz
x3 = -xx; y3 = -yy; z3 = -zz
x4 = xx; y4 = -yy; z4 = -zz
#current
ic = 10							#current, A
nturns = 1						#number of turns
ic_dir = 1						#current flows from beginning to the end of the wire (not reverse)

[bx1,by1,bz1] = FiniteWire(x1,y1,z1, x2, y2, z2, ic, nturns, ic_dir)
[bx2,by2,bz2] = FiniteWire(x2,y2,z2, x3, y3, z3, ic, nturns, ic_dir)
[bx3,by3,bz3] = FiniteWire(x3,y3,z3, x4, y4, z4, ic, nturns, ic_dir)
[bx4,by4,bz4] = FiniteWire(x4,y4,z4, x1, y1, z1, ic, nturns, ic_dir)

bx = bx1 + bx2 + bx3 + bx4
by = by1 + by2 + by3 + by4
bz = bz1 + bz2 + bz3 + bz4

WriteIntoNetCDF(newdata, 'Bx_rect4', 'f8', ('x','y','z'), 'T', 'x-component_of_Bfield_due_to_rectangular_loop_at_negative_y0', bx)
WriteIntoNetCDF(newdata, 'By_rect4', 'f8', ('x','y','z'), 'T', 'y-component_of_Bfield_due_to_rectangular_loop_at_negative_y0', by)
WriteIntoNetCDF(newdata, 'Bz_rect4', 'f8', ('x','y','z'), 'T', 'z-component_of_Bfield_due_to_rectangular_loop_at_negative_y0', bz)
WriteIntoNetCDF(newdata, 'y_rect4', 'f8', 'single', 'm', 'y-coordinate_of_rectangular_loop_at_negative_y0', y1)
print('Finished\n')

##############



sys.exit()			