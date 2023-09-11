import numpy as np
import scipy.constants as cts
import matplotlib.pyplot as plt

mu0 = cts.mu_0

I = 1 # current
R = 2.150*25.4*1e-3 # minmum radius of the mot coils
L = (4.800/2)*25.4*1e-3 # distance between bottoms of mot coils

#def ring_b_field_axial(current=I,radius=R,z_offset=z_offset):
#    mag_z_field = (mu0*current*radius**2)/(2*(z_offset**2+radius**2))
#    return mag_z_field

z = np.linspace(-4.800,4.800,100)*25.4*1e-3/2
mag_z_field1 = (mu0*I*R**2)/(2*((z+L/2)**2+R**2))
mag_z_field2 = (mu0*I*R**2)/(2*((L/2-z)**2+R**2))
mot_total_b_field = mag_z_field1 + mag_z_field2

plt.plot(z/(25.4*1e-3/2),mot_total_b_field*10000,'r.')
plt.xlabel("inches")
plt.ylabel("Gauss")
plt.show()
