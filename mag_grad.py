import numpy as np
import scipy.constants as cts
import matplotlib.pyplot as plt
import json

mu0 = cts.mu_0

# import and choose wire gauge to analyze
with open('AWG_wire_gauges.json') as json_file:
    AWG_DIAMETERS = json.load(json_file)
wire_gauge = AWG_DIAMETERS["15"]

I0 = 1 # current of a single wire
R = 2.150*25.4*1e-3 # minmum radius of the mot coils
L = (4.800/2)*2.54*1e-2 # distance between bottoms of mot coils
AREA_CROSS_SECTION_COILS_POCKET= 2.3689 # inches squared

z = np.linspace(-4.800,4.800,200)*2.54*1e-2/2

def calc_num_wires(diameter,area_wires_current_cross_section = AREA_CROSS_SECTION_COILS_POCKET):
    area_wire_cross_section = np.pi * (diameter / 2) ** 2
    left_over= AREA_CROSS_SECTION_COILS_POCKET/area_wire_cross_section
    N = np.floor(left_over)
    return N
def ring_b_field_axial(z_offset,current=1,radius=R):
    mag_z_field = (mu0*current*radius**2)/(2*(z_offset**2+radius**2))
    return mag_z_field
def anti_helmholz_congfig(mag_top,mag_bot):
    # zeroing the center of the magnetic field
    mot_total_b_field = mag_top - mag_bot
    return mot_total_b_field

# TODO put this data in a numpy array
for wire in AWG_DIAMETERS:
    print(wire)


# define the current from the initial set current multiply by N
N = calc_num_wires(AWG_DIAMETERS["17"],area_wires_current_cross_section=AREA_CROSS_SECTION_COILS_POCKET) # number of loops

# Define the magnetic fields
mag_field_top_coil = ring_b_field_axial(L/2-z,current=N*I0,radius=R)
mag_field_bot_coil = ring_b_field_axial(L/2+z,current=N*I0,radius=R)
mot_total_b_field = anti_helmholz_congfig(mag_field_top_coil,mag_field_bot_coil)

# the magnetic field gradient, derivative of z function
dB_dz = np.diff(mot_total_b_field*10000)/np.diff(z)



plt.plot(z*1e2,mot_total_b_field*10000,'r')
plt.plot(z[:-1]*1e2,dB_dz,'b')
plt.grid()
plt.xlabel("cm")
plt.ylabel("Gauss")
plt.show()
