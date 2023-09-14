"""
Rough magnetic field calculation
"""


import numpy as np
import scipy.constants as cts
import matplotlib.pyplot as plt
import json

mu0 = cts.mu_0
INCHES_TO_METERS = 0.0254

WIRE_GAUGE = "10"
DISTANCE_BETWEEN_COILS = (4.800)*INCHES_TO_METERS
Z_OFFSET_TOP_COIL = (DISTANCE_BETWEEN_COILS/2)*INCHES_TO_METERS
Z_OFFSET_BOTTOM_COIL = -(DISTANCE_BETWEEN_COILS/2)*INCHES_TO_METERS
CURRENT_SINGLE_LOOP = 1

# import and choose wire gauge to analyze
with open('AWG_wire_gauges.json') as json_file:
    AWG_DIAMETERS = json.load(json_file)

R = 2.150*INCHES_TO_METERS # minmum radius of the mot coils
L = (4.800)*INCHES_TO_METERS # distance between bottoms of mot coils

print(f"R: {R} ")

AREA_CROSS_SECTION_COILS_POCKET= 2.3689 # inches squared
z = np.linspace(-DISTANCE_BETWEEN_COILS/2,DISTANCE_BETWEEN_COILS/2,200)
def calc_num_wires(diameter,area_wires_current_cross_section = AREA_CROSS_SECTION_COILS_POCKET):
    area_wire_cross_section = np.pi * (diameter / 2) ** 2
    left_over= area_wires_current_cross_section/area_wire_cross_section
    N = np.floor(left_over)
    return N
def ring_b_field_axial(z_offset,current=1,radius=R):
    mag_z_field = (mu0*current*radius**2)/(2*(z_offset**2+radius**2))
    return mag_z_field
def anti_helmholz_congfig(mag_top,mag_bot):
    # zeroing the center of the magnetic field
    mot_total_b_field = mag_top - mag_bot
    return mot_total_b_field

# define the current from the initial set current multiply by N
N = calc_num_wires(AWG_DIAMETERS[WIRE_GAUGE],area_wires_current_cross_section=AREA_CROSS_SECTION_COILS_POCKET) # number of loops

# Define the magnetic fields
mag_field_top_coil = ring_b_field_axial(L/2-z,current=N*CURRENT_SINGLE_LOOP,radius=R)
mag_field_bot_coil = ring_b_field_axial(L/2+z,current=N*CURRENT_SINGLE_LOOP,radius=R)
mot_total_b_field = anti_helmholz_congfig(mag_field_top_coil,mag_field_bot_coil)

# the magnetic field gradient, derivative of z function
dB_dz = np.diff(mot_total_b_field)/(np.diff(z))
max_gradient = int(np.max(dB_dz))
y_pos_text = 0.80*max_gradient

plt.title(f"{max_gradient} G/cm using wire gauge {WIRE_GAUGE} and {N} calculated loops",fontsize=10)
#plt.title(f"Using wire gauge {WIRE_GAUGE} corresponding to {N} calculated loops",fontsize=10)
plt.plot(z*1e2,mot_total_b_field*10000,'r',label="Magnetic field of the combined coils")
plt.plot(z[:-1]*1e2,dB_dz,'b',label="Magnetic field gradient of the combined coils")
plt.text(2.5,y_pos_text,f"Current per wire: {CURRENT_SINGLE_LOOP} A")
plt.grid()
plt.xlabel("cm")
plt.ylabel("Gauss")
plt.legend()
plt.show()


print(f"field gradient {max_gradient} in G/cm")
