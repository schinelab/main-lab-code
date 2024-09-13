"""
Rough magnetic field calculation
"""


import numpy as np
import scipy.constants as cts
import matplotlib.pyplot as plt
import json

mu0 = cts.mu_0
INCHES_TO_CM = 2.54

CURRENT_SINGLE_LOOP = 1
R = 0.1 #2.150*INCHES_TO_CM # minmum radius of the mot coils
DISTANCE_BETWEEN_COILS = (4.800)*INCHES_TO_CM
AREA_CROSS_SECTION_COILS_POCKET= 2.3689*(INCHES_TO_CM)**2 # inches squared

WIRE_GAUGE = "10"
Z_OFFSET_TOP_COIL = (DISTANCE_BETWEEN_COILS/2)
Z_OFFSET_BOTTOM_COIL = -(DISTANCE_BETWEEN_COILS/2)

N = 1#calc_num_wires(AWG_DIAMETERS[WIRE_GAUGE],area_wires_current_cross_section=AREA_CROSS_SECTION_COILS_POCKET) # number of loops
z = np.linspace(-1,1,200)#np.linspace(-DISTANCE_BETWEEN_COILS/2,DISTANCE_BETWEEN_COILS/2,200)

def ring_b_field_axial(z_offset,current=1,radius=R):
    mag_z_field = (mu0*current*radius**2)/(2*(z_offset**2+radius**2)**(3/2))
    return mag_z_field
def anti_helmholz_congfig(mag_top,mag_bot):
    # zeroing the center of the magnetic field
    mot_total_b_field = mag_top - mag_bot
    return mot_total_b_field

# Define the magnetic fields
mag_field_top_coil = ring_b_field_axial(DISTANCE_BETWEEN_COILS/2-z,current=N*CURRENT_SINGLE_LOOP,radius=R)
mag_field_bot_coil = ring_b_field_axial(DISTANCE_BETWEEN_COILS/2+z,current=N*CURRENT_SINGLE_LOOP,radius=R)
mot_total_b_field = anti_helmholz_congfig(mag_field_top_coil,mag_field_bot_coil)

# the magnetic field gradient, derivative of z function
dB_dz = np.diff(mot_total_b_field*10000)/(np.diff(z))
max_gradient = int(np.max(dB_dz))
y_pos_text = 0.80*max_gradient

fig = plt.figure()
plt.title(f"{max_gradient} G/cm using wire gauge {WIRE_GAUGE} and {N} calculated loops",fontsize=10)
#plt.title(f"Using wire gauge {WIRE_GAUGE} corresponding to {N} calculated loops",fontsize=10)
plt.plot(z,mag_field_bot_coil,'b',label="one coil")
plt.text(2.5,y_pos_text,f"Current per wire: {CURRENT_SINGLE_LOOP} A")
plt.grid()
plt.xlabel("cm")
plt.ylabel("G/cm")
plt.legend()
plt.show()
fig.show()

print(f"field gradient {max_gradient} in G/cm")
