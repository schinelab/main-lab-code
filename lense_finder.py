import numpy as np

INPUT_WAIST = 3 #mm
WAVE_LEN_NM = 399 #nm
WAVE_LEN_MM = WAVE_LEN_NM*(1e3*1e-9) # wavelength in mm

focal_length = 100 # mm

INPUT_PROBE_POWER = 2 # mW
spot_size_mm = 3 # mm WAVE_LEN*focal_length/(np.pi*INPUT_WAIST)
spot_size_cm =spot_size_mm*1e-1
rayleigh_distance = np.pi*(spot_size_cm**2)/WAVE_LEN_MM
I = INPUT_PROBE_POWER/(np.pi*(spot_size_cm)**2) # mW/cm^2



print(f"The spot size is: {spot_size_mm} in mm" )
print(f"The corresponding rayleigh length of this spot size is: {rayleigh_distance} in mm")
print(f"The corresponding intensity of this spot size: {I} in mW/cm^2")

