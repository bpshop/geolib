import swiss_projection
import numpy as np


# Example from Longitude, Latitude, Height in ETRS89 to Swiss Projection LV95
llh = [np.pi/180*8, np.pi/180*47, 400]

xyz = swiss_projection.llh2xyz(llh, "GRS80")
xyz = swiss_projection.etrs2ch(xyz)
llh = swiss_projection.xyz2llh(xyz, "Bessel1841")
ENU = swiss_projection.lv95_projection(llh)

# Inverse example
llh = swiss_projection.inverse_lv95_projection(ENU)
xyz = swiss_projection.llh2xyz(llh, "Bessel1841")
xyz = swiss_projection.ch2etrs(xyz)
llh = swiss_projection.xyz2llh(xyz, "GRS80")

print(180/np.pi*llh[0])
print(180/np.pi*llh[1])
print(llh[2])