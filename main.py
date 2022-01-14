import swiss_projection
import numpy as np


# Example from Longitude, Latitude, Height in ETRS89 to Swiss Projection LV95
llh = [np.pi/180*8, np.pi/180*47, 400]

xyz = swiss_projection.llh2xyz(llh, "GRS80")
xyz = swiss_projection.etrs2ch(xyz)
llh = swiss_projection.xyz2llh(xyz, "Bessel1841")
ENU = np.array(swiss_projection.lv95_projection(llh))

llh = [8, 47, 400]
ENU2 = np.array(swiss_projection.wgs84_to_lv95(llh))
print(ENU-ENU2)

ENU = [3200000, 1500000, 400]


# Inverse example
llh = swiss_projection.inverse_lv95_projection(ENU)
xyz = swiss_projection.llh2xyz(llh, "Bessel1841")
xyz = swiss_projection.ch2etrs(xyz)
llh = swiss_projection.xyz2llh(xyz, "GRS80")

print(180/np.pi*llh[0])
print(180/np.pi*llh[1])
print(llh[2])

llh = [np.pi/180*8.599312654, np.pi/180*47.425664890, 524.5805000]
xyz = swiss_projection.llh2xyz(llh, "GRS80")
xyz = swiss_projection.etrs2ch(xyz)
llh = swiss_projection.xyz2llh(xyz, "Bessel1841")
ENU_fix = np.array(swiss_projection.lv95_projection(llh))
print(ENU_fix)

llh = [np.pi/180*8.599300807, np.pi/180*47.425664970, 524.5738507]
xyz = swiss_projection.llh2xyz(llh, "GRS80")
xyz = swiss_projection.etrs2ch(xyz)
llh = swiss_projection.xyz2llh(xyz, "Bessel1841")
ENU_mobile = np.array(swiss_projection.lv95_projection(llh))
print(ENU_mobile)

print(ENU_fix - ENU_mobile)

# meridian convergence
ENU = [2679520.05,  1212273.44, 400]
llh, mu = swiss_projection.inverse_lv95_projection(ENU, True)

mu = mu*200/np.pi
print(mu)