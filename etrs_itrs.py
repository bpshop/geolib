import json
import requests
# http://epncb.oma.be/api/production/TransformCoordinates.php

### Parameters ################################################################
X = 4331300.160
Y = 567537.081
Z = 4633133.510

vel_X = 0
vel_Y = 0
vel_Z = 0

origin_system = 'ETRF93'
target_system = 'ITRF2014'

# Implemented ETRS
# ETRF89 ETRF90 ETRF91 ETRF92 ETRF93 ETRF94 ETRF96 ETRF97 ETRF2000 ETRF2005 
# ETRF2014 

# Implemented ITRS
# ITRF88 ITRF89 ITRF90 ITRF91 ITRF92 ITRF93 ITRF94 ITRF96 ITRF97 ITRF2000 
# ITRF2005 ITRF2008 ITRF2014

epoch = '2001.0'
to_epoch = '2009.0'
### End of Parameters #########################################################


X_str = '{:}'.format(X)
Y_str = '{:}'.format(Y)
Z_str = '{:}'.format(Z)
vel_X_str = '{:}'.format(vel_X)
vel_Y_str = '{:}'.format(vel_Y)
vel_Z_str = '{:}'.format(vel_Z)



url = "https://epncb.eu/api/production/TransformCoordinates/{:}/{:}/{:}/{:}/{:}/{:}/{:}/{:}/{:}/{:}".format(
    X_str, Y_str, Z_str, origin_system, target_system, epoch, vel_X_str, vel_Y_str, vel_Z_str, to_epoch)

# response = requests.get("https://epncb.eu/api/production/TransformCoordinates/4594489.939/-678368.073/4357065.9/ITRF2005/ETRF96/2001.0/0.01/0.2/0.03/2009.0")
response = requests.get(url)

if response:
    print("Request is successful.")
    loaded_json = json.loads(response.text)
    print( "Coordinates (X,Y,Z) : %s, %s, %s" %  (loaded_json["output"]["coordinates"]["X"], loaded_json["output"]["coordinates"]["Y"], loaded_json["output"]["coordinates"]["Z"]) )
    if loaded_json["output"]["velocities"]["X"] is not None:
        print( "Velocities (X,Y,Z) : %s, %s, %s" %  (loaded_json["output"]["velocities"]["X"], loaded_json["output"]["velocities"]["Y"], loaded_json["output"]["velocities"]["Z"]) )
else:
    print("Request returned an error.")