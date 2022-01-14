# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 20:55:06 2022

@author: Daniel
"""

import numpy as np

def angle_between_vectors(a, b):
    
    if not isinstance(a, (np.ndarray, np.generic)):
        a = np.array(a)

    if not isinstance(b, (np.ndarray, np.generic)):
        b = np.array(b)

    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    
    my_cross = np.linalg.norm(np.cross(a, b))
    my_dot = np.dot(a, b)
    
    theta = np.arctan2(my_cross, my_dot)
    return theta
