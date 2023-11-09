import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt
from numpy.linalg import det

from typing import Iterable

D2R = np.pi/180
R2D = 180/np.pi

def is_real(x):
    return isinstance(x, float) or isinstance(x, int)

def test_ti():
    # a,b,c,d=1.05888,-4.36731,5.73251,-2.13837
    a,b,c,d=0.0302983,0.188219,0.197121,-0.12993
    # a,b,c,d=1.6351,-5.39797,6.07606,-2.13837
    ti = np.roots([a,b,c,d])
    print(f"ti_before={ti}")
    for i in ti:
        if isinstance(i, complex) and i.imag < 0.01:
            i = i.real
            if i > 0 and i < 1:
                ti = i
                break
    if(isinstance(ti, float) == False):
        raise Exception(f"Error compute t_i:{ti}")
    print(f"ti={ti}")
    
if __name__ == "__main__":  
    test_ti()