import pandas as pd
import numpy as np
import math
from math import pi
from matplotlib import pyplot as plt

def read_and_sort():
    Quat = pd.read_excel(io=r'Z:\Learning Materials\Machine Learning Python\MachineLearningPython\ABB_EGM_python\abb_quat2.xlsx',usecols=['Q1','Q2','Q3','Q4'])
    quat = Quat.values
    quat.tolist()
    Joint = pd.read_excel(io=r'Z:\Learning Materials\Machine Learning Python\MachineLearningPython\ABB_EGM_python\abb_quat2.xlsx',usecols=['J1','J2','J3','J4','J5','J6'])
    joint = Joint.values
    joint.tolist()
    return quat,joint

def show_quat(quat):
    l = len(quat)
    # print(l)
    ts = np.linspace(0,1,l)
    plt.subplot(2,1,1)
    plt.title("ABB Wrist Singular Path Quat")
    plt.legend(['x','y','z','w'],loc="upper left")
    plt.plot(ts,quat)
    
    quat_init = quat[0]
    quat_diff = []
    for i in range(l):
        x = abs(quat_init[0]-quat[i][0])
        y = abs(quat_init[1]-quat[i][1])
        z = abs(quat_init[2]-quat[i][2])
        w = abs(quat_init[3]-quat[i][3])
        quad = np.array([x,y,z,w])
        quat_diff.append(quad)
        
    plt.subplot(2,1,2)
    plt.title("ABB Wrist Singular Path Quat Difference")
    plt.legend(['x','y','z','w'],loc="upper left")
    plt.plot(ts,quat_diff)
    plt.show()
            
if __name__ == "__main__":
    quat,joint = read_and_sort()
    show_quat(quat)