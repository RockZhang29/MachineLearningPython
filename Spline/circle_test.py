import numpy as np
import math 
from matplotlib import pyplot as plt

NUMPOINT = 100


def CircleFun(p1: np.array,p2: np.array,p3: np.array):
    x1=p1[0]
    y1=p1[1]
    z1=p1[2]
    x2=p2[0]
    y2=p2[1]
    z2=p2[2]
    x3=p3[0]
    y3=p3[1]
    z3=p3[2]
 
    a1 = (y1*z2 - y2*z1 - y1*z3 + y3*z1 + y2*z3 - y3*z2)
    b1 = -(x1*z2 - x2*z1 - x1*z3 + x3*z1 + x2*z3 - x3*z2)
    c1 = (x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)
    d1 = -(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1)
 
    a2 = 2*(x2 - x1)
    b2 = 2*(y2 - y1)
    c2 = 2*(z2 - z1)
    d2 = x1*x1 + y1*y1 + z1*z1 - x2*x2 - y2*y2 - z2*z2
 
    a3 = 2*(x3 - x1)
    b3 = 2*(y3 - y1)
    c3 = 2*(z3 - z1)
    d3 = x1*x1 + y1*y1 + z1*z1 - x3*x3 - y3*y3 - z3*z3
    
    M=np.array([[a1,b1,c1,0],[a2,b2,c2,0],[a3,b3,c3,0],[0,0,0,1]])
    cp=np.dot(np.linalg.inv(M),np.array([-d1,-d2,-d3,1]))
     
    R=math.sqrt((x1-cp[0])**2+(y1-cp[1])**2+(z1-cp[2])**2)
     
    n=np.array([a1,b1,c1])
    n=n/np.linalg.norm(n) # z-vec
    a=np.array([x1-cp[0],y1-cp[1],z1-cp[2]])
    a=a/np.linalg.norm(a) # x-vec
    b=np.cross(n,a) # y-vec
    c=np.array([x2-cp[0],y2-cp[1],z2-cp[2]])
    c=c/np.linalg.norm(c) # x-vec
    d=np.cross(n,c) # y-vec
    pc = cp[0:3]
    return pc,R,a,b,c,d

def EllipseFun(p1, p2, p3):
    vec1 = p2-p1
    vec2 = p2-p3
    # len1 = math.sqrt(vec1[0]**2+vec1[1]**2+vec1[2]**2)
    # len2 = math.sqrt(vec2[0]**2+vec2[1]**2+vec2[2]**2)
    # cosa = (vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec2[2]*vec2[2])/(len1*len2)
    len1 = math.sqrt(vec1[0]**2+vec1[1]**2)
    len2 = math.sqrt(vec2[0]**2+vec2[1]**2)
    cosa = (vec1[0]*vec2[0]+vec1[1]*vec2[1])/(len1*len2)
    maxa = math.acos(cosa)
    ang = maxa*0.5
    inca = maxa*0.25
    if len1<len2:
        l1 = len2
        l2 = len1
    else:
        l1 = len1
        l2 = len2
        
    for i in range(16):
        theta = ang*0.5
        a = l1*math.sin(theta)
        b = l1*math.cos(theta)
        beta = maxa-theta
        c = l2*math.sin(beta)
        d = l2*math.cos(beta)
        v = (1-d/b)**2+(c**2)/(a**2)
        if v>1:
            ang = ang+inca
        else:
            ang = ang-inca
        inca = inca*0.5
        
    if len1<len2:
        vec = vec2
        len = len2
        pt2 = p3
    else:
        vec = vec1
        len = len1
        pt2 = p1
    
    dir = vec/len
    v = b*b/len
    h = b*a/len
    if (len1<len2 and np.cross(vec1,vec2)>0) or (len1>len2 and np.cross(vec1,vec2)<0):
        # prep = np.array([dir[2],dir[1],dir[0]]) 
        prep = np.array([dir[1],dir[0]])
    else:
        prep = np.array([dir[0],dir[1]])
        # prep = np.array([dir[0],dir[1],dir[2]])
    axis1 = np.dot(-v,dir)-np.dot(h,prep)
    center = p2 - axis1
    axis2 = pt2 - center

    print(axis1,axis2)
    print(center)
    return center
    

def judge_angle(point1,point2,pc):
    angle: float = 0.0
    vec1=np.array(pc)-np.array(point1)
    vec2=np.array(pc)-np.array(point2)
    cos_angle = np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2))
    dots = np.dot(vec1,vec2)
    norms = (np.linalg.norm(vec1)*np.linalg.norm(vec2))
    angle = np.arccos(cos_angle)
    # print(f"夹角:{np.arccos(cos_angle)*180/np.pi}")
    return angle

def generateCircle(pc,R,a,b,p1,p2):
    angle = judge_angle(p1,p2,pc)
    i=np.linspace(0,angle,NUMPOINT)
    # print(f"夹角：{angle*180/np.pi}")
    x1=pc[0] + R*np.cos(i)*a[0] +R*np.sin(i)*b[0]
    y1=pc[1] + R*np.cos(i)*a[1] +R*np.sin(i)*b[1]
    z1=pc[2] + R*np.cos(i)*a[2] +R*np.sin(i)*b[2]
    path1 = [x1,y1,z1]
    return path1

def path_blend(paths: list, points: list):
    theta = np.linspace(0,np.pi/2,NUMPOINT)
    
if __name__ == "__main__":
    p1 = np.array([-2,0])
    p2 = np.array([0,1])
    p3 = np.array([2,0])
    pc = EllipseFun(p1,p2,p3)
    pass