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

def generateCircle(pc,R,a,b,p1,p2,p3):
    angle = judge_angle(p1,p3,pc)
    # 还要判断顺逆时针
    i=np.linspace(0,2*np.pi-angle,NUMPOINT)
    # i=np.linspace(0,angle,NUMPOINT)
    # print(f"夹角：{angle*180/np.pi}")
    x1=pc[0] + R*np.cos(i)*a[0] +R*np.sin(i)*b[0]
    y1=pc[1] + R*np.cos(i)*a[1] +R*np.sin(i)*b[1]
    z1=pc[2] + R*np.cos(i)*a[2] +R*np.sin(i)*b[2]
    path1 = [x1,y1,z1]
    ax = plt.axes(projection='3d')
    ax.plot(x1,y1,z1,lw=1,c='b')
    plt.plot(*pc,c="black", linestyle="none", marker="x")
    plt.plot(*p1,c="red", linestyle="none", marker="o")
    plt.plot(*p2,c="red", linestyle="none", marker="x")
    plt.plot(*p3,c="red", linestyle="none", marker=".")
    ax.set_xlabel('X',labelpad=5)
    ax.set_ylabel('Y',labelpad=5)
    ax.set_zlabel('Z',labelpad=5)
    ax.set_xlim3d([-5,5])
    ax.set_ylim3d([-5,5])
    ax.set_zlim3d([-5,5])
    # plt.tight_layout()
    plt.show()
    
    # return path1

def path_blend(paths: list, points: list):
    theta = np.linspace(0,np.pi/2,NUMPOINT)
    
if __name__ == "__main__":
    p1 = np.array([-6.2,2.4,3])
    p2 = np.array([5,6.3,2])
    p3 = np.array([2,-1,2.2])
    pc,R,a,b,c,d = CircleFun(p1,p2,p3)
    generateCircle(pc,R,a,b,p1,p2,p3)
    pass