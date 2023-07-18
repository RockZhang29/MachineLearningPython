import numpy as np
import math
from matplotlib import pyplot as plt
from skspatial.objects import Plane, Point, Vector
from skspatial.plotting import plot_3d

R2G:float = 180/np.pi

'''
流程1：
1. 获取所有的点位，工具相对于工件的，然后转换为末端（可能是工具，也可能是工件）相对于基坐标系的
2. 获取所有的tcp坐标，相邻的三个分为一组（暂定），通过三点确定一个平面（法向量），也可以确定一个圆形，求圆心坐标
3. 根据前两点跟圆心的夹角，来判断使用椭圆插补还是圆形插补
4. 将当前三点投影到对应平面上，3d->2d，然后进行二维的插补
5. 3d->2d，首先确定三点坐标，然后确定三点平面的方程ax+by+cz+d=0，然后求得平面法向量与平面上的点（建议取三点中的最后一点），根据公式投影出来，得到(x1,y1,z1)
6. 最后得到在当前平面上的插补轨迹，再转换回基坐标系
7. 获得对于基坐标系的所有点位，进行逆解

流程2：
1. 获取所有的点位，工具相对于工件的，然后转换为末端（可能是工具，也可能是工件）相对于基坐标系的
2. 获取所有的tcp坐标，相邻的三个分为一组（暂定），通过三点确定一个平面（法向量），也可以确定一个圆形，求圆心坐标
3. 根据前两点跟圆心的夹角，来判断使用椭圆插补还是圆形插补
4. 通过坐标系转换，以圆心为原点，将所有点通过转移矩阵转移到圆心坐标系下，再进行圆形或者椭圆插补；椭圆插补需要用向量法计算椭圆圆心
5. 获得所有点基于圆心坐标系的位置，然后再反向转移到基坐标系下
6. 获得对于基坐标系的所有点位，进行逆解

融合算法：获得path1(p1p2p3),path2(p2p3p4),path3(p3p4p5)三条路径，其中每条路径由前后两部分构成
最终路径的第一段，直接使用p1p2，第二段路径,使用如下公式: p2p3(theta) = cos(theta)^2*path1(p2p3)+sin(theta)^2*path2(p2p3)
以此类推，最后一段路径p4p5直接使用
'''

# 通过三点确定平面的法向量
def find_plane_normal(point1,point2,point3):
    plane: Plane
    vec1 = np.array(point1)-np.array(point2)
    vec2 = np.array(point1)-np.array(point3)
    a = vec1[1]*vec2[2]-vec2[1]*vec1[2]
    b = vec1[2]*vec2[0]-vec2[2]*vec1[0]
    c = vec1[0]*vec2[1]-vec2[0]*vec1[1]
    normal = [a,b,c]
    point = point3
    plane = Plane(point,normal)
    # print(f"平面法向量:{normal}")
    return normal,plane

# 返回每相邻三个点构成的圆的圆心
def calc_circle_center_and_radius(p1, p2, p3):
    x1 = p1[0]
    y1 = p1[1]
    z1 = p1[2]

    x2 = p2[0]
    y2 = p2[1]
    z2 = p2[2]

    x3 = p3[0]
    y3 = p3[1]
    z3 = p3[2]

    a1 = (y1*z2 - y2*z1 - y1*z3 + y3*z1 + y2*z3 - y3*z2)
    b1 = -(x1*z2 - x2*z1 - x1*z3 + x3*z1 + x2*z3 - x3*z2)
    c1 = (x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)
    d1 = -(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1)

    a2 = 2 * (x2 - x1)
    b2 = 2 * (y2 - y1)
    c2 = 2 * (z2 - z1)
    d2 = x1*x1 + y1*y1 + z1*z1 - x2*x2 - y2*y2 - z2*z2

    a3 = 2 * (x3 - x1)
    b3 = 2 * (y3 - y1)
    c3 = 2 * (z3 - z1)
    d3 = x1*x1 + y1*y1 + z1*z1 - x3*x3 - y3*y3 - z3*z3

    x = -(b1*c2*d3 - b1*c3*d2 - b2*c1*d3 + b2*c3*d1 + b3*c1*d2 - b3*c2*d1) /\
        (a1*b2*c3 - a1*b3*c2 - a2*b1*c3 + a2*b3*c1 + a3*b1*c2 - a3*b2*c1)
    y = (a1*c2*d3 - a1*c3*d2 - a2*c1*d3 + a2*c3*d1 + a3*c1*d2 - a3*c2*d1) /\
        (a1*b2*c3 - a1*b3*c2 - a2*b1*c3 + a2*b3*c1 + a3*b1*c2 - a3*b2*c1)
    z = -(a1*b2*d3 - a1*b3*d2 - a2*b1*d3 + a2*b3*d1 + a3*b1*d2 - a3*b2*d1) /\
        (a1*b2*c3 - a1*b3*c2 - a2*b1*c3 + a2*b3*c1 + a3*b1*c2 - a3*b2*c1)

    pc = np.array([x,y,z])
    r = math.sqrt((x1 - x)**2+(y1 - y)**2+(z1 - z)**2)
    r1 = math.dist(p1,pc)
    r2 = math.dist(p2,pc)
    r3 = math.dist(p3,pc)
    if r1-r2>1e-08 and r1-r3>1e-08:
        print(f"中心点计算错误")
        return -1,-1
    else:
        return pc, r

# 返回每相邻三个点构成的椭圆的圆心
def points2ellipses(p1:np.array,p2:np.array,p3:np.array):
    # 通过一个major axis vertice and minor axis vertice and random point on ellipse to find the center of the ellipse
    if type(p1) is not np.ndarray or type(p2) is not np.ndarray or type(p3) is not np.ndarray:
        p1 = np.array(p1)
        p2 = np.array(p2)
        p3 = np.array(p3)
    # 获得p1p2的中点pc
    pc = (p1+p2)/2
    R = (math.dist(p1,p2))/2
    norm_vec,plane = find_plane_normal(p1,p2,p3)
    # print(f"向量p1p2的中点pc: {pc}")
    print(plane)
    
    a1 = norm_vec[0]
    b1 = norm_vec[1]
    c1 = norm_vec[2]
    m1 = a1*p3[0]+b1*p3[1]+c1*p3[2]
    if a1*p3[0]+b1*p3[1]+c1*p3[2] != a1*p1[0]+b1*p1[1]+c1*p1[2] or a1*p3[0]+b1*p3[1]+c1*p3[2] != a1*p2[0]+b1*p2[1]+c1*p2[2]:
        print(f"数据错误")
        return 
    print(f"a1,b1,c1,m1{a1,b1,c1,m1}")
    
    a2 = (2*pc[0]-p1[0]-p2[0])
    b2 = (2*pc[1]-p1[1]-p2[1])
    c2 = (2*pc[2]-p1[2]-p2[2])
    m2 = pc[0]**2+pc[1]**2+pc[2]**2-R**2-p1[0]*p2[0]-p1[1]*p2[1]-p1[2]*p2[2]
    print(f"a2,b2,c2,m2{a2,b2,c2,m2}")

    
    pe = np.array([0,0,0])
    # print(f"椭圆中心{pe}")
    return pe

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
   
def circle_interpolation(center_point:list, start_point:list, end_point:list, R: float):
    x0 = center_point[0]
    y0 = center_point[1]
    z0 = center_point[2]
    p1,p2,pc = np.array(start_point), np.array(end_point), np.array(center_point)
    angle = judge_angle(p1,p2,pc)
    angle = angle*np.pi/180
    t = np.linspace(0,angle,100)
    x1,y1,z1 = [],[],[]
    x = x0+np.cos(t)*R+np.sin(t)*0
    y = y0+np.cos(t)*0+np.sin(t)*R
    z = z0+np.cos(t)*R+np.sin(t)*1

    x1.append(x)
    y1.append(y)
    z1.append(z)
    ax = plt.axes(projection='3d')
    ax.plot(x,y,z,lw=1,c='r')
    ax.set_xlabel('X',labelpad=5)
    ax.set_ylabel('Y',labelpad=5)
    ax.set_zlabel('Z',labelpad=5)
    # plt.tight_layout()
    plt.show()
    return x1,y1,z1

def ellipsis_interpolate(center_point:list, long:list, short:list):
    x0 = center_point[0]
    y0 = center_point[1]
    z0 = center_point[2]
    
    t = np.linspace(-np.pi,np.pi,100)
    x1,y1,z1 = [],[],[]

    x = x0+long[0]*np.cos(t)+short[0]*np.sin(t)
    y = y0+long[1]*np.cos(t)+short[1]*np.sin(t)
    z = z0+long[2]*np.cos(t)+short[2]*np.sin(t)
    
    x1.append(x)
    y1.append(y)
    z1.append(z)
    ax = plt.axes(projection='3d')
    ax.plot(x,y,z,lw=1,c='b')
    plt.plot(*center_point,c="red", linestyle="none", marker="o")
    plt.plot(*long,c="red", linestyle="none", marker="o")
    plt.plot(*short,c="red", linestyle="none", marker="o")
    ax.set_xlabel('X',labelpad=5)
    ax.set_ylabel('Y',labelpad=5)
    ax.set_zlabel('Z',labelpad=5)
    plt.tight_layout()
    plt.show()

def path_blend(path1: list, path2:list):
    path = []
    if len(path1) != len(path2):
        print("Path point num error")
        return 
    NUM_POINTS = len(path1)
    
    steps = np.linspace(0,np.pi/2,NUM_POINTS)
    for i in range(NUM_POINTS):
        x = np.cos(steps[i])*np.cos(steps[i])*path1[i][0]+(np.sin(steps[i])*np.sin(steps[i]))*path2[i][0]
        y = np.cos(steps[i])*np.cos(steps[i])*path1[i][1]+(np.sin(steps[i])*np.sin(steps[i]))*path2[i][1]
        z = np.cos(steps[i])*np.cos(steps[i])*path1[i][2]+(np.sin(steps[i])*np.sin(steps[i]))*path2[i][2]
        pts = [x,y,z]
        path.append(pts)
    return path

def path_plot(paths: list, POINTS):
    x,y,z = [],[],[]
    for i in range(len(paths)):
        for j in range(len(paths[i])):
            x.append(paths[i][j][0])
            y.append(paths[i][j][1])
            z.append(paths[i][j][2])     
    ax = plt.axes(projection='3d')
    ax.plot(x,y,z,lw=1,c='b')
    ax.set_xlabel('X',labelpad=5)
    ax.set_ylabel('Y',labelpad=5)
    ax.set_zlabel('Z',labelpad=5)
    plt.plot(*zip(*POINTS), c="yellow", linestyle="none", marker="x")
    plt.tight_layout()
    
    var = compare_difference(POINTS,x,y,z)
    plt.figure()
    plt.plot(var, 'r-', linewidth=1)
    plt.grid(True) 
    plt.show()

def test(POINTS: list):
    # 第一步，点位分组
    if len(POINTS)<4:
        print("点位数量过少，用MoveC吧")
        return 
    paths = []
    for i in range(len(POINTS)-2):
        p1 = POINTS[i]
        p2 = POINTS[i+1]
        p3 = POINTS[i+2]
        angle = R2G*judge_angle(p1,p2,p3)
        #if angle == 180 or angle == 0:
        #     print(f"共线")
        #     return 
        # if angle>90 and angle<180:
        #     path = ellipse_interpolation(p1,p2,p3)
        # if angle<90 and angle>0:
        #     path = circle_interpolation(p1,p2,p3)
        pc,r = calc_circle_center_and_radius(p1,p2,p3)
        print(pc)
        path = circle_interpolation(pc,p1,p3,r)
        
        paths.append(path)
    

if __name__ == "__main__":
    p1 = np.array([3,-2,0])
    p2 = np.array([2,2,0])
    p3 = np.array([4,0,0])
    p4 = np.array([-2,1,0])
    POINTS = [p1,p2,p3,p4]
    # test(POINTS)
    
    
    

