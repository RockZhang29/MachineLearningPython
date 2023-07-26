import numpy as np
import random
import math
from matplotlib import pyplot as plt
from skspatial.objects import Plane, Point, Vector
from skspatial.plotting import plot_3d
from sympy import symbols, Eq, solve

D2R = np.pi/180
R2D = 180/np.pi
ITER = 32
NUMPOINT = 100

### 工具函数
# 合并函数
def mix_path(path1,path2):
    x = np.append(path1[0],path2[0])
    y = np.append(path1[1],path2[1])
    z = np.append(path1[2],path2[2])
    path = [x,y,z]
    return path

# 互换函数
def switch(a,b):
    x,y = b,a
    return x,y

#相等函数
def is_equal(a,b):
    if abs(a-b)<1e-07:
        return True
    else:
        # print(abs(a-b))
        return False

# 通过三点确定平面的法向量
def find_plane_normal(point1,point2,point3):
    vec1 = np.array(point1)-np.array(point2)
    vec2 = np.array(point1)-np.array(point3)
    a = vec1[1]*vec2[2]-vec2[1]*vec1[2]
    b = vec1[2]*vec2[0]-vec2[2]*vec1[0]
    c = vec1[0]*vec2[1]-vec2[0]*vec1[1]
    normal = [a,b,c]
    point = point3
    plane = Plane(point,normal)
    # print(f"平面法向量:{normal}")
    return normal

# 返回每相邻三个点构成的圆的圆心
def CircleFun(p1,p2,p3):
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

# 圆弧插补
def circle_interpolate(pc,R,a,b,c,d,p1,p2,p3):
    # 需要判断插补的点,p1,p2,p3
    a1 = find_angle_of_3_pt(p2,p3,p1)
    if abs(a1)<np.pi/2:
        a2 = find_angle_of_3_pt(p1,p2,p3)
        if a2<=np.pi/2:
            angle1 = find_angle_of_3_pt(p1,p2,pc)
            angle2 = find_angle_of_3_pt(p2,p3,pc)
            i=np.linspace(0,angle1,NUMPOINT)
            j=np.linspace(angle1,angle1+angle2,NUMPOINT)
        else:
            angle1 = find_angle_of_3_pt(p1,p2,pc)
            angle3 = find_angle_of_3_pt(p1,p3,pc)
            angle4 = find_angle_of_3_pt(p2,p3,pc)
            angle2 = 2*np.pi-angle3-angle1
            st_angle = 2*np.pi-angle1
            i=np.linspace(0,st_angle,NUMPOINT)
            j=np.linspace(st_angle,st_angle+angle4,NUMPOINT)
    else:
        angle1 = find_angle_of_3_pt(p1,p2,pc)
        angle3 = find_angle_of_3_pt(p1,p3,pc)
        angle2 = 2*np.pi-angle3-angle1
        # angle = 2*np.pi-angle3
        i=np.linspace(0,angle1,NUMPOINT)
        j=np.linspace(angle1,angle1+angle2,NUMPOINT)
        # print(angle*R2D,angle2*R2D,angle3*R2D)
    # 拆解为n和n+1段


    
    x1=pc[0] + R*np.cos(i)*a[0] +R*np.sin(i)*b[0]
    y1=pc[1] + R*np.cos(i)*a[1] +R*np.sin(i)*b[1]
    z1=pc[2] + R*np.cos(i)*a[2] +R*np.sin(i)*b[2]
    path1 = [x1,y1,z1]
    
    
    # print(f"j={j}")
    x2=pc[0] + R*np.cos(j)*a[0] +R*np.sin(j)*b[0]
    y2=pc[1] + R*np.cos(j)*a[1] +R*np.sin(j)*b[1]
    z2=pc[2] + R*np.cos(j)*a[2] +R*np.sin(j)*b[2]
    path2 = [x2,y2,z2]
    # print(path2)
    path = mix_path(path1,path2)
    
    # plot_path(path1,pc,p1,p2,p3)
    # plot_path(path2,pc,p1,p2,p3)
    return path

# 二分法取得角度
def bisect(p1,p2,p3):
    # i = 1,2,3...
    vec1 = p1-p2 # pi-pj, j=i-1+points.length
    vec2 = p3-p2 # pi-pk, k=i+1
    len1 = math.sqrt(vec1[0]**2+vec1[1]**2+vec1[2]**2)
    len2 = math.sqrt(vec2[0]**2+vec2[1]**2+vec2[2]**2)
    cos_angle = (vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2])/(len1*len2)
    theta = np.arccos(cos_angle)
    angle = []
    ang = 0.5*theta
    incA = 0.25*theta

    if len1<len2:
        l1 = len2
        l2 = len1
    else:
        l1 = len1
        l2 = len2

    for i in range(ITER):
        alpha = 0.5*ang
        a = l1*np.sin(alpha)
        b = l1*np.cos(alpha)
        beta = theta-alpha
        c = l2*np.sin(beta)
        d = l2*np.cos(beta)
        v = (1-d/b)**2+(c/a)**2
        if(v>1):
            ang = ang+incA
        else:
            ang = ang-incA
        incA = incA*0.5
        angle.append(alpha)
        
    #print(f"theta1={angle[-1]*R2D}")
    theta_approximate = angle[-1]
    return theta_approximate

# 判断夹角角度
def find_angle_of_3_pt(p1,p2,pc):
    vec1=np.array(pc)-np.array(p1)
    vec2=np.array(pc)-np.array(p2)
    # print(f"vec1={vec1},vec2={vec2}")
    cos_angle = np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2))
    if is_equal(abs(cos_angle),1.0):
        cos_angle = round(cos_angle)
    dots = np.dot(vec1,vec2)
    norms = (np.linalg.norm(vec1)*np.linalg.norm(vec2))
    # print(f"dots={dots},norms={norms}")
    angle = np.arccos(cos_angle)
    # print(f"夹角:{np.arccos(cos_angle)*180/np.pi}")
    return angle

# 返回相邻三点构成的椭圆圆心
def calc_ellipse_center(p1,p2,p3):
    # 调用三点平面生成的法向量
    n = find_plane_normal(p1,p2,p3)
    ## 方程1
    Index1 = n[0]*p3[0]+n[1]*p3[1]+n[2]*p3[2]
    a1 = n[0]
    a2 = n[1]
    a3 = n[2]
    # print(f"n={n}")
    ## 方程2
    theta = find_angle_of_3_pt(p1,p3,p2) # p1,p2,pc
    theta1 = bisect(p1,p2,p3)
    # print(f"theta={theta*R2D}")
    # print(f"theta1={theta1*R2D},theta2={theta2*R2D}")
    len1 = math.dist(p1,p2)
    len2 = len1*np.cos(theta1)
    len3 = len1*np.sin(theta1)
    mm = len1**2+len2**2-2*len1*len2*np.cos(theta1)
    Index2 = p1[0]**2+p1[1]**2+p1[2]**2-p1[0]*p2[0]-p1[1]*p2[1]-p1[2]*p2[2]-mm
    b1 = p1[0]-p2[0]
    b2 = p1[1]-p2[1]
    b3 = p1[2]-p2[2]
    # print(f"b={b1,b2,b3}")
    
    ## 方程3
    theta2 = theta-theta1
    # print(f"theta={theta*R2D}")
    # print(f"theta1={theta1*R2D}")
    # print(f"theta2={theta2*R2D}")
    len4 = math.dist(p2,p3)
    len5 = len4*np.sin(theta2)
    len6 = len4*np.cos(theta2)
    len7 = len2-len6
    # print(f"L = {len4,len5,len6,len7}")
    c1 = 2*(p1[0]-p3[0])
    c2 = 2*(p1[1]-p3[1])
    c3 = 2*(p1[2]-p3[2])
    Index3 = len7**2+len5**2-len3**2+p1[0]**2+p1[1]**2+p1[2]**2-p3[0]**2-p3[1]**2-p3[2]**2
    
    Left = np.array([[a1,a2,a3],[b1,b2,b3],[c1,c2,c3]])
    Right = np.array([[Index1],[Index2],[Index3]])
    # print(f"左:{Left}")
    # print(f"右:{Right}")
    ## 有可能会形成奇异矩阵
    ans = np.linalg.solve(Left,Right)
    
    po = np.array([ans[0][0],ans[1][0],ans[2][0]])
    # print(f"Solution = {po}")
    return po

# 椭圆插补
def ellipsis_interpolate(po, long, short, ps, pe, p2):
    # angle是pe和p2的夹角，总的角度需要再加上ps和p2的夹角(90°)
    # angle = np.pi/2 + angle
    Left = np.array([[long[0],short[0]],[long[1],short[1]]])
    Right = np.array([[pe[0]-po[0]],[pe[1]-po[1]]])
    ans = np.linalg.solve(Left,Right)
    # print(ans)
    sint = ans[0][0]
    cost = ans[1][0]
    if is_equal(sint**2+cost**2,1.0) is False:
        print(f"计算具体角度出错")
        return False
    
    # 若90°
    a1 = np.arcsin(sint)
    # a2 = np.arccos(cost)
    
    # if a2>np.pi/2:
    #     a2 = np.pi-a2
    # print(a1*R2D,a2*R2D)
    # 判断插补方向
    vec1 = ps-p2
    vec2 = p2-pe
    way = np.cross(vec1,vec2)[0]+np.cross(vec1,vec2)[1]+np.cross(vec1,vec2)[2]
    long,short = switch(long,short)
    if way==0:
        print("三点共线")
        return False
    
    ## 从ps插补到pe
    angle = np.pi-a1
    x1,y1,z1 = [],[],[]  
    x2,y2,z2 = [],[],[]  
    t = np.linspace(0,np.pi/2,NUMPOINT)
    k = np.linspace(np.pi/2,angle,NUMPOINT)

    for i in t:
        x = po[0]+long[0]*np.cos(i)+short[0]*np.sin(i)
        y = po[1]+long[1]*np.cos(i)+short[1]*np.sin(i)
        z = po[2]+long[2]*np.cos(i)+short[2]*np.sin(i)
        x1.append(x)
        y1.append(y)
        z1.append(z)
    
    for i in k:
        x = po[0]+long[0]*np.cos(i)+short[0]*np.sin(i)
        y = po[1]+long[1]*np.cos(i)+short[1]*np.sin(i)
        z = po[2]+long[2]*np.cos(i)+short[2]*np.sin(i)
        x2.append(x)
        y2.append(y)
        z2.append(z)
    
    path1 = [x1,y1,z1]
    # plot_path(path1,po,ps,p2,pe)
    path2 = [x2,y2,z2]
    # plot_path(path2,po,ps,p2,pe)
    # 合并两条path
    path = mix_path(path1,path2)
    # plot_path(path,po,ps,p2,pe)
    return path

# 直线插补
def straight_interploate(p1,p2,p3):
    d1 = math.dist(p1,p2)
    d2 = math.dist(p2,p3)
    vec1 = p2-p1
    vec2 = p3-p2
    xa,ya,za = [],[],[]
    xb,yb,zb = [],[],[]
    for i in range(NUMPOINT):
        xa.append(p1[0]+(i+1)*vec1[0]/NUMPOINT)
        ya.append(p1[1]+(i+1)*vec1[1]/NUMPOINT)
        za.append(p1[2]+(i+1)*vec1[2]/NUMPOINT)
        xb.append(p2[0]+(i+1)*vec2[0]/NUMPOINT)
        yb.append(p2[1]+(i+1)*vec2[1]/NUMPOINT)
        zb.append(p2[2]+(i+1)*vec2[2]/NUMPOINT)

    path1 = [xa,ya,za]
    path2 = [xb,yb,zb]
    path = mix_path(path1,path2)
    return path 
    

# 单一路径展示
def plot_path(path,po,ps,p2,pe):
    if len(path) > 3:
        x = path[0]+path[3]
        y = path[1]+path[4]
        z = path[2]+path[5]
    else:
        x = path[0]
        y = path[1]
        z = path[2]
    
    ax = plt.axes(projection='3d')
    ax.plot(x,y,z,lw=1,c='b')
    plt.plot(*po,c="black", linestyle="none", marker=".",label="ellipse center")
    plt.plot(*ps,c="red", linestyle="none", marker=".",label="start point")
    plt.plot(*p2,c="red", linestyle="none", marker="x",label="middle point")
    plt.plot(*pe,c="red", linestyle="none", marker="+",label="end point")
    ax.set_xlabel('X',labelpad=5)
    ax.set_ylabel('Y',labelpad=5)
    ax.set_zlabel('Z',labelpad=5)
    ax.set_xlim3d([-5,5])
    ax.set_ylim3d([-5,5])
    ax.set_zlim3d([-5,5])
    plt.legend(loc="best",fontsize=6)
    plt.show()

# 融合
def path_blend(path1,path2):
    ts = np.linspace(0,np.pi/2,NUMPOINT)
    path_blend = []
    px,py,pz = [],[],[]
    for i in range(NUMPOINT):
        x  = path1[0][NUMPOINT+i]*np.cos(ts[i])**2+path2[0][i]*np.sin(ts[i])**2
        y  = path1[1][NUMPOINT+i]*np.cos(ts[i])**2+path2[1][i]*np.sin(ts[i])**2
        z  = path1[2][NUMPOINT+i]*np.cos(ts[i])**2+path2[2][i]*np.sin(ts[i])**2
        # point = path1[NUMPOINT+i]*np.sin(ts[i])**2+path2[i]*np.cos(ts[i])**2
        px.append(x)
        py.append(y)
        pz.append(z)
        # point = [x,y,z]
        # path_blend.append(point)
    path_blend = [px,py,pz]
    if len(path_blend[0]) != NUMPOINT:
        return False
    
    return path_blend

# 路径重整
def path_reform(paths,points):
    # 如果4个点，总共三段线段，i=0时得到p0p1,i=1时融合p1p2,i=2时得到p2p3
    traj = []
    pn = len(points)
    for i in range(pn-1):
        if i==0:
            # p1->p2直接使用，paths[0]的前半段
            new_path = paths[0] # [x,y,z]
            x = new_path[0][:NUMPOINT]
            y = new_path[1][:NUMPOINT]
            z = new_path[2][:NUMPOINT]
            ini_path = [x,y,z]
            traj.append(ini_path)
            # print(f"开始路径{new_path}")
        elif i==pn-2:
            # pn-1->pn直接使用，paths[-1]的后半段
            new_path = paths[-1]
            x = new_path[0][NUMPOINT:]
            y = new_path[1][NUMPOINT:]
            z = new_path[2][NUMPOINT:]
            end_path = [x,y,z]
            traj.append(end_path)
            # print(f"结束路径{new_path}")
        else:
            # path1的后半段和path2的前半段blend
            path1 = paths[i-1]
            path2 = paths[i]
            new_path = path_blend(path1,path2)
            traj.append(new_path)
            # print(f"融合路径{new_path}")
    return traj

# 路径展示比较   
def path_plot_compare(traj1,traj2,pt1,pt2):
    ax = plt.axes(projection='3d')
    for i in range(len(traj1)):
        ax.plot(traj1[i][0],traj1[i][1],traj1[i][2],lw=1,c='b')
        ax.plot(traj2[i][0],traj2[i][1],traj2[i][2],lw=1,c='r')
    plt.plot(*zip(*pt1), c="red",linestyle="none", marker=".")
    plt.plot(*zip(*pt2), c="red",linestyle="none", marker="x")
    ax.set_xlabel('X',labelpad=5)
    ax.set_ylabel('Y',labelpad=5)
    ax.set_zlabel('Z',labelpad=5)
    ax.set_xlim3d([-10,10])
    ax.set_ylim3d([-10,10])
    ax.set_zlim3d([-10,10])   
    plt.show()

# 路径展示  
def path_plot_num(traj1,points):
    ax = plt.axes(projection='3d')
    for i in range(len(traj1)):
        ax.plot(traj1[i][0],traj1[i][1],traj1[i][2],lw=1,c='b')
    # print(f"终点={traj1[-1][0],traj1[-1][1],traj1[-1][2]}")
    plt.plot(*zip(*points), c="red",linestyle="none", marker=".")
    ax.set_xlabel('X',labelpad=5)
    ax.set_ylabel('Y',labelpad=5)
    ax.set_zlabel('Z',labelpad=5)
    ax.set_xlim3d([-100,100])
    ax.set_ylim3d([-100,100])
    ax.set_zlim3d([-100,100])   
    # ax.set_xlim3d([-10,10])
    # ax.set_ylim3d([-10,10])
    # ax.set_zlim3d([-10,10]) 
    plt.show()
    
# 混合路径生成
def hybird_path(points):
    # 应该有n个点，n-1条线段，首条和末条线段直接使用，其余线段需要blend
    paths = []
    if len(points)<3:
        return False
    
    # 角度小于90，用椭圆;大于等于90,用圆
    for i in range(len(points)-2):
        p1 = points[i]
        p2 = points[i+1]
        p3 = points[i+2]
        rad = find_angle_of_3_pt(p1,p3,p2)
        # print(f"夹角角度={rad*R2D}")
        # 两个点构成一个segment，获取segment,然后对n-1和n segment进行blend
        
        if abs(rad) <= np.pi/2 and abs(rad) != 0:
            # print(f"当前线段使用圆形轨迹")
            pc,R,a,b,c,d = CircleFun(p1,p2,p3)
            # path 包含两个segment
            path = circle_interpolate(pc,R,a,b,c,d,p1,p2,p3)
        # 如果三点共线，直线插补
        elif abs(abs(rad)-np.pi)<1e-06 or abs(rad)<1e-06:
            path = straight_interploate(p1,p2,p3)
        else:
            # print(f"当前线段使用椭圆轨迹")
            l1 = math.dist(p1,p2)
            l2 = math.dist(p3,p2)
            # 令p1和p3中相距p2最远的点作为椭圆的短轴端点
            if l1>l2:
                ps = p1
                pe = p3
                revert = False
                # print(f"不需要反转路径")
            else:
                # 因为默认顺序是p1p2p3，所以如果ps=p3的话，需要倒置一下path
                ps = p3
                pe = p1
                revert = True
                # print(f"需要反转路径")
            po = calc_ellipse_center(ps,p2,pe)
            # 因为可以确定使用椭圆的时候，∠p1p2p3小于90°，于是长短轴是固定的
            long = p2-po
            short = ps-po
            path = ellipsis_interpolate(po,long,short,ps,pe,p2)
            if revert is True:
                # print(f"反转路径")
                # print(path)
                path[0] = path[0][::-1]
                path[1] = path[1][::-1]
                path[2] = path[2][::-1]
        paths.append(path)

    traj = path_reform(paths,points)
    return traj
    
# 圆弧路径生成 
def circle_path(points):
    paths = []
    # 角度小于90，用椭圆;大于等于90,用圆
    for i in range(len(points)-2):
        p1 = points[i]
        p2 = points[i+1]
        p3 = points[i+2]
        # print(i)
        pc,R,a,b,c,d = CircleFun(p1,p2,p3)
        if i == len(points)-3:
            path = circle_interpolate(pc,R,a,b,c,d,p1,p2,p3)
            path_new = path
        else:
            path = circle_interpolate(pc,R,a,b,c,d,p1,p2,p3)
            x = path[0][:NUMPOINT]
            y = path[1][:NUMPOINT]
            z = path[2][:NUMPOINT]
            path_new = [x,y,z]
        paths.append(path_new)
        
    ax = plt.axes(projection='3d')
    for i in range(len(paths)):
        ax.plot(paths[i][0],paths[i][1],paths[i][2],lw=1,c='r')
    ax.set_xlabel('X',labelpad=5)
    ax.set_ylabel('Y',labelpad=5)
    ax.set_zlabel('Z',labelpad=5)
    ax.set_xlim3d([-10,10])
    ax.set_ylim3d([-10,10])
    ax.set_zlim3d([-10,10])
    plt.show()
    return paths

## 单元测试：
# 1.local support (done)
# 2.随机数据 (done)
# 3.三种混合插补方式测试 (done)
# 4.极端情况 (done)
# 5.连续性验证 (done,二阶连续)
# 6.曲率计算 (done,最大curvature出现在control point附近)
# 7.机器人结果验证(pyrob)

def local_support_test():
    p1 = np.array([-2,0,0]) 
    p2 = np.array([0,2,0])  
    p3 = np.array([2,0,0]) 
    p4 = np.array([5,-1,0]) 
    p5 = np.array([3,2,0])
    p6 = np.array([-4,4,0])
    ps = np.array([6,-1.5,0])
    points1 = [p1,p2,p3,p4,p5,p6]
    pathsh1 = hybird_path(points1)
    points2 = [p1,p2,p3,ps,p5,p6,p1]
    pathsh2 = hybird_path(points2)
    path_plot_compare(pathsh1,pathsh2,points1,points2)
    
def random_point_test():
    # 随机测试数据
    points = []
    points.append(np.array([0,0,0]))
    for i in range(10):
        x = random.randint(-50,50)
        y = random.randint(-50,50)
        z = random.randint(-50,50)
        p = np.array([x,y,z])
        points.append(p)
    points.append(np.array([0,0,0]))
    # print(f"最后一点={points[-1]}")
    paths = hybird_path(points)
    path_plot_num(paths,points)
    
def three_kind_interp_test():
    p1 = np.array([-2,0,0]) 
    p2 = np.array([0,2,0])  
    p3 = np.array([2,0,0]) 
    p4 = np.array([4,-1,0]) 
    # p4 = np.array([4,1,0])
    p5 = np.array([-3,-3,0]) # 探讨此点的问题，某些情况下数据混乱
    p6 = np.array([-1,-1,0])
    p7 = np.array([2,2,0])
    p8 = np.array([3,3,0])
    p9 = np.array([5,5,0])
    p10 = np.array([7,5,0])
    pts = [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10]
    paths = hybird_path(pts)
    path_plot_num(paths,pts)

def extreme_test():
    p1 = np.array([-2,0,0]) 
    p2 = np.array([0,2,0])  
    p3 = np.array([2,0,0]) 
    p4 = np.array([9,-12,0]) 
    p5 = np.array([10,-11,0])
    p6 = np.array([-15,15,0])
    p7 = np.array([-12,12,0])
    pts = [p1,p2,p3,p4,p5,p6,p7]
    paths = hybird_path(pts)
    path_plot_num(paths,pts)
    
def Curvature_test():
    points = []
    points.append(np.array([0,0,0]))
    for i in range(5):
        x = random.randint(-10,10)
        y = random.randint(-10,10)
        z = random.randint(-10,10)
        # z = 0
        p = np.array([x,y,z])
        points.append(p)
    paths = hybird_path(points)
    curve = []
    for i in range(len(paths)):
        dx_dt = np.gradient(paths[i][0])
        dy_dt = np.gradient(paths[i][1])
        dz_dt = np.gradient(paths[i][2])
        # curve
        d2x_t = np.gradient(dx_dt)
        d2y_t = np.gradient(dy_dt)
        d2z_t = np.gradient(dz_dt)
        # (a2b3-a3b2,a3b1-a1b3,a1b2-a2b1
        curvature_val = np.abs(d2y_t*dz_dt-d2z_t*dy_dt + d2z_t*dx_dt-d2x_t*dz_dt + d2x_t*dy_dt-d2y_t*dx_dt) / (dx_dt**2 + dy_dt**2 + dz_dt**2)**1.5
        for i in range(NUMPOINT):
            curve.append(curvature_val[i])
    l = len(curve)
    t = np.linspace(0,l,l)
    plt.plot(t,curve,'b-*')
    plt.title('Curvature')
    plt.show()
    path_plot_num(paths,points)
       
def Gradient_test():
    points = []
    points.append(np.array([0,0,0]))
    for i in range(5):
        x = random.randint(-10,10)
        y = random.randint(-10,10)
        z = random.randint(-10,10)
        p = np.array([x,y,z])
        points.append(p)
    paths = hybird_path(points)
    c2 = []
    for i in range(len(paths)):
        dx_dt = np.gradient(paths[i][0])
        dy_dt = np.gradient(paths[i][1])
        dz_dt = np.gradient(paths[i][2])
        ds_dt = np.sqrt(dx_dt**2+dy_dt**2+dz_dt**2)
        d2s_dt = np.gradient(ds_dt)
        c2.append(d2s_dt)
    C2 = []
    for i in range(len(c2)):
        for j in range(len(c2[i])):
            c = c2[i][j]
            C2.append(c)
    l = len(C2)
    t = np.linspace(0,l,l)
    plt.plot(t,C2,'b-*')
    plt.title('Continuity')
    plt.show()
    path_plot_num(paths,points)


        

if __name__ == "__main__":  
    # local_support_test()  
    # random_point_test()
    # three_kind_interp_test()
    # extreme_test()
    Curvature_test()
    # Gradient_test()
    # p1 = np.array([5.2,-1.5,2.6])
    # p2 = np.array([-2.3,-2.8,1.9])
    # p3 = np.array([0.7,-1.0,-0.3])
    # pc,R,a,b,c,d = CircleFun(p1,p2,p3)
    # print(f"pc={pc},R={R}")
    # po = calc_ellipse_center(p1,p2,p3)
    # print(f"po = {po}")
    pass