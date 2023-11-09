import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt
from numpy.linalg import det

from typing import Iterable

D2R = np.pi/180
R2D = 180/np.pi

def is_real(x):
    return isinstance(x, float) or isinstance(x, int)

class QuadraticBezier(object):
    """
    quadratic Bezier interpolation
    Usage:
        interp = QuadraticBezier(x, y, z)
        interp(t) # 0<=t<=1
    Example:
        x, y, z = np.array([[0, 0, 0], [1, 1, 2], [2, 3, 1]])
        interp = QuadraticBezier(x, y, z)
        X = np.linspace(0, 1, 100)
        Y = interp(X)
        plt.plot(Y[:, 0], Y[:, 1],Y[:, 2])
        plt.show()
    """
    def __init__(self, p0, p1, p2):
        self.pts = np.asarray((p0, p1, p2))
        # print(f"pts = {p0,p1,p2}")
        ti = np.roots([np.linalg.norm(p2-p0)**2,3*np.dot(p2-p0, p0-p1),np.dot(3*p0-2*p1-p2, p0-p1),-np.linalg.norm(p0-p1)**2])
        # ti = np.roots([a,b,c,d])
        a = np.linalg.norm(p2-p0)**2
        b = 3*np.dot(p2-p0, p0-p1)
        # print(f"p2-p0={p2-p0}")
        # print(f"p0-p1={p0-p1}")
        c = np.dot(3*p0-2*p1-p2, p0-p1)
        d = -np.linalg.norm(p0-p1)**2
        # print(f"para = {a,b,c,d}")
        # print(f"ti_before={ti}")
        for i in ti:
            # print(f"i in ti:{i}")
            if isinstance(i, complex) and i.imag < 0.01:
                # print(f"i={i}")
                i = i.real
                if i > 0 and i < 1:
                    ti = i
                    break
        if(isinstance(ti, float) == False):
            raise Exception(f"Error compute t_i:{ti}")
        
        #print(f"ti_after={ti}")
        self.b1 = (p1-(1-ti)**2*p0-ti**2*p2)/(2*(1-ti)*ti)
        # print(f"b1={self.b1}")

    def __call__(self, t: float):
        p0 = self.pts[0]
        p2 = self.pts[2]
        if is_real(t):
            # print(f"t={t}")
            p = (1-t)**2 * p0+2*(1-t)*t * self.b1 + t**2*p2
            print(f"path point = {p}")
            return (1-t)**2 * p0+2*(1-t)*t * self.b1 + t**2*p2
        elif isinstance(t, Iterable):
            return np.asarray([self(i) for i in t])
        else:
            raise NotImplementedError("type error")

class YukselC2Interpolation(object):
    """implementation of SIGGRAPH 2020 paper 'A class of C2 interpolating Splines - cem yuksel'"""
    def __init__(self, points, circle_join=False, interp=QuadraticBezier):
        if(len(points) < 3):
            raise Exception("Too few points")

        self.circle_join = circle_join
        self.interp = interp
        self.points = points
        # num_points-1 segments, num_points-2 interp_func
        self.num_points = len(self.points)
        self.interp_func = [interp(*self.points[i:i+3]) for i in range(len(self.points)-2)]
        if circle_join:
            self.interp_func.append(interp(points[-2], points[-1], points[1]))

    def call(self, n: int, t: float):
        assert (0 <= n and ((not self.circle_join and n < self.num_points-2) or (self.circle_join and n < self.num_points-1))), \
            f"Has only {self.num_points-2} interpolation segments"
        return self.interp_func[n](t)

    def scaled_call(self, n: int, theta: float):
        if isinstance(theta, Iterable):
            t = np.asarray(theta)/np.pi
        elif is_real(theta):
            t = theta / np.pi
        else:
            raise Exception("Type Error")

        return self.call(n, t)

    def blend(self, n: int, theta: float):
        """blend n-th and (n+1)-th segment, 0<theta<np.pi/2"""
        assert (0 <= n < self.num_points-1), f"Has only {self.num_points-1} blended segments"

        def _blend(n1, n2):
            if is_real(theta):
                # print(f"theta={theta}")
                # if(theta == np.pi/4):
                #     print(f"theta={theta},points={np.cos(theta)**2*self.scaled_call(n1, theta+np.pi/2) + np.sin(theta)**2 * self.scaled_call(n2, theta)}")
                return np.cos(theta)**2*self.scaled_call(n1, theta+np.pi/2) + np.sin(theta)**2 * self.scaled_call(n2, theta)
            elif isinstance(theta, Iterable):
                return np.asarray([self.blend(n, i) for i in theta])
            else:
                Exception("type Error")
        
        if (n == 0):
            if not self.circle_join:
                return self.scaled_call(0, theta)
            else:
                return _blend(self.num_points-2, 0)
        if (not self.circle_join) and (n == self.num_points - 2):
            return self.scaled_call(n-1, theta+np.pi/2)

        return _blend(n-1, n)

def Nd_Gradient_test(points):
    # x,y,z,j4,j5,j6
    j1,j2,j3,j4,j5,j6 = [],[],[],[],[],[]
    for i in range(len(points)):
        j1.append(points[i][0])
        j2.append(points[i][1])
        j3.append(points[i][2])
        j4.append(points[i][3])
        j5.append(points[i][4])
        j6.append(points[i][5])
    dx_dt = np.gradient(j4)
    dy_dt = np.gradient(j5)
    dz_dt = np.gradient(j6)
    d2x_t = np.gradient(dx_dt)
    d2y_t = np.gradient(dy_dt)
    d2z_t = np.gradient(dz_dt)
    ds = np.sqrt(dx_dt**2 + dy_dt**2 + dz_dt**2)
    dx = dx_dt/ds
    dy = dy_dt/ds
    dz = dz_dt/ds
    curvature_val = np.abs(d2y_t*dz_dt-d2z_t*dy_dt + d2z_t*dx_dt-d2x_t*dz_dt + d2x_t*dy_dt-d2y_t*dx_dt) / (dx_dt**2 + dy_dt**2 + dz_dt**2)**1.5

    d1_dt = np.gradient(j1)
    d2_dt = np.gradient(j2)
    d3_dt = np.gradient(j3)
    ds = np.sqrt(d1_dt**2 + d2_dt**2 + d3_dt**2)
    d1 = d1_dt/ds
    d2 = d2_dt/ds
    d3 = d3_dt/ds
    d1_2t = np.gradient(d1_dt)
    d2_2t = np.gradient(d2_dt)
    d3_2t = np.gradient(d3_dt)
    acc_val = np.abs(d2_2t*d3_dt-d3_2t*d2_dt + d3_2t*d1_dt-d1_2t*d3_dt + d1_2t*d2_dt-d2_2t*d1_dt) / (d1_dt**2 + d2_dt**2 + d3_dt**2)**1.5
    
    ## 加速度
    plt.plot(np.cumsum(ds),d1,'y-*',label='x')
    plt.plot(np.cumsum(ds),d2,'k-*',label='y')
    plt.plot(np.cumsum(ds),d3,'m-*',label='z')
    plt.plot(np.cumsum(ds),dx,'r-*',label='joint1')
    plt.plot(np.cumsum(ds),dy,'b-*',label='joint2')
    plt.plot(np.cumsum(ds),dz,'g-*',label='joint3')
    plt.title('6D Bezier Gradient')
    plt.legend()
    plt.show()
    
def testQ():
    PI = np.pi
    pts = np.array([[1,2,3,1,-2,1.2], [0,1,2,-1.2,2.3,0], [-2,2.5,-1.5,0.4,-1.5,3]])
    # pts = np.array([[0,0,0],[4,2,5],[-1,3,5],[5,-4,2]])
    interp = YukselC2Interpolation(pts)
    X = np.linspace(0, np.pi/2, 100)

    x,y,z = [],[],[]
    j1,j2,j3 = [],[],[]
    for i in range(len(pts)-1):
        Y = interp.blend(i, X)
        x.append(Y[:,0])
        y.append(Y[:,1])
        z.append(Y[:,2])
        j1.append(Y[:,3])
        j2.append(Y[:,4])
        j3.append(Y[:,5])
        
    # ax = plt.axes(projection='3d')
    tx,ty,tz = [],[],[]
    for i in range(len(x)):
        for j in range(len(x[i])):
            tx.append(x[i][j])
            ty.append(y[i][j])
            tz.append(z[i][j])    
    # ax.plot(tx,ty,tz,'b-.')
    # plt.show()
    
    # jx = plt.axes(projection='3d')
    ja,jb,jc = [],[],[]
    for i in range(len(j1)):
        for j in range(len(j1[i])):
            ja.append(j1[i][j])
            jb.append(j2[i][j])
            jc.append(j3[i][j])    
    # jx.plot(ja,jb,jc,'r--')
    # plt.show()

    # 将两个3d->一个6d
    points = []
    for i in range(len(tx)):
        point = np.array([tx[i],ty[i],tz[i],ja[i],jb[i],jc[i]])
        points.append(point)
    # print(points)
    Nd_Gradient_test(points)

def test_cubic_root():
    a = np.float64(1)
    b = np.float64(1)
    c = np.float64(1)
    d = np.float64(-3)
    ti = np.roots([a,b,c,d])
    
    print(f"ti={ti}")

def normal_test():
    PI = np.pi
    pts = np.array([[0,0,1],[1,3,2],[-2,-1,-1],[5,-4,0]])
    # plt.scatter(pts[:, 0], pts[:, 1])
    interp = YukselC2Interpolation(pts, circle_join=False)
    X = np.linspace(0, np.pi/2, 100)
    ax = plt.axes(projection='3d')
    x,y,z=[],[],[]
    for i in range(len(pts)-1):
        Y = interp.blend(i, X)
        ax.plot(Y[:,0], Y[:,1], Y[:,2])
        x.append(Y[:,0])
        y.append(Y[:,1])
        z.append(Y[:,2])
    # print(f"len(X)={len(x[0])}")
    # print(f"x中间点:{Y[:,0]}")
    ax.scatter(pts[:,0],pts[:,1],pts[:,2]) 
    plt.show()

if __name__ == "__main__":  
    # testQ()
    # test_cubic_root()
    normal_test()