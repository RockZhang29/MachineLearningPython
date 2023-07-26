import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt
from numpy.linalg import det

from typing import Iterable

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
        # print(f"ti={ti}")
        for i in ti:
            if isinstance(i, complex) and i.imag < 0.01:
                i = i.real
                if i > 0 and i < 1:
                    ti = i
                    break
        if(isinstance(ti, float) == False):
            raise Exception(f"Error compute t_i:{ti}")

        self.b1 = (p1-(1-ti)**2*p0-ti**2*p2)/(2*(1-ti)*ti)

    def __call__(self, t: float):
        p0 = self.pts[0]
        p2 = self.pts[2]
        if is_real(t):
            # print(f"t={t}")
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

    
def testNd():
    x, y, z = np.array([[0, 0, 0], [-1, 0, 2], [2, 0, 1]])
    interp = QuadraticBezier(x, y, z)
    X = np.linspace(0, 1, 100)
    Y = interp(X)
    ax = plt.axes(projection='3d')
    ax.plot(Y[:, 0],Y[:, 1],Y[:, 2],'b--')
    plt.show()
    
def testQ():
    pts = np.array([[0,0,0,0,0,0], [4,5,2,0,90,75], [3,1,-1,42,-12,22], [4,0,2,-14,37,29]])
    # plt.scatter(pts[:, 0], pts[:, 1], pts[:,2])
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
        
    ax = plt.axes(projection='3d')
    tx,ty,tz = [],[],[]
    for i in range(len(x)):
        for j in range(len(x[i])):
            tx.append(x[i][j])
            ty.append(y[i][j])
            tz.append(z[i][j])    
    ax.plot(tx,ty,tz,'b-.')
    plt.show()
    
    jx = plt.axes(projection='3d')
    ja,jb,jc = [],[],[]
    for i in range(len(j1)):
        for j in range(len(j1[i])):
            ja.append(j1[i][j])
            jb.append(j2[i][j])
            jc.append(j3[i][j])    
    jx.plot(ja,jb,jc,'r--')
    plt.show()

    
    
if __name__ == "__main__":  
    # testNd()
    testQ()