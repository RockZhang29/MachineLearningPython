import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# 设置四个控制点，使用三阶贝塞尔曲线
P0,P1,P2,P3 = ((0,0,0),(100,100,150),(30,80,90),(200,150,200))
t = np.arange(0,1,0.01)
a0 = (1-t)**3
a1 = 3*(1-t)**2*t
a2 = 3*t**2*(1-t)
a3 = t**3

x = a0*P0[0]+a1*P1[0]+a2*P2[0]+a3*P3[0]
y = a0*P0[1]+a1*P1[1]+a2*P2[1]+a3*P3[1]
z = a0*P0[2]+a1*P1[2]+a2*P2[2]+a3*P3[2]

ax = plt.axes(projection='3d')
ax.plot(x,y,z,lw=2,c='b')
ax.set_xlabel('X',labelpad=5)
ax.set_ylabel('Y',labelpad=5)
ax.set_zlabel('Z',labelpad=5)
ax.set_title('3D Bezier Curve')
plt.plot(*P0,c="red",linestyle="none",marker="o")
plt.plot(*P1,c="red",linestyle="none",marker="o")
plt.plot(*P2,c="red",linestyle="none",marker="o")
plt.plot(*P3,c="red",linestyle="none",marker="o")
plt.tight_layout()
plt.show()

