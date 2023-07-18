from sympy import Point3D, Plane
import math
import numpy as np
from matplotlib import pyplot as plt

D2R = np.pi/180
R2D = 180/np.pi
ITER = 32
# a = Plane(Point3D(1, 1, 2), Point3D(2, 4, 7), Point3D(3, 5, 1))
# print(a.equation())
# print(type(a.normal_vector))

# theta已知,a1=0.25*theta,a2= theta-a1
# theta = 71.565*D2R
# 第一组数据
p1 = np.array([2,0,1]) # pj
p2 = np.array([0,4,1])  # pi
p3 = np.array([-1,math.sqrt(12),1]) # pk

# 第二组数据
# p1 = np.array([4,0]) # pj
# p2 = np.array([0,2])  # pi
# p3 = np.array([-1,math.sqrt(15)/2])
# p3 = np.array([-2,math.sqrt(3)]) # pk
# p3 = np.array([-4,0])

# 第三组数据
# p1 = np.array([3,0]) # pj
# p2 = np.array([0,2])  # pi
# # p3 = np.array([-2,math.sqrt(20)/3]) # pk
# p3 = np.array([-math.sqrt(27)/2,1])

# 第四组数据
# p1 = np.array([-2,0]) # pj
# p2 = np.array([0,3])  # pi
# # p3 = np.array([-1,math.sqrt(27)/2])
# p3 = np.array([math.sqrt(20)/3,2])

# i = 1,2,3...
vec1 = p1-p2 # pi-pj, j=i-1+points.length
vec2 = p3-p2 # pi-pk, k=i+1
print(f"vecs={vec1,vec2}")
# len1 = math.sqrt(vec1[0]**2+vec1[1]**2)
# len2 = math.sqrt(vec2[0]**2+vec2[1]**2)
len1 = math.sqrt(vec1[0]**2+vec1[1]**2+vec1[2]**2)
len2 = math.sqrt(vec2[0]**2+vec2[1]**2+vec2[2]**2)
# print(f"lens={len1,len2}")
# cos_angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
cos_angle = (vec1[0]*vec2[0]+vec1[1]*vec2[1])/(len1*len2)
theta = np.arccos(cos_angle)
print(f"theta_between_p1p2&p3p2={theta*R2D}")

y = []
x = []
angle = []
ang = 0.5*theta
incA = 0.25*theta
inca = []

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
        
    y.append(v)
    x.append(i)
    # print(f"theta={alpha*R2D}")
    angle.append(alpha*R2D)
    inca.append(incA*R2D)
    # print(f"angle={ang*R2D}")
    

# print(f"a,b,c,d={a,b,c,d}")
# print(f"ac={(c/a)**2},bd={(1-d/b)**2}")
print(f"theta1={angle[-1]}")
plt.plot(x,y,'r-')
plt.plot(x,angle,'b-')
plt.plot(x,inca,'g-')

# print(f"min_v={min(y)}")
# print(f"final_ang={angle[-1]}")

# if len1<len2:
#     vec = vec2
#     l = len2 # l is len
#     pt2 = p3 # pk->p3
# else:
#     vec = vec1
#     l = len1 
#     pt2 = p1 # pj->p1
    
# dir = vec/l
# #print(f"dir={dir}")
# cross = np.cross(vec1,vec2)
# #print(f"cross={cross}")
# if (len1<len2 and cross>0) or (len2>len1 and cross<0):
#     perp = np.array([dir[1],-dir[0]])
# else:
#     perp = np.array([-dir[1],dir[0]])
# #print(f"perp={perp}")
# v = b*b/l
# h = b*a/l
# # print(f"v,h={v,h}")
# # axis1 = np.dot(h,perp)-np.dot(-v,dir)
# # center = axis1-p1
# axis1 = np.dot(-v,dir)-np.dot(h,perp)
# center = p1-axis1
# print(f"axis1={axis1}")
# print(f"center={center}")
plt.show()

