import numpy as np
import matplotlib.pyplot as plt
QUADRUPLE_SIZE:int = 4
def num_segments(point_chain:tuple) -> int:
    return len(point_chain)-(QUADRUPLE_SIZE-1)

def flatten(list_of_lists) -> list:
    return [elem for lst in list_of_lists for elem in lst]

def catmull_rom_spline(P0: tuple, P1: tuple, P2: tuple, P3: tuple, num_points: int, alpha: float = 0.5):
    def tj(ti: float, pi: tuple, pj: tuple):
        xi,yi,zi = pi
        xj,yj,zj = pj
        dx,dy,dz = xj-xi,yj-yi,zj-zi
        l = (dx**2+dy**2+dz**2)**0.5
        return ti+l**alpha
    t0: float = 0.0
    t1: float = tj(t0,P0,P1)
    t2: float = tj(t1,P1,P2)
    t3: float = tj(t2,P2,P3)
    t = np.linspace(t1,t2,num_points).reshape(num_points,1)
    
    A1 = (t1-t) / (t1-t0) * P0 + (t-t0) / (t1-t0) * P1
    A2 = (t2-t) / (t2-t1) * P1 + (t-t1) / (t2-t1) * P2
    A3 = (t3-t) / (t3-t2) * P2 + (t-t2) / (t3-t2) * P3
    B1 = (t2-t) / (t2-t0) * A1 + (t-t0) / (t2-t0) * A2
    B2 = (t3-t) / (t3-t1) * A2 + (t-t1) / (t3-t1) * A3
    points = (t2-t) / (t2-t1) * B1 + (t-t1) / (t2-t1) * B2
    return points

def catmull_rom_chain(points: tuple, num_points: int) -> list:
    point_quadruples = (
        (points[idx_segment_start + d] for d in range(QUADRUPLE_SIZE))
        for idx_segment_start in range(num_segments(points))
    )
    all_splines = (catmull_rom_spline(*pq, num_points) for pq in point_quadruples)
    return flatten(all_splines)

if __name__ == "__main__":
    POINTS: tuple = ((-20,80,0),(0,100,0),(40,120,0),(30,75,0),(50,30,0),(65,70,0),(90,110,0),(120,70,0),(140,25,0),(165,25,0),(180,25,0))  # Red points
    POINTS1: tuple = ((-20,80,0),(0,100,0),(40,120,0),(30,75,0),(45,50,0),(65,70,0),(90,110,0),(120,70,0),(140,25,0),(165,25,0),(180,25,0))  # Red points
    NUM_POINTS: int = 200  # Density of blue chain points between two red points

    chain_points: list = catmull_rom_chain(POINTS, NUM_POINTS)
    assert len(chain_points) == num_segments(POINTS) * NUM_POINTS  
    chain_points1: list = catmull_rom_chain(POINTS1, NUM_POINTS)
    assert len(chain_points1) == num_segments(POINTS1) * NUM_POINTS 

    x1,y1,z1,x2,y2,z2 = [],[],[],[],[],[]
    for p1 in chain_points:
        x1.append(p1[0])
        y1.append(p1[1])
        z1.append(p1[2])
        
    for p2 in chain_points1:
        x2.append(p2[0])
        y2.append(p2[1])
        z2.append(p2[2])
    
    ax = plt.axes(projection='3d')
    ax.plot(x1,y1,z1,lw=1,c='b')
    ax.plot(x2,y2,z2,lw=1,c='r')
    ax.set_xlabel('X',labelpad=5)
    ax.set_ylabel('Y',labelpad=5)
    ax.set_zlabel('Z',labelpad=5)
    ax.set_title('3D Cutmull-Rom Curve')
    plt.plot(*zip(*chain_points), c="blue")
    plt.plot(*zip(*POINTS), c="yellow", linestyle="none", marker="o")
    plt.plot(*zip(*POINTS1), c="orange", linestyle="none", marker="x")
    plt.tight_layout()
    plt.show()