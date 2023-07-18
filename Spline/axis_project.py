# 计算基于世界坐标系的新坐标系的旋转矩阵
def calPoseFrom3Points(Oab, Pxb, Pyb):
    '''
    功能：已知坐标系a的原点、x轴正半轴上任一点和y轴正半轴上任一点在坐标系b下的坐标，
        求解坐标系a到坐标系b的旋转矩阵R和平移矩阵T；
    输入：坐标系a的原点在坐标系b下的坐标:Oab(x1,y1,z1);
        坐标系a的x轴正半轴上任一点在坐标系b下的坐标:Pxb(x2,y2,z2);
        坐标系a的y轴正半轴上任一点在坐标系b下的坐标:Pyb(x3,y3,z3);
    输出：坐标系n到坐标系s的旋转矩阵Rns，输出格式为矩阵;
    '''
    x = (Pxb - Oab) / np.linalg.norm(Pxb - Oab)
    y = (Pyb - Oab) / np.linalg.norm(Pyb - Oab)
    z = np.cross(x, y)
    length = np.linalg.norm(z)
    z = z / length
    Rab = np.matrix([x, y, z]).transpose()
    Tab = np.matrix(Oab).transpose()
    Mat = np.eye(4)
    for i in range(3):
        Mat[0][i] = x[i]
        Mat[1][i] = y[i]
        Mat[2][i] = z[i]
        Mat[i][3] = Oab[i]
    return Mat

# 将对应的两点投影到相对应平面上，从3d变2d，方便插补
def point_project2plane(plane: Plane, point: Point):
    # plane = Plane(point=[0, 0, 2], normal=[1, 0, 2])
    # point = Point([5, 9, 3])
    point_projected = plane.project_point(point)
    vector_projection = Vector.from_points(point, point_projected)

    plot_3d(
        plane.plotter(lims_x=(0, 10), lims_y=(0, 10), alpha=0.3),
        point.plotter(s=75, c='k'),
        point_projected.plotter(c='r', s=75, zorder=3),
        vector_projection.plotter(point=point, c='k', linestyle='--'),
    )
    plt.show()
    return point_projected