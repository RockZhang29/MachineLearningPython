import numpy as np
import matplotlib.pyplot as plt
## 标准方程的实现
# 随机生成的线性数据集,函数为y=4+3x+高斯噪声
x = 2*np.random.rand(100,1)
y = 4+3*x+np.random.rand(100,1)
X_b = np.c_[np.ones((100,1)), x] # 按行连接两个矩阵
theta_best = np.linalg.inv(X_b.T.dot(X_b)).dot(X_b.T).dot(y)
print(f"实际结果为：{theta_best[0],theta_best[1]},期待的theta0和theta1分别为{4,3}")
# 绘制预测模型
X_new = np.array([[0],[2]])
X_new_b = np.c_[np.ones((2,1)),X_new]
y_predict = X_new_b.dot(theta_best)

plt.plot(X_new,y_predict, "r-")
plt.plot(x,y,"b.")
plt.axis([0,2,0,15])
plt.show()

## 梯度下降（Gradient Decent）
eta = 0.1 # learning rate, 学习率
iterations = 1000
m = 100 # 训练实例数量
theta = np.random.randn(2,1)
for i in range(iterations):
    gradients = 2/m*X_b.T.dot(X_b.dot(theta)-y)
    theta = theta-eta*gradients
    
# SGD
n_epochs = 50
t0,t1 = 5,50
def learning_schedule(t):
    return t0/(t+t1)
theta = np.random.randn(2,1)
for j in range(n_epochs):
    for i in range(m):
        random_index = np.random.randint(m)
        xi = X_b[random_index:random_index+1]
        yi = y[random_index:random_index+1] 
        gradients - 2*xi.T.dot(xi.dot(theta)-yi)
        eta = learning_schedule(j*m+i)
        theta = theta-eta*gradients

## 多项式线性回归,拟合复杂曲线模型

