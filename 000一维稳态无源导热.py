import numpy as np

# 域离散
L = 0.5       # 棒总长
N = 5       # 控制单元数目
dx = L / N    # 控制单元长度
centriod = [i * dx + dx / 2 for i in range(N)]   # 控制单元中心坐标
A = [0.001] * N                                  # 控制单元中心对应的界面面积
k = [1000] * N                                   # 控制体单元中心处的导热系数
T = [0] * N                                      # 控制体单元中心温度初始值。可以任意给定
A = np.zeros((N, N))                             # 整体系数矩阵
b = np.zeros(N)                                  # 整体方程右端项
TA = 100
TB = 500


# 方程离散和整体组装
for ii in range(1, N-1):  # 用来组装整体系数矩阵的整体循环，注意边界单元没有循环
    A[ii, ii-1] = -(k[ii-1] + k[ii]) / 2 / dx
    A[ii, ii+1] = -(k[ii]+k[ii+1]) / 2 / dx
    A[ii, ii] = -(A[ii, ii-1] + A[ii, ii+1])
    b_C = 0
    b[ii] = b_C

# 第一个边界单元修正
A[0, 1] = -(k[0] + k[1]) / 2 / dx
b[0] = 2 * k[0] / dx * TA
A[0, 0] = -(A[0, 1]) + b[0] / TA

# 第二个边界单元修正
A[N-1, N-2] = -(k[N-2]+k[N-1]) / 2 / dx
b[N-1] = k[N-1] * 2 / dx * TB
A[N-1, N-1] = -(A[N-1, N-2]) + b[N-1] / TB

# 方程求解
T = np.linalg.solve(A, b)

# 结果可视化
import matplotlib.pyplot as plt
fig = plt.figure()
plt.plot(np.arange(N)/N*L, T, 'o')
plt.show()


# 结果分析
import matplotlib.pyplot as plt
x = np.array([0, *centriod, L])
T_all = np.array([TA, *T, TB])
exact = 800 * x + 100
fig = plt.figure()
plt.plot(x, T_all, 'o', label='FVM')
plt.plot(x, exact, '-', label='Exact')
plt.legend()
plt.show()
