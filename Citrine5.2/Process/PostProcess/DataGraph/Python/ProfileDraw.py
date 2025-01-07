import numpy as np
import matplotlib.pyplot as plt

# 数据导入
data = np.loadtxt('newoutput1.csv', delimiter=',')
rows, cols = data.shape

# 检查列数是否为10的倍数
if cols % 10 == 0:
    Tck = np.zeros((rows, (cols // 10) * 3))
    indices = np.arange(0, cols, 10)
    for i, idx in enumerate(indices):
        Tck[:, 3 * i: 3 * i + 3] = data[:, idx:idx + 3]
else:
    Tck = data.copy()

# 新数据矩阵初始化
data_new = np.zeros((10000, 150))

# 计算 v_int_Sim 和 v_int_tra 的向量化版本
def v_int_Sim_vectorized(t1, t2, row, data):
    h = 0.0001
    dis = np.sum((data[t1:t2-2, row] + 4 * data[t1+1:t2-1, row] + data[t1+2:t2, row]) * (h / 6))
    return dis

def v_int_tra_vectorized(t1, t2, row, data):
    h = 0.0001
    dis = np.sum((data[t1:t2-1, row] + data[t1+1:t2, row]) * (h / 2))
    return dis

# 向量化填充 data_new
rows = np.arange(3, 10000)[:, None]
cols = np.arange(150)

for i in cols:
    data_new[3:10000, i] = np.array([v_int_Sim_vectorized(1, j, i, Tck) for j in range(3, 10000)])
    data_new[2, i] = v_int_tra_vectorized(1, 2, i, Tck)

# 修改 data_new 的第二列值，使用广播提升速度
K = np.arange(0, -10.2, -0.2)
j_indices = np.arange((data_new.shape[1] // 3) - 1)  # 避免超界
data_new[:, 3 * j_indices + 2] += K[j_indices]

# 定义点的集合，向量化
Xr = data_new[np.arange(0, 2001, 100)[:, None], 3 * j_indices + 2]
Yr = data_new[np.arange(0, 2001, 100)[:, None], 3 * j_indices + 3]
Zr = data_new[np.arange(0, 2001, 100)[:, None], 3 * j_indices + 1]

# 定义列集合
Xc, Yc, Zc = data_new[:1000, 2], data_new[:1000, 3], data_new[:1000, 1]

# 绘制3D图
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for i in range(Xr.shape[0]):
    ax.plot3D(Xr[i], Yr[i], Zr[i])

# 绘制中心点
ax.plot3D(Xc, Yc, Zc, linestyle=":", color="#A2142F", linewidth=1.2)

# 图形设置
ax.set_xlim([-12, 1])
ax.set_ylim([-1, 1])
ax.set_zlim([-1, 1])
ax.set_xticks(range(-50, 2, 10))
ax.set_yticks(range(-25, 26, 10))
ax.set_zticks(range(-25, 26, 10))
ax.set_title("Profile in Three-dimensional space at Anytime", fontsize=16, fontname='Times New Roman')
ax.set_xlabel('Z', fontname='Times New Roman')
ax.set_ylabel('Y', fontname='Times New Roman')
ax.set_zlabel('X', fontname='Times New Roman')
ax.view_init(60, 30)
ax.grid(True)
plt.show()
