import numpy as np
import matplotlib.pyplot as plt

# 数据导入
data = np.loadtxt('newoutput1.csv', delimiter=',')
rows, cols = data.shape

# 处理数据的矩阵重组
if cols % 10 == 0:
    Tck = np.zeros((rows, (cols // 10) * 3))
    indices = np.arange(0, cols, 10)
    for i, idx in enumerate(indices):
        Tck[:, 3 * i: 3 * i + 3] = data[:, idx: idx + 3]
else:
    Tck = np.zeros((rows, cols))

# 初始化结果矩阵
data_new = np.zeros((10000, 150))

# Sim函数向量化
def v_int_Sim_vectorized(t1, t2, row, data):
    h = 0.0005
    idx = np.arange(t1 - 1, t2 - 2)
    dis = np.sum((data[idx, row] + 4 * data[idx + 1, row] + data[idx + 2, row]) * (h / 6))
    return dis

# Tra函数向量化
def v_int_tra_vectorized(t1, t2, row, data):
    h = 0.0005
    idx = np.arange(t1 - 1, t2 - 1)
    dis = np.sum((data[idx, row] + data[idx + 1, row]) * (h / 2))
    return dis

# 填充 data_new，向量化优化
for i in range(150):
    data_new[3:10000, i] = [v_int_Sim_vectorized(1, j, i, Tck) for j in range(3, 10000)]
    data_new[2, i] = v_int_tra_vectorized(1, 2, i, Tck)

# 时间轴生成
t = np.arange(0.0005, 5.0005, 0.0005)

# 提取数据绘图
Y1, Y50 = data_new[:, 0], data_new[:, 147]
X1, X50 = data_new[:, 1], data_new[:, 148]
Z1, Z50 = data_new[:, 2], data_new[:, 149]

# 统一绘图参数
def plot_displacement(t, data, title, subplot_pos):
    plt.subplot(2, 3, subplot_pos)
    plt.plot(t, data, '-', linewidth=1.8)
    plt.xlim([0, 10])
    plt.ylim([-0.5, 0.5])
    plt.xticks(np.arange(0, 11, 1))
    plt.yticks(np.arange(-0.5, 0.6, 0.1))
    plt.title(title, fontsize=16, fontname='Times New Roman')
    plt.xlabel("Time(s)", fontname='Times New Roman')
    plt.ylabel("Displacement(m)", fontname='Times New Roman')
    plt.grid(True)

# 绘图
plt.figure(figsize=(12, 6))

plot_displacement(t, Y1, "Displacement at Point1 in Y-direction", 1)
plot_displacement(t, Y50, "Displacement at Point50 in Y-direction", 4)
plot_displacement(t, X1, "Displacement at Point1 in X-direction", 2)
plot_displacement(t, X50, "Displacement at Point50 in X-direction", 5)
plot_displacement(t, Z1, "Displacement at Point1 in Z-direction", 3)
plot_displacement(t, Z50, "Displacement at Point50 in Z-direction", 6)

plt.tight_layout()
plt.show()
