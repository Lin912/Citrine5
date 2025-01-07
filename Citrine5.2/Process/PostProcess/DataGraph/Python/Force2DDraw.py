import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

# 加载数据
data = pd.read_csv('newoutput2.csv', header=None).values

# 初始化 Fx, Fy, Fz 矩阵
num_rows, num_cols = data.shape
Fx = np.zeros((num_rows, num_cols // 10))
Fy = np.zeros((num_rows, num_cols // 10))
Fz = np.zeros((num_rows, num_cols // 10))

# 提取全部 Force
if num_cols % 10 == 0:
    flag = 1
else:
    flag = 0

if flag == 1:
    for j in range(1, num_cols // 10 + 1):
        for i in range(num_rows):
            F = data[i, 10 * (j - 1) + 3:10 * (j - 1) + 6]
            Fx[i, j - 1] = F[0]
            Fy[i, j - 1] = F[1]
            Fz[i, j - 1] = F[2]

# 时间数组
t = np.arange(0.0001, 1.0001, 0.0001)

# 设置 Matplotlib 使用嵌入字体
mpl.rcParams['pdf.fonttype'] = 42  # 嵌入字体，而不是转换为路径
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Times New Roman'


# 绘图
fig, axs = plt.subplots(3, 3, figsize=(15, 8))

# 定义绘图函数
def plot_subplot(ax, x, y, title, xlabel, ylabel, xlim, ylim, xticks, yticks, color='b'):
    ax.plot(x, y, '-', linewidth=1.5, color=color)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_title(title, fontsize=16, fontname='Times New Roman')
    ax.set_xlabel(xlabel, fontname='Times New Roman')
    ax.set_ylabel(ylabel, fontname='Times New Roman')
    ax.grid(True)

# 绘制每个子图
plot_subplot(axs[0, 0], t, Fy[:, 0], "Tension at Point1", "Time(s)", "Force(N)", [0, 10], [0, 150], range(0, 11, 1), range(0, 151, 15), color='#56B4E9')
plot_subplot(axs[1, 0], t, Fx[:, 0], "Fy at Point1", "Time(s)", "Force(N)", [0, 10], [-15, 15], range(0, 11, 1), range(-15, 16, 3))
plot_subplot(axs[2, 0], t, Fz[:, 0], "Fz at Point1", "Time(s)", "Force(N)", [0, 10], [-15, 15], range(0, 11, 1), range(-15, 16, 3))

plot_subplot(axs[0, 1], t, Fy[:, 24], "Tension at Point25", "Time(s)", "Force(N)", [0, 10], [0, 150], range(0, 11, 1), range(0, 151, 15))
plot_subplot(axs[1, 1], t, Fx[:, 24], "Fy at Point25", "Time(s)", "Force(N)", [0, 10], [-15, 15], range(0, 11, 1), range(-15, 16, 3))
plot_subplot(axs[2, 1], t, Fz[:, 24], "Fz at Point25", "Time(s)", "Force(N)", [0, 10], [-15, 15], range(0, 11, 1), range(-15, 16, 3))

plot_subplot(axs[0, 2], t, Fy[:, 49], "Tension at Point50", "Time(s)", "Force(N)", [0, 10], [0, 150], range(0, 11, 1), range(0, 151, 15))
plot_subplot(axs[1, 2], t, Fx[:, 49], "Fy at Point50", "Time(s)", "Force(N)", [0, 10], [-15, 15], range(0, 11, 1), range(-15, 16, 3))
plot_subplot(axs[2, 2], t, Fz[:, 49], "Fz at Point50", "Time(s)", "Force(N)", [0, 10], [-15, 15], range(0, 11, 1), range(-15, 16, 3))

plt.tight_layout()
plt.show()

# 保存为 PDF 或 EPS
# plt.savefig("figure.pdf", format="pdf")
