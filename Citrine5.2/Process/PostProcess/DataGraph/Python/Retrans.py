import numpy as np
import pandas as pd

def process_data(data, transform_cols):
    """
    """
    for i in range(data.shape[0]):
        for j in range(data.shape[1] // 10):
            start_idx = 10 * j
            rotation_matrix = np.linalg.inv(np.array([
                [np.cos(data[i, start_idx + 7]) * np.cos(data[i, start_idx + 6]), -np.sin(data[i, start_idx + 7]), np.sin(data[i, start_idx + 6]) * np.cos(data[i, start_idx + 7])],
                [np.cos(data[i, start_idx + 6]) * np.sin(data[i, start_idx + 7]),  np.cos(data[i, start_idx + 7]), np.sin(data[i, start_idx + 7]) * np.sin(data[i, start_idx + 6])],
                [-np.sin(data[i, start_idx + 6]), 0, np.cos(data[i, start_idx + 6])]
            ]))
            data[i, start_idx + transform_cols[0]:start_idx + transform_cols[1]] = \
                np.dot(data[i, start_idx + transform_cols[0]:start_idx + transform_cols[1]], rotation_matrix)
    return data

# 主程序
if __name__ == "__main__":
    # 加载数据
    data = pd.read_csv('output.csv', header=None).fillna(0).values
    
    # 检查列数是否为10的倍数
    if data.shape[1] % 10 != 0:
        raise ValueError("数据列数必须是10的倍数。")

    # 处理前三列逻辑（对应 Retans.m）
    data1 = process_data(data.copy(), transform_cols=(0, 3))
    np.savetxt('newoutput1.csv', data1, delimiter=',')

    # 处理第4-6列逻辑（对应 Retans2.m）
    data2 = process_data(data.copy(), transform_cols=(3, 6))
    np.savetxt('newoutput2.csv', data2, delimiter=',')
