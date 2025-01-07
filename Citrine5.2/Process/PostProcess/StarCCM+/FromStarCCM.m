format long;

% 定义旋转角度
phi = -0.000101;    % 绕Z轴旋转
theta = 0.126918;   % 绕Y轴旋转
psi = 0.000014;     % 绕X轴旋转

% 构建旋转矩阵
Rz = [cos(phi), -sin(phi), 0;
      sin(phi), cos(phi), 0;
      0, 0, 1];

Ry = [cos(theta), 0, sin(theta);
      0, 1, 0;
      -sin(theta), 0, cos(theta)];

Rx = [1, 0, 0;
      0, cos(psi), -sin(psi);
      0, sin(psi), cos(psi)];

% 计算整体旋转矩阵 E
E = Rx * (Ry * Rz);

% 随体坐标系中的速度向量
v_body = [0.677722, 0.000059, 0.092987];

% 计算外部坐标系中的速度向量
v_inertial = E * v_body';

% 显示结果
disp('重心在外部坐标系下的速度为:');
disp(v_inertial);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 计算角速度向量的外部坐标系表示
omega_body = [0.000024; 0.103177; -0.000084];
r_body = [1; 0; 0];

% 使用前面计算好的旋转矩阵 E
omega_inertial = E * omega_body;
r_inertial = E * r_body;

% 计算叉积得到旋转引起的速度
v_rotation = cross(omega_inertial, r_inertial);

% 外部坐标系中重心的速度已计算为 v_inertial
% 非重心位置的总速度
v_non_cm_inertial = v_inertial + v_rotation;

% 显示结果
disp('非重心位置在外部坐标系下的速度为:');
disp(v_non_cm_inertial);


