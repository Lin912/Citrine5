clc;clear;
syms Theta Phi Vx Vy Vz u v w;
syms Omega2 Omega3;
syms T Sn Sb;

Trans = [cos(Phi)*cos(Theta),   cos(Theta)*sin(Phi),    -sin(Theta);
         -sin(Phi),             cos(Phi),               0; 
         cos(Phi)*sin(Theta),   sin(Phi)*sin(Theta),    cos(Theta)];
%Trans from X-Y-Z -> t-n-b

V0 = [Vx;Vy;Vz];
u0 = [u;v;w];

ff = Trans* V0;
f = u0-ff;
f(4) = Omega2;
f(5) = Omega3;


jac = jacobian(f,[u,v,w,T,Sn,Sb,Theta,Phi,Omega2,Omega3]);