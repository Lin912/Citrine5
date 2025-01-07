clc;clear;
syms Theta Phi u v w real;
syms Omega2 Omega3 real;
syms T Sn Sb real;
syms Gx Gy Gz real;
syms rho Cdn Cdb Sdn Sdb real;

Trans = [cos(Phi)*cos(Theta),   cos(Theta)*sin(Phi),    -sin(Theta);
         -sin(Phi),             cos(Phi),               0; 
         cos(Phi)*sin(Theta),   sin(Phi)*sin(Theta),    cos(Theta)];
%Trans from X-Y-Z -> t-n-b

g = [Gx;Gy;Gz];
t = [T;Sn;Sb];

ff = Trans *g;
f = t-ff;
f(4) = Omega2;
f(5) = Omega3;

jac = jacobian(f,[u,v,w,T,Sn,Sb,Theta,Phi,Omega2,Omega3]);