clear;
clc;
syms u0old v0old w0old T0old Sn0old Sb0old Phi0old Theta0old O2mega0old O3mega0old real;
syms u0new v0new w0new T0new Sn0new Sb0new Phi0new Theta0new O2mega0new O3mega0new real;
syms u1old v1old w1old T1old Sn1old Sb1old Phi1old Theta1old O2mega1old O3mega1old real;
syms u1new v1new w1new T1new Sn1new Sb1new Phi1new Theta1new O2mega1new O3mega1new real;
syms M ma E A I w0 pi rho d0 Cdt Cdn Cdb real;
syms deltaT deltaS real;
syms Vx Vy Vz real;

N0old = str2sym('[0,0,0,-1,0,0,0,0,0,0;0,0,0,0,-1,0,0,0,0,0;0,0,0,0,0,-1,0,0,0,0;-1,0,0,0,0,0,0,0,0,0;0,-1,0,0,0,0,0,0,0,0;0,0,1,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,E*I,0;0,0,0,0,0,0,0,0,0,E*I;0,0,0,0,0,0,1,0,0,0;0,0,0,0,0,0,0,cos(Theta0old),0,0;]');

N1old = str2sym('[0,0,0,-1,0,0,0,0,0,0;0,0,0,0,-1,0,0,0,0,0;0,0,0,0,0,-1,0,0,0,0;-1,0,0,0,0,0,0,0,0,0;0,-1,0,0,0,0,0,0,0,0;0,0,1,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,E*I,0;0,0,0,0,0,0,0,0,0,E*I;0,0,0,0,0,0,1,0,0,0;0,0,0,0,0,0,0,cos(Theta1old),0,0;]');

N1new = str2sym('[0,0,0,-1,0,0,0,0,0,0;0,0,0,0,-1,0,0,0,0,0;0,0,0,0,0,-1,0,0,0,0;-1,0,0,0,0,0,0,0,0,0;0,-1,0,0,0,0,0,0,0,0;0,0,1,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,E*I,0;0,0,0,0,0,0,0,0,0,E*I;0,0,0,0,0,0,1,0,0,0;0,0,0,0,0,0,0,cos(Theta1new),0,0;]');

N0new = str2sym('[0,0,0,-1,0,0,0,0,0,0;0,0,0,0,-1,0,0,0,0,0;0,0,0,0,0,-1,0,0,0,0;-1,0,0,0,0,0,0,0,0,0;0,-1,0,0,0,0,0,0,0,0;0,0,1,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,E*I,0;0,0,0,0,0,0,0,0,0,E*I;0,0,0,0,0,0,1,0,0,0;0,0,0,0,0,0,0,cos(Theta0new),0,0;]');