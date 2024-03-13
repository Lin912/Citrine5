syms Theta Phi V1 V2 V3 u v w;
syms Omega2 Omega3;
syms T Sn Sb;

Trans = [cos(Theta)*cos(Phi),-cos(Theta)*sin(Phi),-sin(Phi);sin(Phi),cos(Phi),0;sin(Theta)*cos(Phi),-sin(Theta)*sin(Phi),cos(Theta)];
V0 = [V1;V2;V3];
u0 = [u;v;w];
f = u0-Trans*V0;
f(4) = Omega2;
f(5) = Omega3;


jac = jacobian(f,[u,v,w,T,Sn,Sb,Theta,Phi,Omega2,Omega3])