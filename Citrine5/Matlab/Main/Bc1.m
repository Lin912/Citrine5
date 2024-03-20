syms Theta Phi Vx Vy Vz u v w;
syms Omega2 Omega3;
syms T Sn Sb;

Trans = [cos(Theta)*cos(Phi),-sin(Phi),cos(Phi)*sin(Theta);cos(Theta)*sin(Phi),cos(Phi),sin(Phi)*sin(Theta);-sin(Theta),0,cos(Theta)];
V0 = [Vz;Vx;Vy];
u0 = [u;v;w];
f = u0-Trans*V0;
f(4) = Omega2;
f(5) = Omega3;


jac = jacobian(f,[u,v,w,T,Sn,Sb,Theta,Phi,Omega2,Omega3])