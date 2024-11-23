syms Theta Phi Gx Gy Gz u v w;
syms Omega2 Omega3;
syms T Sn Sb;
syms Gx Gy Gz G;
syms rho Cdn Cdb Sdn Sdb;

Trans = [cos(Theta)*cos(Phi),-sin(Phi),cos(Phi)*sin(Theta);cos(Theta)*sin(Phi),cos(Phi),sin(Phi)*sin(Theta);-sin(Theta),0,cos(Theta)];

g = [Gz,Gx,Gy]
t = [T,Sn,Sb]

gg = g*Trans;
hdr1 = 0.5*rho*Cdn*Sdn*v*(v^2+w^2)^(1/2);
hdr2 = 0.5*rho*Cdb*Sdb*w*(v^2+w^2)^(1/2);

f = t - gg;
f(2) = f(2) - hdr1;
f(3) = f(3) - hdr2;
f(4) = Omega2;
f(5) = Omega3;

jac = jacobian(f.',[u,v,w,T,Sn,Sb,Theta,Phi,Omega2,Omega3])
