syms Theta Phi V1 V2 V3 u v w;
syms Omega2 Omega3;
syms T Sn Sb;
syms Gx Gy Gz G;
syms rho Cdn Cdt Sd;

Trans = [cos(Theta)*cos(Phi),-cos(Theta)*sin(Phi),-sin(Phi);sin(Phi),cos(Phi),0;sin(Theta)*cos(Phi),-sin(Theta)*sin(Phi),cos(Theta)];

g = [Gx;Gy;Gz]
t = [T;Sn;Sb]

gg = Trans*g;
hdr1 = 0.5*rho*Cdn*Sd*v*(v^2+w^2)^(1/2);
hdr2 = 0.5*rho*Cdn*Sd*w*(v^2+w^2)^(1/2);

f = t + gg;
f(2) = f(2) + hdr1;
f(3) = f(3) + hdr2;
f(4) = Omega2;
f(5) = Omega3;

jacobian(f,[u,v,w,T,Sn,Sb,Theta,Phi,Omega2,Omega3])
