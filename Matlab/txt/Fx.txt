syms u0old v0old w0old T0old Sn0old Sb0old Phi0old Theta0old O2mega0old O3mega0old;
syms u0new v0new w0new T0new Sn0new Sb0new Phi0new Theta0new O2mega0new O3mega0new;
syms u1old v1old w1old T1old Sn1old Sb1old Phi1old Theta1old O2mega1old O3mega1old;
syms u1new v1new w1new T1new Sn1new Sb1new Phi1new Theta1new O2mega1new O3mega1new;
syms M ma E A  I g pi rho d0 Cdt Cdn Cdb ;
syms deltaT deltaS;
syms V1 V2 V3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ynew1 = str2sym('[u1new; v1new; w1new; T1new; Sn1new; Sb1new; Phi1new; Theta1new; O2mega1new; O3mega1new;]');
Yold1 = str2sym('[u1old; v1old; w1old; T1old; Sn1old; Sb1old; Phi1old; Theta1old; O2mega1old; O3mega1old;]');
Ynew0 = str2sym('[u0new; v0new; w0new; T0new; Sn0new; Sb0new; Phi0new; Theta0new; O2mega0new; O3mega0new;]');
Yold0 = str2sym('[u0old; v0old; w0old; T0old; Sn0old; Sb0old; Phi0old; Theta0old; O2mega0old; O3mega0old;]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M0old = str2sym('[M 0 0 0 0 0 M*w0old -M*v0old*cos(Theta0old) 0 0;0 M-ma 0 0 0 0 0 M*u0old*cos(Theta0old)+M*w0old*sin(Theta0old)+ma*V1*cos(Phi0old)-ma*V2*sin(Phi0old) 0 0;0 0 M-ma 0 0 0 -M*u0old+ma*V1*cos(Theta0old)*cos(Phi0old)-ma*V2*cos(Theta0old)*sin(Phi0old)-ma*V3*sin(Theta0old) -M*v0old*sin(Theta0old)-ma*V1*sin(Theta0old)*sin(Phi0old)-ma*V2*sin(Theta0old)*cos(Phi0old) 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 1/(E*A) 0 0 0 0 0 0;0 0 0 0 0 0 0 (1+T0old/(E*A))*cos(Theta0old) 0 0;0 0 0 0 0 0 (1+T0old/(E*A)) 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;]');

M0new = str2sym('[M 0 0 0 0 0 M*w0new -M*v0new*cos(Theta0new) 0 0;0 M-ma 0 0 0 0 0 M*u0new*cos(Theta0new)+M*w0new*sin(Theta0new)+ma*V1*cos(Phi0new)-ma*V2*sin(Phi0new) 0 0;0 0 M-ma 0 0 0 -M*u0new+ma*V1*cos(Theta0new)*cos(Phi0new)-ma*V2*cos(Theta0new)*sin(Phi0new)-ma*V3*sin(Theta0new) -M*v0new*sin(Theta0new)-ma*V1*sin(Theta0new)*sin(Phi0new)-ma*V2*sin(Theta0new)*cos(Phi0new) 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 1/(E*A) 0 0 0 0 0 0;0 0 0 0 0 0 0 (1+T0new/(E*A))*cos(Theta0new) 0 0;0 0 0 0 0 0 (1+T0new/(E*A)) 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;]');

M1old = str2sym('[M 0 0 0 0 0 M*w1old -M*v1old*cos(Theta1old) 0 0;0 M-ma 0 0 0 0 0 M*u1old*cos(Theta1old)+M*w1old*sin(Theta1old)+ma*V1*cos(Phi1old)-ma*V2*sin(Phi1old) 0 0;0 0 M-ma 0 0 0 -M*u1old+ma*V1*cos(Theta1old)*cos(Phi1old)-ma*V2*cos(Theta1old)*sin(Phi1old)-ma*V3*sin(Theta1old) -M*v1old*sin(Theta1old)-ma*V1*sin(Theta1old)*sin(Phi1old)-ma*V2*sin(Theta1old)*cos(Phi1old) 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 1/(E*A) 0 0 0 0 0 0;0 0 0 0 0 0 0 (1+T1old/(E*A))*cos(Theta1old) 0 0;0 0 0 0 0 0 (1+T1old/(E*A)) 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;]');

M1new = str2sym('[M 0 0 0 0 0 M*w1new -M*v1new*cos(Theta1new) 0 0;0 M-ma 0 0 0 0 0 M*u1new*cos(Theta1new)+M*w1new*sin(Theta1new)+ma*V1*cos(Phi1new)-ma*V2*sin(Phi1new) 0 0;0 0 M-ma 0 0 0 -M*u1new+ma*V1*cos(Theta1new)*cos(Phi1new)-ma*V2*cos(Theta1new)*sin(Phi1new)-ma*V3*sin(Theta1new) -M*v1new*sin(Theta1new)-ma*V1*sin(Theta1new)*sin(Phi1new)-ma*V2*sin(Theta1new)*cos(Phi1new) 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 1/(E*A) 0 0 0 0 0 0;0 0 0 0 0 0 0 (1+T1new/(E*A))*cos(Theta1new) 0 0;0 0 0 0 0 0 (1+T1new/(E*A)) 0 0 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0;]');

M1 = M1new + M1old;
M2 = M0new + M0old;

Y1 = Ynew1 - Yold1;
Y2 = Ynew0 - Yold0; 

M = (M1 * Y1 + M2 * Y2)*deltaS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N0old = str2sym('[0 0 0 -1 0 0 0 0 0 0;0 0 0 0 -1 0 0 0 0 0;0 0 0 0 0 -1 0 0 0 0;0 0 0 0 0 0 0 0 E*I 0;0 0 0 0 0 0 0 0 0 E*I;-1 0 0 0 0 0 0 0 0 0;0 -1 0 0 0 0 0 0 0 0;0 0 1 0 0 0 0 0 0 0;0 0 0 0 0 0 -1 0 0 0;0 0 0 0 0 0 0 -1*cos(Theta0old) 0 0;]');

N0new = str2sym('[0 0 0 -1 0 0 0 0 0 0;0 0 0 0 -1 0 0 0 0 0;0 0 0 0 0 -1 0 0 0 0;0 0 0 0 0 0 0 0 E*I 0;0 0 0 0 0 0 0 0 0 E*I;-1 0 0 0 0 0 0 0 0 0;0 -1 0 0 0 0 0 0 0 0;0 0 1 0 0 0 0 0 0 0;0 0 0 0 0 0 -1 0 0 0;0 0 0 0 0 0 0 -1*cos(Theta0new) 0 0;]');

N1old = str2sym('[0 0 0 -1 0 0 0 0 0 0;0 0 0 0 -1 0 0 0 0 0;0 0 0 0 0 -1 0 0 0 0;0 0 0 0 0 0 0 0 E*I 0;0 0 0 0 0 0 0 0 0 E*I;-1 0 0 0 0 0 0 0 0 0;0 -1 0 0 0 0 0 0 0 0;0 0 1 0 0 0 0 0 0 0;0 0 0 0 0 0 -1 0 0 0;0 0 0 0 0 0 0 -1*cos(Theta1old) 0 0;]');

N1new = str2sym('[0 0 0 -1 0 0 0 0 0 0;0 0 0 0 -1 0 0 0 0 0;0 0 0 0 0 -1 0 0 0 0;0 0 0 0 0 0 0 0 E*I 0;0 0 0 0 0 0 0 0 0 E*I;-1 0 0 0 0 0 0 0 0 0;0 -1 0 0 0 0 0 0 0 0;0 0 1 0 0 0 0 0 0 0;0 0 0 0 0 0 -1 0 0 0;0 0 0 0 0 0 0 -1*cos(Theta1new) 0 0;]');

N1 =  N1new + N0new;
N2 = N1old + N0old;

Y3 = Ynew1 - Ynew0;
Y4 = Yold1 - Yold0;

N = (N1 * Y3+ N2 * Y4)*deltaT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q0old = str2sym('[-Sb0old*O2mega0old+Sn0old*O3mega0old-M*g*cos(Phi0old)*cos(Theta0old)-0.5*pi*rho*Cdt*(u0old-(V1*cos(Phi0old)*cos(Theta0old)-V2*cos(Theta0old)*sin(Phi0old)-V3*sin(Theta0old)))*abs(u0old-(V1*cos(Phi0old)*cos(Theta0old)-V2*cos(Theta0old)*sin(Phi0old)-V3*sin(Theta0old)))*sqrt(1+(T0old/(E*A)));-T0old*O3mega0old-Sb0old*O3mega0old*tan(Theta0old)-M*g*sin(Phi0old)-0.5*rho*Cdn*(v0old-(V1*sin(Phi0old)+V2*cos(Phi0old)))*sqrt((v0old-(V1*sin(Phi0old)+V2*cos(Phi0old)))^2+(w0old-(V1*cos(Phi0old)*sin(Theta0old)-V2*sin(Phi0old)*sin(Theta0old)+V3*cos(Theta0old)))^2)*sqrt(1+(T0old/(E*A)));-Sn0old*O3mega0old*tan(Theta0old)+T0old*O2mega0old-M*g*cos(Phi0old)*sin(Theta0old)-0.5*rho*Cdb*(w0old-(V1*cos(Phi0old)*sin(Theta0old)-V2*sin(Phi0old)*sin(Theta0old)+V3*cos(Theta0old)))*sqrt((v0old-(V1*sin(Phi0old)+V2*cos(Phi0old)))^2+(w0old-(V1*cos(Phi0old)*sin(Theta0old)-V2*sin(Phi0old)*sin(Theta0old)+V3*cos(Theta0old)))^2)*sqrt(1+(T0old/(E*A)));E*I*O3mega0old*O3mega0old*tan(Theta0old)-Sb0old*(1+T0old/(E*A))*(1+T0old/(E*A))*(1+T0old/(E*A));-E*I*O2mega0old*O3mega0old*tan(Theta0old)+Sn0old(1+T0old/(E*A))*(1+T0old/(E*A))*(1+T0old/(E*A));-O2mega0old*w0old+O3mega0old*v0old;-O3mega0old*(u0old+w0old*tan(Theta0old));-v0old*O3mega0old*tan(Theta0old)-O2mega0old*u0old;O2mega0old;O3mega0old;]');

Q0new = str2sym('[-Sb0new*O2mega0new+Sn0new*O3mega0new-M*g*cos(Phi0new)*cos(Theta0new)-0.5*pi*rho*Cdt*(u0new-(V1*cos(Phi0new)*cos(Theta0new)-V2*cos(Theta0new)*sin(Phi0new)-V3*sin(Theta0new)))*abs(u0new-(V1*cos(Phi0new)*cos(Theta0new)-V2*cos(Theta0new)*sin(Phi0new)-V3*sin(Theta0new)))*sqrt(1+(T0new/(E*A)));-T0new*O3mega0new-Sb0new*O3mega0new*tan(Theta0new)-M*g*sin(Phi0new)-0.5*rho*Cdn*(v0new-(V1*sin(Phi0new)+V2*cos(Phi0new)))*sqrt((v0new-(V1*sin(Phi0new)+V2*cos(Phi0new)))^2+(w0new-(V1*cos(Phi0new)*sin(Theta0new)-V2*sin(Phi0new)*sin(Theta0new)+V3*cos(Theta0new)))^2)*sqrt(1+(T0new/(E*A)));-Sn0new*O3mega0new*tan(Theta0new)+T0new*O2mega0new-M*g*cos(Phi0new)*sin(Theta0new)-0.5*rho*Cdb*(w0new-(V1*cos(Phi0new)*sin(Theta0new)-V2*sin(Phi0new)*sin(Theta0new)+V3*cos(Theta0new)))*sqrt((v0new-(V1*sin(Phi0new)+V2*cos(Phi0new)))^2+(w0new-(V1*cos(Phi0new)*sin(Theta0new)-V2*sin(Phi0new)*sin(Theta0new)+V3*cos(Theta0new)))^2)*sqrt(1+(T0new/(E*A)));E*I*O3mega0new*O3mega0new*tan(Theta0new)-Sb0new*(1+T0new/(E*A))*(1+T0new/(E*A))*(1+T0new/(E*A));-E*I*O2mega0new*O3mega0new*tan(Theta0new)+Sn0new(1+T0new/(E*A))*(1+T0new/(E*A))*(1+T0new/(E*A));-O2mega0new*w0new+O3mega0new*v0new;-O3mega0new*(u0new+w0new*tan(Theta0new));-v0new*O3mega0new*tan(Theta0new)-O2mega0new*u0new;O2mega0new;O3mega0new;]');

Q1new = str2sym('[-Sb1new*O2mega1new+Sn1new*O3mega1new-M*g*cos(Phi1new)*cos(Theta1new)-0.5*pi*rho*Cdt*(u1new-(V1*cos(Phi1new)*cos(Theta1new)-V2*cos(Theta1new)*sin(Phi1new)-V3*sin(Theta1new)))*abs(u1new-(V1*cos(Phi1new)*cos(Theta1new)-V2*cos(Theta1new)*sin(Phi1new)-V3*sin(Theta1new)))*sqrt(1+(T1new/(E*A)));-T1new*O3mega1new-Sb1new*O3mega1new*tan(Theta1new)-M*g*sin(Phi1new)-0.5*rho*Cdn*(v1new-(V1*sin(Phi1new)+V2*cos(Phi1new)))*sqrt((v1new-(V1*sin(Phi1new)+V2*cos(Phi1new)))^2+(w1new-(V1*cos(Phi1new)*sin(Theta1new)-V2*sin(Phi1new)*sin(Theta1new)+V3*cos(Theta1new)))^2)*sqrt(1+(T1new/(E*A)));-Sn1new*O3mega1new*tan(Theta1new)+T1new*O2mega1new-M*g*cos(Phi1new)*sin(Theta1new)-0.5*rho*Cdb*(w1new-(V1*cos(Phi1new)*sin(Theta1new)-V2*sin(Phi1new)*sin(Theta1new)+V3*cos(Theta1new)))*sqrt((v1new-(V1*sin(Phi1new)+V2*cos(Phi1new)))^2+(w1new-(V1*cos(Phi1new)*sin(Theta1new)-V2*sin(Phi1new)*sin(Theta1new)+V3*cos(Theta1new)))^2)*sqrt(1+(T1new/(E*A)));E*I*O3mega1new*O3mega1new*tan(Theta1new)-Sb1new*(1+T1new/(E*A))*(1+T1new/(E*A))*(1+T1new/(E*A));-E*I*O2mega1new*O3mega1new*tan(Theta1new)+Sn1new(1+T1new/(E*A))*(1+T1new/(E*A))*(1+T1new/(E*A));-O2mega1new*w1new+O3mega1new*v1new;-O3mega1new*(u1new+w1new*tan(Theta1new));-v1new*O3mega1new*tan(Theta1new)-O2mega1new*u1new;O2mega1new;O3mega1new;]');

Q1old = str2sym('[-Sb1old*O2mega1old+Sn1old*O3mega1old-M*g*cos(Phi1old)*cos(Theta1old)-0.5*pi*rho*Cdt*(u1old-(V1*cos(Phi1old)*cos(Theta1old)-V2*cos(Theta1old)*sin(Phi1old)-V3*sin(Theta1old)))*abs(u1old-(V1*cos(Phi1old)*cos(Theta1old)-V2*cos(Theta1old)*sin(Phi1old)-V3*sin(Theta1old)))*sqrt(1+(T1old/(E*A)));-T1old*O3mega1old-Sb1old*O3mega1old*tan(Theta1old)-M*g*sin(Phi1old)-0.5*rho*Cdn*(v1old-(V1*sin(Phi1old)+V2*cos(Phi1old)))*sqrt((v1old-(V1*sin(Phi1old)+V2*cos(Phi1old)))^2+(w1old-(V1*cos(Phi1old)*sin(Theta1old)-V2*sin(Phi1old)*sin(Theta1old)+V3*cos(Theta1old)))^2)*sqrt(1+(T1old/(E*A)));-Sn1old*O3mega1old*tan(Theta1old)+T1old*O2mega1old-M*g*cos(Phi1old)*sin(Theta1old)-0.5*rho*Cdb*(w1old-(V1*cos(Phi1old)*sin(Theta1old)-V2*sin(Phi1old)*sin(Theta1old)+V3*cos(Theta1old)))*sqrt((v1old-(V1*sin(Phi1old)+V2*cos(Phi1old)))^2+(w1old-(V1*cos(Phi1old)*sin(Theta1old)-V2*sin(Phi1old)*sin(Theta1old)+V3*cos(Theta1old)))^2)*sqrt(1+(T1old/(E*A)));E*I*O3mega1old*O3mega1old*tan(Theta1old)-Sb1old*(1+T1old/(E*A))*(1+T1old/(E*A))*(1+T1old/(E*A));-E*I*O2mega1old*O3mega1old*tan(Theta1old)+Sn1old(1+T1old/(E*A))*(1+T1old/(E*A))*(1+T1old/(E*A));-O2mega1old*w1old+O3mega1old*v1old;-O3mega1old*(u1old+w1old*tan(Theta1old));-v1old*O3mega1old*tan(Theta1old)-O2mega1old*u1old;O2mega1old;O3mega1old;]');

Q = (Q0new + Q0old + Q1old + Q1new)*deltaT  *deltaS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = M + N + Q;

Fj0 = jacobian(F, [u0new, v0new, w0new, T0new, Sn0new, Sb0new, Theta0new,Phi0new, O2mega0new, O3mega0new]);
Fj1 = jacobian(F, [u1new, v1new, w1new, T1new, Sn1new, Sb1new, Theta1new, Phi1new, O2mega1new, O3mega1new]);