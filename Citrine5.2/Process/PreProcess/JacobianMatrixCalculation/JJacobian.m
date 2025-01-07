clear;
clc;
syms u0old v0old w0old T0old Sn0old Sb0old Phi0old Theta0old O2mega0old O3mega0old real;
syms u0new v0new w0new T0new Sn0new Sb0new Phi0new Theta0new O2mega0new O3mega0new real;
syms u1old v1old w1old T1old Sn1old Sb1old Phi1old Theta1old O2mega1old O3mega1old real;
syms u1new v1new w1new T1new Sn1new Sb1new Phi1new Theta1new O2mega1new O3mega1new real;
syms M ma E A I pi rho d0 Cdt Cdn Cdb real;
syms deltaT deltaS real;
syms Vx Vy Vz real;
syms Gx Gy Gz real;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ynew1 = str2sym('[u1new; v1new; w1new; T1new; Sn1new; Sb1new; Theta1new; Phi1new; O2mega1new; O3mega1new;]');
Yold1 = str2sym('[u1old; v1old; w1old; T1old; Sn1old; Sb1old; Theta1old; Phi1old; O2mega1old; O3mega1old;]');
Ynew0 = str2sym('[u0new; v0new; w0new; T0new; Sn0new; Sb0new; Theta0new; Phi0new; O2mega0new; O3mega0new;]');
Yold0 = str2sym('[u0old; v0old; w0old; T0old; Sn0old; Sb0old; Theta0old; Phi0old; O2mega0old; O3mega0old;]');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M0old = str2sym(['[M, 0, 0, 0, 0, 0, M*w0old, -M*cos(Theta0old)*v0old, 0, 0; ' ...
                 '0, M+ma, 0, 0, 0, 0, 0, M*u0old*cos(Theta0old) + M*w0old*sin(Theta0old) + ma*Vy*sin(Phi0old) + ma*Vx*cos(Phi0old), 0, 0; ' ...
                 '0, 0, M+ma, 0, 0, 0, -M*u0old - ma*Vx*cos(Phi0old)*cos(Theta0old) - ma*Vy*sin(Phi0old)*cos(Theta0old) + ma*Vz*sin(Theta0old), -M*v0old*sin(Theta0old) + ma*Vx*sin(Phi0old)*sin(Theta0old) - ma*Vy*cos(Phi0old)*sin(Theta0old), 0, 0; ' ...
                 '0, 0, 0, (1/(E*A)), 0, 0, 0, 0, 0, 0; ' ...
                 '0, 0, 0, 0, 0, 0, 0, (1+(T0old/(E*A)))*cos(Theta0old), 0, 0; ' ...
                 '0, 0, 0, 0, 0, 0, (1+(T0old/(E*A))), 0, 0, 0; ' ...
                 '0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ' ...
                 '0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ' ...
                 '0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ' ...
                 '0, 0, 0, 0, 0, 0, 0, 0, 0, 0]']);

M1old = str2sym(['[M, 0, 0, 0, 0, 0, M*w1old, -M*cos(Theta1old)*v1old, 0, 0; ' ...
                 '0, M+ma, 0, 0, 0, 0, 0, M*u1old*cos(Theta1old) + M*w1old*sin(Theta1old) + ma*Vy*sin(Phi1old) + ma*Vx*cos(Phi1old), 0, 0; ' ...
                 '0, 0, M+ma, 0, 0, 0, -M*u1old - ma*Vx*cos(Phi1old)*cos(Theta1old) - ma*Vy*sin(Phi1old)*cos(Theta1old) + ma*Vz*sin(Theta1old), -M*v1old*sin(Theta1old) + ma*Vx*sin(Phi1old)*sin(Theta1old) - ma*Vy*cos(Phi1old)*sin(Theta1old), 0, 0; ' ...
                 '0, 0, 0, (1/(E*A)), 0, 0, 0, 0, 0, 0; ' ...
                 '0, 0, 0, 0, 0, 0, 0, (1+(T1old/(E*A)))*cos(Theta1old), 0, 0; ' ...
                 '0, 0, 0, 0, 0, 0, (1+(T1old/(E*A))), 0, 0, 0; ' ...
                 '0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ' ...
                 '0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ' ...
                 '0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ' ...
                 '0, 0, 0, 0, 0, 0, 0, 0, 0, 0]']);

M0new = str2sym(['[M, 0, 0, 0, 0, 0, M*w0new, -M*cos(Theta0new)*v0new, 0, 0; ' ...
                 '0, M+ma, 0, 0, 0, 0, 0, M*u0new*cos(Theta0new) + M*w0new*sin(Theta0new) + ma*Vy*sin(Phi0new) + ma*Vx*cos(Phi0new), 0, 0; ' ...
                 '0, 0, M+ma, 0, 0, 0, -M*u0new - ma*Vx*cos(Phi0new)*cos(Theta0new) - ma*Vy*sin(Phi0new)*cos(Theta0new) + ma*Vz*sin(Theta0new), -M*v0new*sin(Theta0new) + ma*Vx*sin(Phi0new)*sin(Theta0new) - ma*Vy*cos(Phi0new)*sin(Theta0new), 0, 0; ' ...
                 '0, 0, 0, (1/(E*A)), 0, 0, 0, 0, 0, 0; ' ...
                 '0, 0, 0, 0, 0, 0, 0, (1+(T0new/(E*A)))*cos(Theta0new), 0, 0; ' ...
                 '0, 0, 0, 0, 0, 0, (1+(T0new/(E*A))), 0, 0, 0; ' ...
                 '0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ' ...
                 '0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ' ...
                 '0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ' ...
                 '0, 0, 0, 0, 0, 0, 0, 0, 0, 0]']);

M1new = str2sym(['[M, 0, 0, 0, 0, 0, M*w1new, -M*cos(Theta1new)*v1new, 0, 0; ' ...
                 '0, M+ma, 0, 0, 0, 0, 0, M*u1new*cos(Theta1new) + M*w1new*sin(Theta1new) + ma*Vy*sin(Phi1new) + ma*Vx*cos(Phi1new), 0, 0; ' ...
                 '0, 0, M+ma, 0, 0, 0, -M*u1new - ma*Vx*cos(Phi1new)*cos(Theta1new) - ma*Vy*sin(Phi1new)*cos(Theta1new) + ma*Vz*sin(Theta1new), -M*v1new*sin(Theta1new) + ma*Vx*sin(Phi1new)*sin(Theta1new) - ma*Vy*cos(Phi1new)*sin(Theta1new), 0, 0; ' ...
                 '0, 0, 0, (1/(E*A)), 0, 0, 0, 0, 0, 0; ' ...
                 '0, 0, 0, 0, 0, 0, 0, (1+(T1new/(E*A)))*cos(Theta1new), 0, 0; ' ...
                 '0, 0, 0, 0, 0, 0, (1+(T1new/(E*A))), 0, 0, 0; ' ...
                 '0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ' ...
                 '0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ' ...
                 '0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ' ...
                 '0, 0, 0, 0, 0, 0, 0, 0, 0, 0]']);

M1 = M1new + M1old;
M2 = M0new + M0old;
Y1 = Ynew1 - Yold1;
Y2 = Ynew0 - Yold0;

MM = (M1 * Y1 + M2 * Y2)*deltaS;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N0old = str2sym('[0,0,0,-1,0,0,0,0,0,0;0,0,0,0,-1,0,0,0,0,0;0,0,0,0,0,-1,0,0,0,0;-1,0,0,0,0,0,0,0,0,0;0,-1,0,0,0,0,0,0,0,0;0,0,1,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,E*I,0;0,0,0,0,0,0,0,0,0,E*I;0,0,0,0,0,0,1,0,0,0;0,0,0,0,0,0,0,cos(Theta0old),0,0;]');
N1old = str2sym('[0,0,0,-1,0,0,0,0,0,0;0,0,0,0,-1,0,0,0,0,0;0,0,0,0,0,-1,0,0,0,0;-1,0,0,0,0,0,0,0,0,0;0,-1,0,0,0,0,0,0,0,0;0,0,1,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,E*I,0;0,0,0,0,0,0,0,0,0,E*I;0,0,0,0,0,0,1,0,0,0;0,0,0,0,0,0,0,cos(Theta1old),0,0;]');
N1new = str2sym('[0,0,0,-1,0,0,0,0,0,0;0,0,0,0,-1,0,0,0,0,0;0,0,0,0,0,-1,0,0,0,0;-1,0,0,0,0,0,0,0,0,0;0,-1,0,0,0,0,0,0,0,0;0,0,1,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,E*I,0;0,0,0,0,0,0,0,0,0,E*I;0,0,0,0,0,0,1,0,0,0;0,0,0,0,0,0,0,cos(Theta1new),0,0;]');
N0new = str2sym('[0,0,0,-1,0,0,0,0,0,0;0,0,0,0,-1,0,0,0,0,0;0,0,0,0,0,-1,0,0,0,0;-1,0,0,0,0,0,0,0,0,0;0,-1,0,0,0,0,0,0,0,0;0,0,1,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,E*I,0;0,0,0,0,0,0,0,0,0,E*I;0,0,0,0,0,0,1,0,0,0;0,0,0,0,0,0,0,cos(Theta0new),0,0;]');

N1 = N1new + N0new;
N2 = N1old + N0old;
Y3 = Ynew1 - Ynew0;
Y4 = Yold1 - Yold0;

NN = (N1 * Y3+ N2 * Y4)*deltaT;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q0old = str2sym('[-O2mega0old*Sb0old+O3mega0old*Sn0old+(-Gx*cos(Phi0old)*cos(Theta0old) - Gy*cos(Theta0old)*sin(Phi0old) + Gz*sin(Theta0old))+0.5*pi*rho*d0*Cdt*(u0old - Vx*cos(Phi0old)*cos(Theta0old) - Vy*cos(Theta0old)*sin(Phi0old) + Vz* sin(Theta0old))*abs(u0old - Vx*cos(Phi0old)*cos(Theta0old) - Vy*cos(Theta0old)*sin(Phi0old) + Vz* sin(Theta0old))*sqrt(1+(T0old/(E*A)));-O3mega0old*T0old-tan(Theta0old)*O3mega0old*Sb0old+(-Gy*cos(Phi0old) + Gx*sin(Phi0old))+0.5*rho*d0*Cdn*(v0old - Vy*cos(Phi0old) + Vx*sin(Phi0old))*sqrt((v0old - Vy*cos(Phi0old) + Vx*sin(Phi0old))^2+(w0old - Vx*cos(Phi0old)*sin(Theta0old) -Vy*sin(Phi0old)*sin(Theta0old) - Vz*cos(Theta0old))^2)*sqrt(1+(T0old/(E*A)));O3mega0old*tan(Theta0old)*Sn0old+O2mega0old*T0old+(-Gx*cos(Phi0old)*sin(Theta0old) - Gy*sin(Phi0old)*sin(Theta0old) - Gz*cos(Theta0old))+0.5*rho*d0*Cdb*(w0old - Vx*cos(Phi0old)*sin(Theta0old) -Vy*sin(Phi0old)*sin(Theta0old) - Vz*cos(Theta0old))*sqrt((v0old - Vy*cos(Phi0old) + Vx*sin(Phi0old))^2+(w0old - Vx*cos(Phi0old)*sin(Theta0old) -Vy*sin(Phi0old)*sin(Theta0old) - Vz*cos(Theta0old))^2)*sqrt(1+(T0old/(E*A)));-O2mega0old*w0old+O3mega0old*v0old;-O3mega0old*u0old-tan(Theta0old)*O3mega0old*w0old;-tan(Theta0old)*O3mega0old*v0old-O2mega0old*u0old;E*I*O3mega0old*O3mega0old*tan(Theta0old)-Sb0old*(1+(T0old/(E*A)))^3;-E*I*O2mega0old*O3mega0old*tan(Theta0old)+Sn0old*(1+(T0old/(E*A)))^3;-O2mega0old;-O3mega0old]');
Q0new = str2sym('[-O2mega0new*Sb0new+O3mega0new*Sn0new+(-Gx*cos(Phi0new)*cos(Theta0new) - Gy*cos(Theta0new)*sin(Phi0new) + Gz*sin(Theta0new))+0.5*pi*rho*d0*Cdt*(u0new - Vx*cos(Phi0new)*cos(Theta0new) - Vy*cos(Theta0new)*sin(Phi0new) + Vz* sin(Theta0new))*abs(u0new - Vx*cos(Phi0new)*cos(Theta0new) - Vy*cos(Theta0new)*sin(Phi0new) + Vz* sin(Theta0new))*sqrt(1+(T0new/(E*A)));-O3mega0new*T0new-tan(Theta0new)*O3mega0new*Sb0new+(-Gy*cos(Phi0new) + Gx*sin(Phi0new))+0.5*rho*d0*Cdn*(v0new - Vy*cos(Phi0new) + Vx*sin(Phi0new))*sqrt((v0new - Vy*cos(Phi0new) + Vx*sin(Phi0new))^2+(w0new - Vx*cos(Phi0new)*sin(Theta0new) -Vy*sin(Phi0new)*sin(Theta0new) - Vz*cos(Theta0new))^2)*sqrt(1+(T0new/(E*A)));O3mega0new*tan(Theta0new)*Sn0new+O2mega0new*T0new+(-Gx*cos(Phi0new)*sin(Theta0new) - Gy*sin(Phi0new)*sin(Theta0new) - Gz*cos(Theta0new))+0.5*rho*d0*Cdb*(w0new - Vx*cos(Phi0new)*sin(Theta0new) -Vy*sin(Phi0new)*sin(Theta0new) - Vz*cos(Theta0new))*sqrt((v0new - Vy*cos(Phi0new) + Vx*sin(Phi0new))^2+(w0new - Vx*cos(Phi0new)*sin(Theta0new) -Vy*sin(Phi0new)*sin(Theta0new) - Vz*cos(Theta0new))^2)*sqrt(1+(T0new/(E*A)));-O2mega0new*w0new+O3mega0new*v0new;-O3mega0new*u0new-tan(Theta0new)*O3mega0new*w0new;-tan(Theta0new)*O3mega0new*v0new-O2mega0new*u0new;E*I*O3mega0new*O3mega0new*tan(Theta0new)-Sb0new*(1+(T0new/(E*A)))^3;-E*I*O2mega0new*O3mega0new*tan(Theta0new)+Sn0new*(1+(T0new/(E*A)))^3;-O2mega0new;-O3mega0new]');
Q1old = str2sym('[-O2mega1old*Sb1old+O3mega1old*Sn1old+(-Gx*cos(Phi1old)*cos(Theta1old) - Gy*cos(Theta1old)*sin(Phi1old) + Gz*sin(Theta1old))+0.5*pi*rho*d0*Cdt*(u1old - Vx*cos(Phi1old)*cos(Theta1old) - Vy*cos(Theta1old)*sin(Phi1old) + Vz* sin(Theta1old))*abs(u1old - Vx*cos(Phi1old)*cos(Theta1old) - Vy*cos(Theta1old)*sin(Phi1old) + Vz* sin(Theta1old))*sqrt(1+(T1old/(E*A)));-O3mega1old*T1old-tan(Theta1old)*O3mega1old*Sb1old+(-Gy*cos(Phi1old) + Gx*sin(Phi1old))+0.5*rho*d0*Cdn*(v1old - Vy*cos(Phi1old) + Vx*sin(Phi1old))*sqrt((v1old - Vy*cos(Phi1old) + Vx*sin(Phi1old))^2+(w1old - Vx*cos(Phi1old)*sin(Theta1old) -Vy*sin(Phi1old)*sin(Theta1old) - Vz*cos(Theta1old))^2)*sqrt(1+(T1old/(E*A)));O3mega1old*tan(Theta1old)*Sn1old+O2mega1old*T1old+(-Gx*cos(Phi1old)*sin(Theta1old) - Gy*sin(Phi1old)*sin(Theta1old) - Gz*cos(Theta1old))+0.5*rho*d0*Cdb*(w1old - Vx*cos(Phi1old)*sin(Theta1old) -Vy*sin(Phi1old)*sin(Theta1old) - Vz*cos(Theta1old))*sqrt((v1old - Vy*cos(Phi1old) + Vx*sin(Phi1old))^2+(w1old - Vx*cos(Phi1old)*sin(Theta1old) -Vy*sin(Phi1old)*sin(Theta1old) - Vz*cos(Theta1old))^2)*sqrt(1+(T1old/(E*A)));-O2mega1old*w1old+O3mega1old*v1old;-O3mega1old*u1old-tan(Theta1old)*O3mega1old*w1old;-tan(Theta1old)*O3mega1old*v1old-O2mega1old*u1old;E*I*O3mega1old*O3mega1old*tan(Theta1old)-Sb1old*(1+(T1old/(E*A)))^3;-E*I*O2mega1old*O3mega1old*tan(Theta1old)+Sn1old*(1+(T1old/(E*A)))^3;-O2mega1old;-O3mega1old]');
Q1new = str2sym('[-O2mega1new*Sb1new+O3mega1new*Sn1new+(-Gx*cos(Phi1new)*cos(Theta1new) - Gy*cos(Theta1new)*sin(Phi1new) + Gz*sin(Theta1new))+0.5*pi*rho*d0*Cdt*(u1new - Vx*cos(Phi1new)*cos(Theta1new) - Vy*cos(Theta1new)*sin(Phi1new) + Vz* sin(Theta1new))*abs(u1new - Vx*cos(Phi1new)*cos(Theta1new) - Vy*cos(Theta1new)*sin(Phi1new) + Vz* sin(Theta1new))*sqrt(1+(T1new/(E*A)));-O3mega1new*T1new-tan(Theta1new)*O3mega1new*Sb1new+(-Gy*cos(Phi1new) + Gx*sin(Phi1new))+0.5*rho*d0*Cdn*(v1new - Vy*cos(Phi1new) + Vx*sin(Phi1new))*sqrt((v1new - Vy*cos(Phi1new) + Vx*sin(Phi1new))^2+(w1new - Vx*cos(Phi1new)*sin(Theta1new) -Vy*sin(Phi1new)*sin(Theta1new) - Vz*cos(Theta1new))^2)*sqrt(1+(T1new/(E*A)));O3mega1new*tan(Theta1new)*Sn1new+O2mega1new*T1new+(-Gx*cos(Phi1new)*sin(Theta1new) - Gy*sin(Phi1new)*sin(Theta1new) - Gz*cos(Theta1new))+0.5*rho*d0*Cdb*(w1new - Vx*cos(Phi1new)*sin(Theta1new) -Vy*sin(Phi1new)*sin(Theta1new) - Vz*cos(Theta1new))*sqrt((v1new - Vy*cos(Phi1new) + Vx*sin(Phi1new))^2+(w1new - Vx*cos(Phi1new)*sin(Theta1new) -Vy*sin(Phi1new)*sin(Theta1new) - Vz*cos(Theta1new))^2)*sqrt(1+(T1new/(E*A)));-O2mega1new*w1new+O3mega1new*v1new;-O3mega1new*u1new-tan(Theta1new)*O3mega1new*w1new;-tan(Theta1new)*O3mega1new*v1new-O2mega1new*u1new;E*I*O3mega1new*O3mega1new*tan(Theta1new)-Sb1new*(1+(T1new/(E*A)))^3;-E*I*O2mega1new*O3mega1new*tan(Theta1new)+Sn1new*(1+(T1new/(E*A)))^3;-O2mega1new;-O3mega1new]');

QQ = (Q0new + Q0old + Q1old + Q1new)*deltaT*deltaS;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = MM + NN + QQ;

Fj0 = jacobian(F, [u0new, v0new, w0new, T0new, Sn0new, Sb0new, Theta0new, Phi0new, O2mega0new, O3mega0new]);
Fj1 = jacobian(F, [u1new, v1new, w1new, T1new, Sn1new, Sb1new, Theta1new, Phi1new, O2mega1new, O3mega1new]);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outputFolder = 'FJ0 output_txt_files';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

[numRows, numCols] = size(Fj0);
for row = 1:numRows
    for col = 1:numCols
        data = Fj0(row, col);
        filename = fullfile(outputFolder, sprintf('Fj0_%d_%d.txt', row, col));
        fid = fopen(filename, 'w');
        fprintf(fid, '%s', char(data));
        fclose(fid);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outputFolder = 'FJ1 output_txt_files';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

[numRows, numCols] = size(Fj1);
for row = 1:numRows
    for col = 1:numCols
        data = Fj1(row, col);
        filename = fullfile(outputFolder, sprintf('Fj1_%d_%d.txt', row, col));
        fid = fopen(filename, 'w');
        fprintf(fid, '%s', char(data));
        fclose(fid);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('All files have been saved to the output_txt_files folder');









