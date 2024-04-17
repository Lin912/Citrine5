clear;
clc;

data = importdata('newoutput2.csv');
F = size(data);
Fx = [5000,50];
Fy = [5000,50];
Fz = [5000,50];

%提取全部Force
if mod(size(data,2),10) ==0
    flag = 1;
else
    flag = 0;
end
if flag == 1
    for j = 1:size(data,2)/10
        for i = 1:size(data,1)
            F(i,3*(j-1)+1:3*(j-1)+3)= data(i,10*(j-1)+4:10*(j-1)+6);
        end
    end
else
end


t = 0.001:0.001:5;

%提取ForceX(Tension)
 for j = 1:size(F,2)/3
        for i = 1:5000
            Fx(i,1*(j-1)+1)= F(i,3*(j-1)+1);
        end
 end

%提取ForceY
 for j = 1:size(F,2)/3
        for i = 1:5000
            Fy(i,1*(j-1)+1)= F(i,3*(j-1)+2);
        end
 end
 
%提取ForceZ
 for j = 1:size(F,2)/3
        for i = 1:5000
            Fz(i,1*(j-1)+1)= F(i,3*(j-1)+3);
        end
 end

%提取(i)号节点Tension数据
%Ti = Fx(:,i);
%plot(t,Ti);

subplot(3,3,1);
T11 = Fy(:,1);
plot(t,T11,"-",'LineWidth',1.2);
xlim([0 10]);
ylim([-10 10]);
set(gca,'XTick',0:1:10);
set(gca,'YTick',-10:2:10);
title("Tension at Point1",'fontsize',16,'fontname','times new roman');
xlabel("Time(s)",'fontname','times new roman');
ylabel("Force(N)",'fontname','times new roman');
grid on;
hold on;

subplot(3,3,4);
T12 = Fx(:,1);
plot(t,T12,"-",'LineWidth',1.2);
xlim([0 10]);
ylim([-10 10]);
set(gca,'XTick',0:1:10);
set(gca,'YTick',-10:4:10);
title("Fy at Point1",'fontsize',16,'fontname','times new roman');
xlabel("Time(s)",'fontname','times new roman');
ylabel("Force(N)",'fontname','times new roman');
grid on;

subplot(3,3,7);
T13 = Fz(:,1);
plot(t,T13,"-",'LineWidth',1.2);
xlim([0 10]);
ylim([-10 10]);
set(gca,'XTick',0:1:10);
set(gca,'YTick',-10:4:10);
title("Fz at Point1",'fontsize',16,'fontname','times new roman');
xlabel("Time(s)",'fontname','times new roman');
ylabel("Force(N)",'fontname','times new roman');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,3,2);
T21 = Fy(:,25);
plot(t,T21,"-",'LineWidth',1.2);
xlim([0 10]);
ylim([-10 10]);
set(gca,'XTick',0:1:10);
set(gca,'YTick',-10:2:10);
title("Tension at Point25",'fontsize',16,'fontname','times new roman');
xlabel("Time(s)",'fontname','times new roman');
ylabel("Force(N)",'fontname','times new roman');
grid on;
hold on;

subplot(3,3,5);
T22 = Fx(:,25);
plot(t,T22,"-",'LineWidth',1.2);
xlim([0 10]);
ylim([-10 10]);
set(gca,'XTick',0:1:10);
set(gca,'YTick',-10:4:10);
title("Fy at Point25",'fontsize',16,'fontname','times new roman');
xlabel("Time(s)",'fontname','times new roman');
ylabel("Force(N)",'fontname','times new roman');
grid on;

subplot(3,3,8);
T23 = Fz(:,25);
plot(t,T23,"-",'LineWidth',1.2);
xlim([0 10]);
ylim([-10 10]);
set(gca,'XTick',0:1:10);
set(gca,'YTick',-10:4:10);
title("Fz at Point25",'fontsize',16,'fontname','times new roman');
xlabel("Time(s)",'fontname','times new roman');
ylabel("Force(N)",'fontname','times new roman');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,3,3);
T31 = Fy(:,50);
plot(t,T31,"-",'LineWidth',1.2);
xlim([0 10]);
ylim([-10 10]);
set(gca,'XTick',0:1:10);
set(gca,'YTick',-10:2:10);
title("Tension at Point50",'fontsize',16,'fontname','times new roman');
xlabel("Time(s)",'fontname','times new roman');
ylabel("Force(N)",'fontname','times new roman');
grid on;
hold on;

subplot(3,3,6);
T32 = Fx(:,50);
plot(t,T32,"-",'LineWidth',1.2);
xlim([0 10]);
ylim([-10 10]);
set(gca,'XTick',0:1:10);
set(gca,'YTick',-10:4:10);
title("Fy at Point50",'fontsize',16,'fontname','times new roman');
xlabel("Time(s)",'fontname','times new roman');
ylabel("Force(N)",'fontname','times new roman');
grid on;

subplot(3,3,9);
T33 = Fz(:,50);
plot(t,T33,"-",'LineWidth',1.2);
xlim([0 10]);
ylim([-10 10]);
set(gca,'XTick',0:1:10);
set(gca,'YTick',-10:4:10);
title("Fz at Point50",'fontsize',16,'fontname','times new roman');
xlabel("Time(s)",'fontname','times new roman');
ylabel("Force(N)",'fontname','times new roman');
grid on;

