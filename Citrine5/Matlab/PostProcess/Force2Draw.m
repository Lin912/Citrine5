clear;
clc;

data = importdata('newoutput2.csv');
F = size(data);
Fx = [2000,50];
Fy = [2000,50];
Fz = [2000,50];

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


t = 0.002:0.002:4;

%提取ForceX(Tension)
 for j = 1:size(F,2)/3
        for i = 1:2000
            Fx(i,1*(j-1)+1)= F(i,3*(j-1)+1);
        end
 end

%提取ForceY
 for j = 1:size(F,2)/3
        for i = 1:2000
            Fy(i,1*(j-1)+1)= F(i,3*(j-1)+2);
        end
 end
 
%提取ForceZ
 for j = 1:size(F,2)/3
        for i = 1:2000
            Fz(i,1*(j-1)+1)= F(i,3*(j-1)+3);
        end
 end

%提取(i)号节点Tension数据
%Ti = Fx(:,i);
%plot(t,Ti);

subplot(2,1,1);
T1 = Fx(:,1);
plot(t,T1,"-",'LineWidth',1.2);
xlim([0 4]);
ylim([-6 6]);
set(gca,'XTick',0:0.2:4);
set(gca,'YTick',-6:2:6);
title("Force at Point1",'fontsize',16,'fontname','times new roman');
xlabel("Time(s)",'fontname','times new roman');
ylabel("Force(N)",'fontname','times new roman');
grid on;
hold on;

subplot(2,1,2);
T2 = Fy(:,1);
plot(t,T2,"-",'LineWidth',1.2);
xlim([0 4]);
ylim([-10 10]);
set(gca,'XTick',0:0.2:4);
set(gca,'YTick',-10:4:10);
title("Tension at Point1",'fontsize',16,'fontname','times new roman');
xlabel("Time(s)",'fontname','times new roman');
ylabel("Force(N)",'fontname','times new roman');
grid on;

