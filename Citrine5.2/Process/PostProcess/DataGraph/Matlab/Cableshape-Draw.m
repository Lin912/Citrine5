%%%%%%%%%%%%%%%%
%DataInput
data = importdata('newoutput1.csv');
Tck = size(data);
%DataIntProcess
if mod(size(data,2),10)==0
    flag=1;
else
    flag=0;
end
if flag==1
    for j=1:size(data,2)/10
        for i=1:size(data,1)
            Tck(i,3*(j-1)+1:3*(j-1)+3)= data(i,10*(j-1)+1:10*(j-1)+3);
        end
    end
else
end

data_new=zeros(10000,150);

for i=1:150
    for j=3:10000
        data_new(j,i)=v_int_Sim(1,j,i,Tck);
    end
end

for i=1:150
    for j=2:3
        data_new(j,i)=v_int_tra(1,j,i,Tck);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


K = 0:-0.13:-6.5;

for j=1:size(data_new,2)/3
   for i=1:size(data_new,1)
       data_new(i,3*(j-1)+1)= data_new(i,3*(j-1)+1) + K(j);
   end
end

%Row
Xr0 = zeros(50,1);
Yr0 = zeros(50,1);
Zr0 = zeros(50,1);

Xr1 = zeros(50,1);
Yr1 = zeros(50,1);
Zr1 = zeros(50,1);

Xr2 = zeros(50,1);
Yr2 = zeros(50,1);
Zr2 = zeros(50,1);

Xr3 = zeros(50,1);
Yr3 = zeros(50,1);
Zr3 = zeros(50,1);

Xr4 = zeros(50,1);
Yr4 = zeros(50,1);
Zr4 = zeros(50,1);

Xr5 = zeros(50,1);
Yr5 = zeros(50,1);
Zr5 = zeros(50,1);

Xr6 = zeros(50,1);
Yr6 = zeros(50,1);
Zr6 = zeros(50,1);

Xr7 = zeros(50,1);
Yr7 = zeros(50,1);
Zr7 = zeros(50,1);

Xr8 = zeros(50,1);
Yr8 = zeros(50,1);
Zr8 = zeros(50,1);

Xr9 = zeros(50,1);
Yr9 = zeros(50,1);
Zr9 = zeros(50,1);

Xr10 = zeros(50,1);
Yr10 = zeros(50,1);
Zr10 = zeros(50,1);

Xr11 = zeros(50,1);
Yr11 = zeros(50,1);
Zr11 = zeros(50,1);

Xr12 = zeros(50,1);
Yr12 = zeros(50,1);
Zr12 = zeros(50,1);

Xr13 = zeros(50,1);
Yr13 = zeros(50,1);
Zr13 = zeros(50,1);

Xr14 = zeros(50,1);
Yr14 = zeros(50,1);
Zr14 = zeros(50,1);

Xr15 = zeros(50,1);
Yr15 = zeros(50,1);
Zr15 = zeros(50,1);

Xr16 = zeros(50,1);
Yr16 = zeros(50,1);
Zr16 = zeros(50,1);

Xr17 = zeros(50,1);
Yr17 = zeros(50,1);
Zr17 = zeros(50,1);

Xr18 = zeros(50,1);
Yr18 = zeros(50,1);
Zr18 = zeros(50,1);

Xr19 = zeros(50,1);
Yr19 = zeros(50,1);
Zr19 = zeros(50,1);

Xr20 = zeros(50,1);
Yr20 = zeros(50,1);
Zr20 = zeros(50,1);

for j=1:size(data_new,2)/3 
       Xr0(j) = data_new(1,3*(j-1)+1);
       Yr0(j) = data_new(1,3*(j-1)+2);
       Zr0(j) = data_new(1,3*(j-1)+3);
end

for j=1:size(data_new,2)/3 
       Xr1(j) = data_new(100,3*(j-1)+1);
       Yr1(j) = data_new(100,3*(j-1)+2);
       Zr1(j) = data_new(100,3*(j-1)+3);
end

for j=1:size(data_new,2)/3 
       Xr2(j) = data_new(200,3*(j-1)+1);
       Yr2(j) = data_new(200,3*(j-1)+2);
       Zr2(j) = data_new(200,3*(j-1)+3);
end

for j=1:size(data_new,2)/3 
       Xr3(j) = data_new(300,3*(j-1)+1);
       Yr3(j) = data_new(300,3*(j-1)+2);
       Zr3(j) = data_new(300,3*(j-1)+3);
end

for j=1:size(data_new,2)/3 
       Xr4(j) = data_new(400,3*(j-1)+1);
       Yr4(j) = data_new(400,3*(j-1)+2);
       Zr4(j) = data_new(400,3*(j-1)+3);
end

for j=1:size(data_new,2)/3 
       Xr5(j) = data_new(500,3*(j-1)+1);
       Yr5(j) = data_new(500,3*(j-1)+2);
       Zr5(j) = data_new(500,3*(j-1)+3);
end

for j=1:size(data_new,2)/3 
       Xr6(j) = data_new(600,3*(j-1)+1);
       Yr6(j) = data_new(600,3*(j-1)+2);
       Zr6(j) = data_new(600,3*(j-1)+3);
end

for j=1:size(data_new,2)/3 
       Xr7(j) = data_new(700,3*(j-1)+1);
       Yr7(j) = data_new(700,3*(j-1)+2);
       Zr7(j) = data_new(700,3*(j-1)+3);
end

for j=1:size(data_new,2)/3 
       Xr8(j) = data_new(800,3*(j-1)+1);
       Yr8(j) = data_new(800,3*(j-1)+2);
       Zr8(j) = data_new(800,3*(j-1)+3);
end

for j=1:size(data_new,2)/3 
       Xr9(j) = data_new(900,3*(j-1)+1);
       Yr9(j) = data_new(900,3*(j-1)+2);
       Zr9(j) = data_new(900,3*(j-1)+3);
end

for j=1:size(data_new,2)/3 
       Xr10(j) = data_new(1000,3*(j-1)+1);
       Yr10(j) = data_new(1000,3*(j-1)+2);
       Zr10(j) = data_new(1000,3*(j-1)+3);
end

for j=1:size(data_new,2)/3 
       Xr11(j) = data_new(1100,3*(j-1)+1);
       Yr11(j) = data_new(1100,3*(j-1)+2);
       Zr11(j) = data_new(1100,3*(j-1)+3);
end

for j=1:size(data_new,2)/3 
       Xr12(j) = data_new(1200,3*(j-1)+1);
       Yr12(j) = data_new(1200,3*(j-1)+2);
       Zr12(j) = data_new(1200,3*(j-1)+3);
end

for j=1:size(data_new,2)/3 
       Xr13(j) = data_new(1300,3*(j-1)+1);
       Yr13(j) = data_new(1300,3*(j-1)+2);
       Zr13(j) = data_new(1300,3*(j-1)+3);
end

for j=1:size(data_new,2)/3 
       Xr14(j) = data_new(1400,3*(j-1)+1);
       Yr14(j) = data_new(1400,3*(j-1)+2);
       Zr14(j) = data_new(1400,3*(j-1)+3);
end

for j=1:size(data_new,2)/3 
       Xr15(j) = data_new(1500,3*(j-1)+1);
       Yr15(j) = data_new(1500,3*(j-1)+2);
       Zr15(j) = data_new(1500,3*(j-1)+3);
end

for j=1:size(data_new,2)/3 
       Xr16(j) = data_new(1600,3*(j-1)+1);
       Yr16(j) = data_new(1600,3*(j-1)+2);
       Zr16(j) = data_new(1600,3*(j-1)+3);
end

for j=1:size(data_new,2)/3 
       Xr17(j) = data_new(1700,3*(j-1)+1);
       Yr17(j) = data_new(1700,3*(j-1)+2);
       Zr17(j) = data_new(1700,3*(j-1)+3);
end

for j=1:size(data_new,2)/3 
       Xr18(j) = data_new(1800,3*(j-1)+1);
       Yr18(j) = data_new(1800,3*(j-1)+2);
       Zr18(j) = data_new(1800,3*(j-1)+3);
end

for j=1:size(data_new,2)/3 
       Xr19(j) = data_new(1900,3*(j-1)+1);
       Yr19(j) = data_new(1900,3*(j-1)+2);
       Zr19(j) = data_new(1900,3*(j-1)+3);
end

for j=1:size(data_new,2)/3 
       Xr20(j) = data_new(2000,3*(j-1)+1);
       Yr20(j) = data_new(2000,3*(j-1)+2);
       Zr20(j) = data_new(2000,3*(j-1)+3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Col
Xc = zeros(1000,1);
Yc = zeros(1000,1);
Zc = zeros(1000,1);
for j=1:size(data_new,1) 
       Xc(j) = data_new(j,1);
       Yc(j) = data_new(j,2);
       Zc(j) = data_new(j,3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot
plot3(Xr0,Yr0,Zr0);
hold on;
plot3(Xr1,Yr1,Zr1);
hold on;
plot3(Xr2,Yr2,Zr2);
hold on;
plot3(Xr3,Yr3,Zr3);
hold on;
plot3(Xr4,Yr4,Zr4);
hold on;
plot3(Xr5,Yr5,Zr5);
hold on;
plot3(Xr6,Yr6,Zr6);
hold on;
plot3(Xr7,Yr7,Zr7);
hold on;
plot3(Xr8,Yr8,Zr8);
hold on;
plot3(Xr9,Yr9,Zr9);
hold on;
plot3(Xr10,Yr10,Zr10);
hold on;
plot3(Xr11,Yr11,Zr11);
hold on;
plot3(Xr12,Yr12,Zr12);
hold on;
plot3(Xr13,Yr13,Zr13);
hold on;
plot3(Xr14,Yr14,Zr14);
hold on;
plot3(Xr15,Yr15,Zr15);
hold on;
plot3(Xr16,Yr16,Zr16);
hold on;
plot3(Xr17,Yr17,Zr17);
hold on;
plot3(Xr18,Yr18,Zr18);
hold on;
plot3(Xr19,Yr19,Zr19);
hold on;
plot3(Xr20,Yr20,Zr20);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%图例设置
xlim([-10 1]);
ylim([-5 5]);
zlim([-5 5]);
set(gca,'XTick',-10:1:1);
set(gca,'YTick',-5:1:5);
set(gca,'ZTick',-5:1:5);
hold on;

%顶部
p = plot3(Xc,Yc,Zc,":");
p.Color = "#A2142F";
p.LineWidth = 1.5;

%标题坐标轴标注
title("Profile in 3-Dimensional",'fontsize',16,'fontname','times new roman');
xlabel('X','fontname','times new roman');
ylabel('Y','fontname','times new roman');
zlabel('Z','fontname','times new roman');
view(60,30);
grid on;
axis square;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%function
function dis=v_int_Sim(t1,t2,row,data)  %辛普森,t1,t2为积分区间，row为第几列
    h=0.001;
    dis=0;
    for t=t1:1:t2-2
        int_1=data(t,row);
        int_2=data(t+1,row);
        int_3=data(t+2,row);
        dis=dis+(int_1+4*int_2+int_3)*(h/6);
    end
end

function dis=v_int_tra(t1,t2,row,data)  %梯形法，t1,t2为积分区间，row为第几列
    h=0.001;
    dis=0;
    for t=t1:1:t2-1
        int_1=data(t,row);
        int_2=data(t+1,row);
        dis=dis+(int_1+int_2)*h/2;
    end
end

