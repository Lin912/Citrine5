
data = importdata('newoutput1.csv');
Tck = size(data);

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



for j=1:size(data_new,2)/3
   for i=1:size(data_new,1)
       data_new(i,3*(j-1)+1)= data_new(i,3*(j-1)+1);
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0.001:0.001:10;

Y1 = data_new(:,1);
Y50 = data_new(:,148);

subplot(2,3,1);
plot(t,Y1,"-",'LineWidth',1.2);
xlim([0 10]);
ylim([-0.5 0.5]);
set(gca,'XTick',0:1:10);
set(gca,'YTick',-0.5:0.1:0.5);
title("Displacement at Point1 in Y-direction",'fontsize',16,'fontname','times new roman');
xlabel("Time(s)",'fontname','times new roman');
ylabel("Displacement(m)",'fontname','times new roman');
grid on;
hold on;


subplot(2,3,4);
plot(t,Y50,"-",'LineWidth',1.2);
xlim([0 10]);
ylim([-0.5 0.5]);
set(gca,'XTick',0:1:10);
set(gca,'YTick',-0.5:0.1:0.5);
title("Displacement at Point50 in Y-direction",'fontsize',16,'fontname','times new roman');
xlabel("Time(s)",'fontname','times new roman');
ylabel("Displacement(m)",'fontname','times new roman');
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X1 = data_new(:,2);
X50 = data_new(:,149);

subplot(2,3,2);
plot(t,X1,"-",'LineWidth',1.2);
xlim([0 10]);
ylim([-0.5 0.5]);
set(gca,'XTick',0:1:10);
set(gca,'YTick',-0.5:0.1:0.5);
title("Displacement at Point1 in X-direction",'fontsize',16,'fontname','times new roman');
xlabel("Time(s)",'fontname','times new roman');
ylabel("Displacement(m)",'fontname','times new roman');
grid on;
hold on;


subplot(2,3,5);
plot(t,X50,"-",'LineWidth',1.2);
xlim([0 10]);
ylim([-0.5 0.5]);
set(gca,'XTick',0:1:10);
set(gca,'YTick',-0.5:0.1:0.5);
title("Displacement at Point50 in X-direction",'fontsize',16,'fontname','times new roman');
xlabel("Time(s)",'fontname','times new roman');
ylabel("Displacement(m)",'fontname','times new roman');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z1 = data_new(:,3);
Z50 = data_new(:,150);

subplot(2,3,3);
plot(t,Z1,"-",'LineWidth',1.2);
xlim([0 10]);
ylim([-0.5 0.5]);
set(gca,'XTick',0:1:10);
set(gca,'YTick',-0.5:0.1:0.5);
title("Displacement at Point1 in Z-direction",'fontsize',16,'fontname','times new roman');
xlabel("Time(s)",'fontname','times new roman');
ylabel("Displacement(m)",'fontname','times new roman');
grid on;
hold on;


subplot(2,3,6);
plot(t,Z50,"-",'LineWidth',1.2);
xlim([0 10]);
ylim([-0.5 0.5]);
set(gca,'XTick',0:1:10);
set(gca,'YTick',-0.5:0.1:0.5);
title("Displacement at Point50 in Z-direction",'fontsize',16,'fontname','times new roman');
xlabel("Time(s)",'fontname','times new roman');
ylabel("Displacement(m)",'fontname','times new roman');
grid on;





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

