
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

data_new=zeros(2000,150);

for i=1:150
    for j=3:2000
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
t = 0.002:0.002:4;

Y1 = data_new(:,1);
Y50 = data_new(:,148);

subplot(2,1,1);
plot(t,Y1,"-",'LineWidth',1.2);
xlim([0 4]);
ylim([-0.2 0.2]);
set(gca,'XTick',0:0.2:4);
set(gca,'YTick',-0.2:0.1:0.2);
title("Displacement at Point1",'fontsize',16,'fontname','times new roman');
xlabel("Time(s)",'fontname','times new roman');
ylabel("Displacement(m)",'fontname','times new roman');
grid on;
hold on;





subplot(2,1,2);
plot(t,Y50,"-",'LineWidth',1.2);
xlim([0 4]);
ylim([-0.2 0.2]);
set(gca,'XTick',0:0.2:4);
set(gca,'YTick',-0.2:0.1:0.2);
title("Displacement at Point50",'fontsize',16,'fontname','times new roman');
xlabel("Time(s)",'fontname','times new roman');
ylabel("Displacement(m)",'fontname','times new roman');
grid on;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%function
function dis=v_int_Sim(t1,t2,row,data)  %辛普森,t1,t2为积分区间，row为第几列
    h=0.002;
    dis=0;
    for t=t1:1:t2-2
        int_1=data(t,row);
        int_2=data(t+1,row);
        int_3=data(t+2,row);
        dis=dis+(int_1+4*int_2+int_3)*(h/6);
    end
end

function dis=v_int_tra(t1,t2,row,data)  %梯形法，t1,t2为积分区间，row为第几列
    h=0.002;
    dis=0;
    for t=t1:1:t2-1
        int_1=data(t,row);
        int_2=data(t+1,row);
        dis=dis+(int_1+int_2)*h/2;
    end
end

