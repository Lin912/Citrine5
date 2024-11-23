clear
clc
data=importdata('output.csv');

if mod(size(data,2),10)==0
    flag=1;
else
    flag=0;
end
if flag==1
    for i=1:size(data,1)
        for j=1:size(data,2)/10
            data(i,10*(j-1)+1:10*(j-1)+3)=data(i,10*(j-1)+1:10*(j-1)+3)*inv( ...
                [cos(data(i,10*(j-1)+7))*cos(data(i,10*(j-1)+8)),cos(data(i,10*(j-1)+7))*sin(data(i,10*(j-1)+8)),-sin(data(i,10*(j-1)+7)); ...
                -sin(data(i,10*(j-1)+8)),cos(data(i,10*(j-1)+8)),0; ...
                cos(data(i,10*(j-1)+8))*sin(data(i,10*(j-1)+7)),sin(data(i,10*(j-1)+8))*sin(data(i,10*(j-1)+7)),cos(data(i,10*(j-1)+7))]);
        end
    end
else
end
dlmwrite('newoutput1.csv',data);
data;