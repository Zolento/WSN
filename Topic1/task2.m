clc;
clear all;
close all;

data = load("lab1data.txt");
x = data(:,1);
y = data(:,2);
n = size(data,1);%行数=个数
num = [];
nodeflag = zeros(1,n);%某位置1表示通信半径内存在其他传感器节点
allflag = 0;%1表示所有节点的通信半径内存在其他传感器节点
f = 0.01;
T = 10-f;
last=zeros(1,2);
while(~allflag)
    T = T + f;%T先加，初始值是10
    for i=1:n-1
        L = find(nodeflag==0);%找到nodeflag还为0的节点比较
        if size(L,2) == 0%如果没有为0的节点了，跳出
            allflag = 1;
            break;
        else
            for p = L%x(i)和x(p)依次比较
                if sqrt((x(i)-x(p))^2+(y(i)-y(p))^2)<T && i~=p%不和自身比较
                    last=[i,p];
                    nodeflag(i) = 1;
                    nodeflag(p) = 1;
                end
            end
        end
    end
    L = find(nodeflag==0);
    num = [num, size(L,2)];
end
d_true=sqrt((x(last(1))-x(last(2)))^2+(y(last(1))-y(last(2)))^2);
plot(10-f:f:T, num);
title("节点通信半径取值与节点数目间的关系");
xlabel("T");
ylabel("num");
