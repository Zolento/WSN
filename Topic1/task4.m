clc;
clear all;
close all;
data = load("lab1data.txt");
x = data(:,1);
y = data(:,2);
[n,m] = size(data);
for i=1:n
    c = num2str(i);
    scatter(x(i),y(i),'d');
    text(x(i),y(i),c);
    hold on;
end

SN = 1;
RM = 2;
mind = x(1)^2+y(1)^2;
maxd = (x(RM)-x(SN))^2+(y(RM)-y(SN))^2;
for i=2:n
    if mind>x(i)^2+y(i)^2
        mind = x(i)^2+y(i)^2;
        SN = i;
    end
end
for i=2:n
    if maxd<(x(i)-x(SN))^2+(y(i)-y(SN))^2
        maxd = (x(i)-x(SN))^2+(y(i)-y(SN))^2;
        RM = i;
    end
end

scatter(x(SN),y(SN), 'r')
scatter(x(RM),y(RM), 'g')
hold on;
a=zeros(n);
for i=1:n
    for j=i:n
        a(i,j) = sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
    end
end

a=a+a';
a(a==0|a>222.1734)=inf; % 自身和超过距离阈值的节点置为inf
a(a~=inf)=1;% 路径全部置1，路径长度等同跳数
SNs=SN;  %SN 查找的源点
found=zeros(1,n);
found(SNs)=1;  %已经求出最小跳数，found置1
found_order=SNs; %存放顺序
path=SNs*ones(1,n); %最短通路中第i节点 前序节点的序号
pre_node=-1*ones(1,n);
d(1:n)=inf;
d(SNs)=0;  %i节点到源点最短通路的值
% run dijkstra
while sum(found)<n  %看是否所有的点都标记为P标号
    target=find(found==0); 
    dtmp=d;
    d(target)=min(d(target),d(SNs)+a(SNs,target));%对还没有存入found的点更新最短通路
    
    minp_nodeidx=find(d(target)==min(d(target)));
    updated=find(d(target)~=dtmp(target));%最短路径被更新过的节点
    if size(updated,2)~=0
        pre_node(target(updated))=SNs;
    end
    SNs=target(minp_nodeidx(1));   %将当前符合最短通路的节点设为源点
    %可能会有多个符合条件的节点，选择其中一个
    target(minp_nodeidx)
    pre_node(target(minp_nodeidx))
    found(SNs)=1;%已找到最小路径 found(i)=1
    found_order=[found_order,SNs];  %更新路径查找顺序
    path(SNs)=pre_node(SNs);
    path(SNs)
    plot([x(SNs),x(path(SNs))],[y(SNs),y(path(SNs))],'m');
    hold on;
end

length = [];
m=0;
route=cell(1,n-1);
for i=1:n
    trace = [];
    curnode = i;
    while(1)
        trace = [trace, curnode];
        if curnode==SN
            break;
        end
        curnode = path(curnode);
    end
    route{i}=trace;
    for i=1:size(trace,2)-1
       plot([x(trace(i)),x(trace(i+1))],[y(trace(i)),y(trace(i+1))],'m'); 
    end
    length=[length,size(trace,2)-1]; 
end

avglength = sum(length)/(n-1);
%avglength = sum(d/(n-1));
