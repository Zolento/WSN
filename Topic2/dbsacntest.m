clear all;
close all;
load 'rppsummary.mat' rppsummary
load 'NwithResults.mat' N
%% 记录mexep轮后的rpp
% 就是rppsummary的最后一列
n=500;
maxep=200; % size(rppsummary,1);
lastrpp=rppsummary(:,maxep);
figure(maxep+1)
title('记录的累计转发率分布')
for i=1:n
    if N(i).AN
        scatter(i,lastrpp(i),'r')
        hold on
    else
        scatter(i,lastrpp(i),'b')
        hold on
    end
    if N(i).isabnormal
        scatter(i,lastrpp(i),'k','+')
        hold on
    end
end
figure(maxep+2)
title('理想的转发率平均值分布')
for i=1:n
    if N(i).AN
        scatter(i,N(i).rppavg,'r')
        hold on
    else
        scatter(i,N(i).rppavg,'b')
        hold on
    end
    if N(i).isabnormal
        scatter(i,N(i).rppavg,'k','+')
        hold on
    end
end
%% dbscan聚类
% for i=1:n
%     if lastrpp(i)>1
%         lastrpp(i)=1;
%     end
% end
figure(maxep+3)
title('聚类结果')
m=max(lastrpp);
X=[(0:m/n:m-m/n)',lastrpp];
scatter(X(:,1),X(:,2))
idx = dbscan(X,0.1,7);
gscatter(X(:,1),X(:,2),idx);
%% 计算误检率
numANc=0;
numFANc=0;
for i=1:n
    N(i).ANc=idx(i);
    if N(i).ANc==-1&&N(i).ispath
        numANc=numANc+1;
        if ~N(i).AN
            numFANc=numFANc+1;
        end
    end
end
error=numFANc/numANc;
fprintf('检测正确率:%f',1-error);