clear all;
close all;
load 'rppsummary.mat' rppsummary
lastrpp=rppsummary(:,end);
n=500;
% for i=1:n
%     if lastrpp(i)>1.5
%         lastrpp(i)=1.5;
%     end
% end
m=max(lastrpp);
X=[(0:m/n:m-m/n)',lastrpp];
scatter(X(:,1),X(:,2))
idx = dbscan(X,0.2,5);
gscatter(X(:,1),X(:,2),idx);
