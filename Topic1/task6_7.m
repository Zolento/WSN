clc;
clear all;
close all;

data = load("lab1data.txt");
x = data(:,1);
y = data(:,2);
[n,m] = size(data);

a=zeros(n);
for i=1:n
    for j=i:n
        a(i,j) = sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
    end
end

%a=a+a';
a(a==0|a>225)=inf;% 自身和超过距离阈值的节点置为inf
nb = cell(1,n);
nball = 0;
for i=1:n
    for j=i:n
        if a(i,j) ~= inf
            nb{i} = [nb{i},j];
            nb{j} = [nb{j},i];
        end  
    end
    nball = nball+size(nb{i},2); 
end
avgnball = nball/n;% 平均邻居数
SN = 3;
ST = SN;% 一开始发送消息的节点
MaxN = 6;
cnt = zeros(1,MaxN);% 记录跳数
avgdata = zeros(1,MaxN);
avgnbcom = zeros(1,MaxN);

% run gossiping
% for t=1:100
for N=1:MaxN
    figure(N);
    for i=1:n
        c = num2str(i);
        scatter(x(i),y(i),'d');
        text(x(i),y(i),c);
        hold on;
    end
    spread_target = cell(1,n);
    allflag = 0;% 停止标志：所有节点都收到了消息，置1
    nbcom = 0;% 通信的邻居数
    rdata = zeros(1,n);% 每个结点接收到的消息条数
    STsent = 0;% ST已向邻居发送
    rumorfront = ST;
    while ~allflag
        if ~STsent
            rdata(nb{ST}) = rdata(nb{ST})+1;
            alltarget = rumorfront;
            nbcom = nbcom+size(nb{ST},2);
            rdata(ST) = 1;% 需要rdata(ST) = 1来让ST持续发送
            STsent = 1;
        else
            alltarget = rumorfront;
        end
        
        for i = 1:size(alltarget,2)
            sp = nb{alltarget(i)};% 扩散目标;随机选取1个邻居
%             sp_notrcd = sp(rdata(sp)==0);% 使用这个来筛选在扩散目标中还未收到消息的邻居
%             sp_rcd = sp(rdata(sp)~=0);
%             if size(sp,2)>=N
%                 if size(sp_notrcd,2)>=N
%                     final_sp = sp_notrcd(randperm(size(sp_notrcd,2),N));
%                 else
%                     final_sp = sp_notrcd;
%                     final_sp = [final_sp, sp_rcd(randperm(size(sp_rcd,2),N-size(sp_notrcd,2)))];%随机选取N个邻居    
%                 end
%             else
%                 final_sp = sp;
%             end
            if size(sp,2)>=N
                final_sp = sp(randperm(size(sp,2),N));
            else
                final_sp = sp;
            end
            spread_target{alltarget(i)} = final_sp;
        end
        
        for i = 1:size(alltarget,2)
            for j = 1:size(spread_target{alltarget(i)},2)
                sp_target = spread_target{alltarget(i)};
                rdata(sp_target(j))=rdata(sp_target(j))+1;% 发送数据
                nbcom = nbcom+size(sp_target(j),2);% 计算发送数
                plot([x(alltarget(i)),x(sp_target(j))],[y(alltarget(i)),y(sp_target(j))],'m');%画出target点到它的邻居的路径
                hold on;
            end
        end
        tmp = zeros(1,n);
        for i=1:n
            for j = spread_target{i}
                tmp(j) = 1;
            end
        end
        spread_target = cell(1,n);
        rumorfront = find(tmp==1);
        cnt(N) = cnt(N)+1;
        allflag = size(find(rdata(2:n)==0),2)==0;
        title(sprintf('传播邻居数=%d 运行轮数=%d',N,cnt(N)));
    end
    
    rdata(ST) = rdata(ST)-1;% 把rdata(ST)还原成0
    avgdata(N) = sum(rdata)/n;
    avgnbcom(N) = nbcom/n;
end
% times(t)=avgdata(MaxN);
% cnt = zeros(1,MaxN);
% end
figure(N+1);
b = bar(1:MaxN, cnt);
title('单位时间数');
text(b.XEndPoints,b.YEndPoints,string(b.YData),'HorizontalAlignment','center','VerticalAlignment','bottom');
figure(N+2);
b = bar(1:MaxN, avgdata*n);
title('数据包转发总数')
text(b.XEndPoints,b.YEndPoints,string(b.YData),'HorizontalAlignment','center','VerticalAlignment','bottom');