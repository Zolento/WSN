clc;
clear all;

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
avgnball = nball/n;

rdata = zeros(1,n);% 每个结点接收到的消息条数
sent = zeros(1,n);% 已向邻居发送过数据的标识
prenode = zeros(1,n);% 前序节点，防止回传(同跳节点间依旧存在回传)
spread_target = cell(1,n);
allflag = 0;% 停止标志：所有节点都收到了消息，置1
cnt = 0;% 记录跳数
nbcom = 0;% 通信的邻居数
SN = 3;
ST = SN;% 一开始发送消息的节点
N = 1;% 随机选取的邻居个数

% run gossiping
while ~allflag
    if ~sent(ST)
        rdata(nb{ST}) = rdata(nb{ST})+1;
        alltarget = ST;
        rdata(ST) = 1;% 需要rdata(ST) = 1来让ST继续传输
        sent(ST) = 1;
    else
        alltarget = find(rdata~=0);
        %sent(alltarget) = 1;% 标志为已发送
    end
    for i = 1:size(alltarget,2)
        sp = nb{alltarget(i)};% 扩散目标;随机选取1个邻居
        % sp = sp(rdata(sp)==0);% 在扩散目标中已经收到消息的邻居
        if size(sp,2)>=N
            sp = sp(randperm(size(sp,2),N));%随机选取N个邻居
        end
        spread_target{alltarget(i)} = sp;
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
    cnt = cnt+1;
    allflag = size(find(rdata==0),2)==0;% 最后一个收到的节点不会回传给邻居了;和v2有点小区别
end
rdata(ST) = rdata(ST)-1;% 把rdata(ST)还原成0
avgdata = sum(rdata)/n;
avgnbcom = nbcom/50;

