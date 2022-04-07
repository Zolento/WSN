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
allflag = 0;% 停止标志：所有节点都收到了消息，置1
cnt = 0;% 记录跳数
nbcom = 0;% 通信的邻居数
SN = 3;
ST = SN;% 一开始发送消息的节点


while ~allflag
    pause(0.1);
    if ~sent(ST)
        %rdata(nb{ST}) = rdata(nb{ST})+1;
        alltarget = ST;
        sent(ST) = 1;
    else
        alltarget = find(rdata~=0&sent~=1);
        sent(alltarget) = 1;% 标志为已发送
    end
    for i = 1:size(alltarget,2)
        spread_target = nb{alltarget(i)};% 扩散目标;找出邻居
        spread_target = spread_target(spread_target~=prenode(alltarget(i)));% 在扩散目标中去除当前目标的前序节点
        spread_target = spread_target(rdata(spread_target)==0);% 在扩散目标中去除已收到消息的节点
        prenode(spread_target(rdata(spread_target)==0)) = alltarget(i);% 更新扩散目标的前序节点
        rdata(spread_target)=rdata(spread_target)+1;% 发送数据
        nbcom = nbcom+size(spread_target,2);% 计算发送数
        for j = 1:size(spread_target,2)
            plot([x(alltarget(i)),x(spread_target(j))],[y(alltarget(i)),y(spread_target(j))],'m');%画出target点到它的邻居的路径
            hold on;
        end
    end
    cnt = cnt+1;
    allflag = size(find(rdata==0),2)==1;% 最后一个收到的节点不会回传给邻居了
end
avgdata = sum(rdata)/(n-1);
avgnbcom = nbcom/n;