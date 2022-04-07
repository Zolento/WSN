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
SN.x = x(3);
SN.y = y(3);% 汇聚节点
x = [x(1:2);x(4:n)];
y = [y(1:2);y(4:n)];
roud=[];

childmean=[];
childvar=[];
headall=[];
% for times=1:50
Eo=1; % J
Eelec=50*10^(-9); % J/bit
ETx=Eelec; % J/bit
ERx=Eelec; % J/bit
% Transmit Amplifier Types %
Eamp=10*10^(-12); % J/bit/m^2
% Data Aggregation Energy %
EDA=5*10^(-9); % J/bit

k=1024; % 每次传输的bit数

p=0.08;% 簇首占比
type = zeros(1,n-1);% 0 簇节点；1 簇首
cluster = zeros(1,n-1);% 节点归属的簇
dtch = zeros(1,n-1);% 当成为簇节点时，离它的簇首的距离
dts = zeros(1,n-1);% 当成为簇首时，离SN的距离
times = zeros(1,n-1);% 有几次被选为簇首
rleft = zeros(1,n-1);% 每个节点离下次选举剩下的轮数
op = zeros(1,n-1);% 可操作标志
E = Eo*ones(1,n-1);% 节点初始能量（J）
dead_nodes = 0;% 死亡节点数
alive = ones(1,n-1);% 0 节点已经死亡
r = 1;% 当前轮数
rmax = zeros(1,n-1);% 节点运行的最大轮数
rnd = zeros(1,n-1);
alive_nodes = [49];

while dead_nodes<n-1
    t=(p/(1-p*(mod(r,1/p))));% 阈值
    tleft=mod(r,round(1/p));
    numCLheads = 0;
    rdatach = zeros(1,n-1);% 簇首接收到的数据量
    
    for i=1:n-1
        cluster(i) = 0;
        type(i) = 0;
        if rleft(i)>0
           rleft(i)=rleft(i)-1;
        end
        if E(i)>0 && rleft(i)==0
            tmp=rand;	
            if tmp< t
                type(i) = 1;
                rleft(i)=round(1/p)-tleft;% 更新到下次选举轮数
                numCLheads=numCLheads+1;% 簇首数加1
                cluster(i)=numCLheads;% 记录下簇首序号
                dts(i)=sqrt((SN.x-x(i))^2 + (SN.y-y(i))^2);
            end
        end
    end
    
    CLheads = find(type==1);% 找出簇首集合
    
    if size(CLheads,2)~=0 % 如果选举出新的簇首，更新簇首信息,计算能量损耗
        figure(r)
%     for pp=1:10
        headall=[headall,size(CLheads,2)];
        childtmp=zeros(1,size(CLheads,2));
        if size(find(E<=0),2)>=49
            break;
        end
        for i=1:n-1
            if  type(i)==0 && E(i)>0
                dtmp = [];
                for m=1:size(CLheads,2)
                    dtmp=[dtmp,sqrt((x(CLheads(m))-x(i))^2 + (y(CLheads(m))-y(i))^2)];
                end
            [v,idx]=min(dtmp(:));% 找到离当前簇节点最近的簇首
            cluster(i)=idx; % 声明属于该簇(第几号簇)
            childtmp(idx)=childtmp(idx)+1;
            dtch(i)= v; % 更新距离
            end
        end
        if sum(childtmp)~=0
            childmean=[childmean,mean(childtmp)];
            childvar=[childvar,var(childtmp)];
        else
            alive;
        end
        for i=1:n-1
            if alive(i)==1 && type(i)==0 % 计算簇节点能量损耗
                if E(i)>0
                    if dtch(i)<=inf
                        ETx=(Eelec + Eamp*dtch(i)^2)*k;
                        E(i)=E(i)-ETx;
                    else
                        ETx=(Eelec + Eamp*dtch(i)^4)*k;
                        E(i)=E(i)-ETx;
                    end

                    if E(CLheads(cluster(i)))>0 && alive(CLheads(cluster(i)))==1 % 计算簇首能量损耗
                        ERx=(Eelec + EDA)*k;
                        E(CLheads(cluster(i)))=E(CLheads(cluster(i))) - ERx;
                        rdatach(CLheads(cluster(i))) = rdatach(CLheads(cluster(i)))+k;
                        if E(CLheads(cluster(i)))<=0  % 如果簇首在接收后死亡
                            alive(CLheads(cluster(i)))=0;
                            dead_nodes = dead_nodes +1;
                            rmax(CLheads(cluster(i)))=r;
                            dead_nodes=dead_nodes +1;
                        end
                    end
                end

                if E(i)<=0 % 如果簇节点在发送后死亡
                    alive(i) = 0;
                    dead_nodes = dead_nodes +1;
                    cluster(i) = -1;
                    rmax(i) = r;
                end
            end
        end

        for i = CLheads % 簇首与汇聚节点通信
            if alive(i)==1
                    if sqrt((x(i)-SN.x)^2+(y(i)-SN.y)^2)<=inf
                        ETx = (Eelec+EDA+Eamp*((x(i)-SN.x)^2+(y(i)-SN.y)^2))*((rdatach(i)+k)*0.3);
                        E(i)=E(i) - ETx;
                    else
                        ETx = (Eelec+EDA+Eamp*((x(i)-SN.x)^2+(y(i)-SN.y)^2)^2)*((rdatach(i)+k)*0.3);
                        E(i)=E(i) - ETx;
                    end
                if E(CLheads(cluster(i)))<=0  % 如果簇首在发射后死亡
                    alive(i) = 0;
                    rmax(i) = r;
                    dead_nodes = dead_nodes + 1;
                end
            end
        end
        r=r+1;    
%     end
                % 画图
%%%%%%%%%%%%%%%%%%
        for i=1:n-1
            if i>=3
                c = num2str(i+1);
            else
                c = num2str(i);
            end
            scatter(x(i),y(i),'b');
            text(x(i),y(i),c);
            hold on;
        end
        scatter(SN.x,SN.y,'r');
        text(SN.x,SN.y,'3');
        hold on;
        for i=1:n-1
            curCH=cluster(i);
            if alive(i)==0||type(i)==-2
                scatter(x(i),y(i),'k');
                hold on;
            end
            if type(i)==0
                plot([x(i),x(CLheads(curCH))],[y(i),y(CLheads(curCH))],'m');
                hold on;
            elseif type(i)==1
                scatter(x(i),y(i),'r');
                plot([x(i),SN.x],[y(i),SN.y],'g');
                hold on;
            end
        end
        xlim([0 1100])
        ylim([0 1100])
%%%%%%%%%%%%%%%%%%
    end
    alive_nodes = [alive_nodes, n-1-dead_nodes];
    if r==5
        break;
    end
end
% figure(2)
% plot(0:size(alive_nodes,2)-1,alive_nodes)
% ylabel("存活节点数")
% xlabel("轮数")
% xlim([0 r])
% ylim([0 n])