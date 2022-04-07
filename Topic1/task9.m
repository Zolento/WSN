clc;
clear all;
close all;

data = load("lab1data.txt");
xall = data(:,1);
yall = data(:,2);

[n,m] = size(data);
r = 1;% 当前轮数
Rc = 225;% 最大通信半径
SN.x = xall(3);
SN.y = yall(3);% 汇聚节点
x = [xall(1:2);xall(4:n)];
y = [yall(1:2);yall(4:n)];

a=zeros(n-1);
for i=1:n-1
    for j=i:n-1
        a(i,j) = sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
    end
end
a=a+a';
nb_Rc = cell(1,n-1);
nb_2Rc = cell(1,n-1);
for i=1:n-1
    for j=i:n-1
        if 0<a(i,j)&&a(i,j)<=Rc
            nb_Rc{i} = [nb_Rc{i},j];
            nb_Rc{j} = [nb_Rc{j},i];
        elseif Rc<a(i,j)&&a(i,j)<=2*Rc
            nb_2Rc{i} = [nb_2Rc{i},j];
            nb_2Rc{j} = [nb_2Rc{j},i];
        end
    end
end

Eo=1; % J
Eelec=50*10^(-9); % J/bit
ETx=Eelec; % J/bit
ERx=Eelec; % J/bit
% Transmit Amplifier Types %
Eamp=10*10^(-12); % J/bit/m^2
% Data Aggregation Energy %
EDA=5*10^(-9); % J/bit

k=1024; % 每次传输的bit数
firstdead=0;
p = 0.1;% 簇首占比
type = -1*ones(1,n-1);% 0 候选节点；1 候选簇头；2 簇节点；3 簇头；-1 普通节点；-2 死亡
cluster = zeros(1,n-1);% 节点归属的簇
chid = zeros(1,n-1);% 簇首编号
dtch = zeros(1,n-1);% 当成为簇节点时，离它的簇首的距离
dts = zeros(1,n-1);% 当成为簇首时，离SN的距离
times = zeros(1,n-1);% 有几次被选为簇首
rleft = zeros(1,n-1);% 每个节点离下次选举剩下的轮数
op = zeros(1,n-1);% 可操作标志
E = Eo*ones(1,n-1);% 节点初始能量（J）
dead_nodes = 0;% 死亡节点数
alive = ones(1,n-1);% 0 节点已经死亡
rmax = zeros(1,n-1);% 节点运行的最大轮数
rnd = zeros(1,n-1);
alive_nodes = [];


while dead_nodes<n-1
    t=(p/(1-p*(mod(r,1/p))));% 阈值
    tleft=mod(r,round(1/p));
    numCLheads = 0;
    type(type~=-2)=-1;
    rdatach = zeros(1,n-1);% 簇首接收到的数据量
    start = 1;% 新一轮选举
    CLheads = [];
    
    while size(find(type==-1),2)>0||size(find(type==0),2)>0
        chtarget = [];
        if start
            target = find(type==-1);
        else
            target = find(type==0);
        end
        if size(find(type==0),2)==0
            target = find(type==-1);
        end
        for i=target % step1 选出候选簇头
            if rleft(i)>0
               rleft(i)=rleft(i)-1;
            end
            if E(i)>0 && rleft(i)==0
                tmp=rand;
                if tmp<t*E(i)/Eo*size(nb_Rc{i})/7.44
                    type(i) = 1;
                    rleft(i)=round(1/p)-tleft;% 更新到下次选举轮数
                    chtarget = [chtarget,i];% 本轮候选簇头集
                    numCLheads = numCLheads+1;
                    chid(i) = numCLheads;
                    start=0;
%                     numCLheads=numCLheads+1;% 簇首数加1
%                     cluster(i)=numCLheads;% 记录下簇首序号
%                     dts(i)=sqrt((SN.x-x(i))^2 + (SN.y-y(i))^2);
                end
            end
        end
        l = size(chtarget,2);
        conflict = zeros(l);
        del = [];
        for i=1:l% step2 <Rc范围 簇头竞争
            for j=i:l
                if size(find(nb_Rc{chtarget(i)}==chtarget(j)),2) ~= 0
                    conflict(i,j)=1;
                end
            end
        end
        while sum(conflict(:))~=0
            [row,col]=find(conflict==1);
            if E(chtarget(row(1)))>=E(chtarget(col(1)))
                conflict(:,col(1))=0;
                del = [del,col(1)];
                type(chtarget(col(1)))=-1;
            else
                conflict(row(1),:)=0;
                del = [del,row(1)];
                type(chtarget(row(1)))=-1;
            end
        end
        for i=1:l
            if sum(i==del)==0
                CLheads = [CLheads, chtarget(i)];
            end
        end
        type(CLheads)=3;
        for i=CLheads
            type(nb_2Rc{i}(type(nb_2Rc{i})==-1))=0;% Rc-2Rc的普通节点一律声明为候选节点
        end
        for i=CLheads
            cluster(nb_Rc{i})=i;
            type(nb_Rc{i}(type(nb_Rc{i})~=-2))=2;
        end
    end
    
    if size(CLheads,2)~=0 % 如果选举出新的簇首，更新簇首信息,计算能量损耗
        figure(r)
        RP=cell(1,n-1);% 中继信息
        RPplan=cell(1,n-1);% 中继方案
        s = size(CLheads,2);
        path=1*ones(1,s);
        d(1:s)=inf;
        SNs=-1;
        chd=zeros(s); % 计算簇首路径
        found=zeros(1,s);
        pre_node=-1*ones(1,s);
        for i=1:s
            if  CLheads(i)==12 && SNs==-1
                SNs=i;
            end
            for j=i+1:s
                for tmp=nb_Rc{CLheads(i)}
                    if  tmp==12 && SNs==-1
                        SNs=i;
                    end
                    if type(tmp)==2
                       for t=nb_Rc{tmp}
                           if cluster(t)==CLheads(j)&&alive(t)&&alive(tmp)
                               chd(i,j)=1;
                               chd(j,i)=1;
                               RP{CLheads(i)}=[RP{CLheads(i)},{[CLheads(j),tmp,t]}];
                               RP{CLheads(j)}=[RP{CLheads(j)},{[CLheads(i),t,tmp]}];
                           end
                       end
                    end
                end
            end
        end
        for tmp=nb_Rc{CLheads(i)}
            if  tmp==12 && SNs==-1
                SNs=i;
            end
        end
        if SNs==-1
            break;
        end
        root=SNs;
        found(SNs)=1;
        d(SNs)=0;
        chd(chd==0)=inf;
        if sum(alive==0)~=0
            pause(0.1);
        end
        % run dijkstra
        while sum(found)<s  %看是否所有的点都标记为P标号
            target=find(found==0); 
            dtmp=d;
            d(target)=min(d(target),d(SNs)+chd(SNs,target));%对还没有存入found的点更新最短通路
            minp_nodeidx=find(d(target)==min(d(target)));
            updated=find(d(target)~=dtmp(target));%最短路径被更新过的节点
            if size(updated,2)~=0
                pre_node(target(updated))=SNs;
            end
            SNs=target(minp_nodeidx(1));   %将当前符合最短通路的节点设为源点
            found(SNs)=1;%已找到最小路径 found(i)=1
            path(SNs)=pre_node(SNs);
        end
        
        if size(find(path==-1),2)% 网络产生分割
            break;
        end
        
        for i=1:s % 求出最短RP方案,画出路径
            pre=CLheads(path(i));
            cur=CLheads(i);
            mind=inf;
            idx=0;
            pos=0;
            if i~=root
                for j=RP{cur}
                    pos=pos+1;
                    if j{1}(1)==pre
                        d1=a(cur,j{1}(2));
                        d2=a(j{1}(2),j{1}(3));
                        d3=a(j{1}(3),pre);
                        if mind>d1+d2+d3
                            mind=d1+d2+d3;
                            idx=pos;
                        end
                    end
                end
                RPplan(cur)=RP{cur}(idx);
                plot([x(cur),x(RPplan{cur}(2))],[y(cur),y(RPplan{cur}(2))],'g');
                hold on;
                plot([x(RPplan{cur}(2)),x(RPplan{cur}(3))],[y(RPplan{cur}(2)),y(RPplan{cur}(3))],'g');
                hold on;
                plot([x(RPplan{cur}(3)),x(pre)],[y(RPplan{cur}(3)),y(pre)],'g');
                hold on;
            else
                if cur==12
                    plot([x(cur),SN.x],[y(cur),SN.y],'g');
                    hold on;
                else
                    plot([x(cur),x(12)],[y(cur),y(12)],'g');
                    hold on;
                    plot([x(12),SN.x],[y(12),SN.y],'g');
                    hold on;
                end
                
            end
        end
        
        for i=1:n-1
            if alive(i)==1 && type(i)==2 % 计算簇节点能量损耗
                curCH=cluster(i);
                dtch(i) = a(i,curCH);
                ETx=(Eelec + Eamp*dtch(i)^2)*k;
                E(i)=E(i)-ETx;
                ERx=(Eelec + EDA)*k;
                E(curCH)=E(curCH) - ERx;
                rdatach(curCH) = rdatach(curCH)+k;
            end
        end
        
        while sum(rdatach(CLheads))~=0
            for i=CLheads % 簇首与汇聚节点通信
                if alive(i)&&i~=CLheads(root)&&rdatach(i)~=0
                    nextCH=RPplan{i}(1);
                    RP1=RPplan{i}(2);
                    RP2=RPplan{i}(3);
                    
                    RPpack=rdatach(i)*0.3;
                    ETx=(Eelec+EDA+Eamp*((x(i)-x(RP1))^2+(y(i)-y(RP1))^2))*RPpack;%发射给中继节点1
                    E(i)=E(i)-ETx;
                    rdatach(i)=0;
                    
                    ERx=(Eelec + EDA)*RPpack;% 中继节点1接收
                    ETx=(Eelec + Eamp*a(RP1,RP2)^2)*RPpack;% 中继节点1发射给2
                    E(RP1)=E(RP1)-ETx-ERx;
                    ERx=(Eelec + EDA)*RPpack;
                    ETx=(Eelec + Eamp*a(RP2,i)^2)*RPpack;% 中继节点2发射给下一个簇首
                    E(RP2)=E(RP2)-ETx-ERx;
                    
                    ERx=(Eelec + EDA)*RPpack;
                    E(nextCH)=E(nextCH)-ERx;
                    rdatach(nextCH)=rdatach(nextCH)+RPpack;
                end
            end
            if find(rdatach~=0)==CLheads(root)%已经全部传到最近的簇首了
                break;
            end
        end
        
        dead_new=find(E<=0&alive==1);
        for i=dead_new %更新死亡节点(放最后显然是不科学的，但是传输途中更新太麻烦了，不想写备用路径方案)
            alive(i)=0;
            rmax(i)=r;
            type(i)=-2;
            dead_nodes = dead_nodes + 1;
        end
        
        if dead_nodes>0 && firstdead==0
            firstdead=r;
        end
        
        alive_nodes = [alive_nodes, n-1-dead_nodes];
        % 画图
%%%%%%%%%%%%%%%%%%
        
        title(num2str(s-1))
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
        para = [SN.x-Rc, SN.y-Rc, 2*Rc, 2*Rc];
        rectangle('Position', para, 'Curvature', [1 1]);
        text(xall(3),yall(3),'3');
        hold on;
        for i=1:n-1
            curCH=cluster(i);
            if alive(i)==0||type(i)==-2
                scatter(x(i),y(i),'k');
                hold on;
            end
            if type(i)==2
                plot([x(i),x(curCH)],[y(i),y(curCH)],'m');
                hold on;
            elseif type(i)==3
                scatter(x(i),y(i),'r');
%                 plot([x(i),SN.x],[y(i),SN.y],'g');
                para = [x(i)-Rc, y(i)-Rc, 2*Rc, 2*Rc];
                rectangle('Position', para, 'Curvature', [1 1]);
                hold on;
            end
        end
        xlim([-200 1200])
        ylim([-200 1200])
%%%%%%%%%%%%%%%%%%
    else
        r = r-1;
    end
    % e = e+1;
    r=r+1;
    if r == 5
        break;
    end
end
% plot(1:size(alive_nodes,2),alive_nodes)%存活节点数曲线
% xlim([0 r])
% ylim([0 n])