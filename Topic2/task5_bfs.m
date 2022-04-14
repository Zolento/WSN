clc;
clear all;
close all;
%%
r=40000/pi;
rch=r/2;% 簇首成簇半径
%%
N=load('N.mat');
N=N.N;
a=0.5;
b=0.5;
c=0.1;
p=0.08; % 簇首占比
n=500;
Rc=30; % 通信半径
Rs=10; % 感应半径
Eo=1; % 初始能量
numAN=0;
Et=0.0003; % 发送一个数据包消耗能量
Er=0.0001; % 接收一个数据包消耗能量
packlen=20; % 20个数据包
dis=zeros(n); % 距离矩阵
EN=[];
onestep=[];
er=sqrt(r*rand(1));
etheta=2*pi*rand(1);
ex=er*cos(etheta);
ey=er*sin(etheta);
%% 初始化1 需要提前固定的字段
AN=rand(1,n);
for i=1:n 
    N(i).type=-1; % 普通节点 -1 事件节点(EN 0) 中继节点（RP 1） 监督节点（IN 2)
    N(i).steps=0; % 跳数
    N(i).E=Eo; % 初始能量
    N(i).IN=0; % 该节点所属的监督节点
    N(i).INrp=[]; % 如果是监督节点，该节点用来中继到下一个监督节点的中继点(可能有一个以上,以target,rp的键值对保存)
    N(i).INn=[]; % 如果是监督节点，该节点所监督的节点(可能有一个以上,以i,rp的键值对保存)
    N(i).INnpklen=0; % 如果是监督节点，监测到的转发包数(可能有一个以上,一次只存一个)
    N(i).INpath=[]; % IN路径
    N(i).credit=1; % 信誉度
    N(i).nb=[]; % Rc范围内邻居
    N(i).nbhf=[]; % Rc/2范围内邻居
    N(i).d=0; % 到SN的距离
    N(i).path=[]; % 到SN的路径
    N(i).r=1; % 实际接收
    N(i).t=2; % 实际发送
    N(i).ANc=0; % 申明为恶意的节点
    N(i).bfspath=[]; % bfs搜索到的路径
    if AN(i)<0.05 % 恶意节点
        N(i).AN=1;
        numAN=numAN+1;
        N(i).rpp=0.7*rand; % 转发率
    else
        N(i).AN=0;
        N(i).rpp=0.8+0.2*rand; % 转发率
    end
end
%% 初始化2
for i=1:n
    N(i).d=sqrt((N(i).x)^2+(N(i).y)^2);
    for j=i:n
        dis(i,j) = sqrt((N(i).x-N(j).x)^2+(N(i).y-N(j).y)^2);
        if 0<dis(i,j)&&dis(i,j)<=Rc
            N(i).nb = [N(i).nb,j];
            N(j).nb = [N(j).nb,i];
        end
        if 0<dis(i,j)&&dis(i,j)<=Rc/2
            N(i).nbhf = [N(i).nbhf,j];
            N(j).nbhf = [N(j).nbhf,i];
        end
    end
    if (N(i).x-ex)^2+(N(i).y-ey)^2<=Rs^2
        N(i).type=0;
        N(i).E=inf; % 假设事件节点的能量很大
        EN=[EN,i];
    end
    if (N(i).x)^2+(N(i).y)^2<=Rc^2
        onestep=[onestep,i];
    end
    N(i).steps=ceil(N(i).d/Rc);
end
N(n+1).x=0;
N(n+1).y=0;
N(n+1).nb=onestep;
N(n+1).type=-1;
N(n+1).steps=0;
N(n+1).E=inf;
dis=dis+dis';
dis(dis>Rc|dis==0)=inf;
%% 寻找每一个节点的RP和IN
s=20;
for ep=1:400
    if mod(ep,s)==0
        figure(ep);
        para = [-sqrt(r), -sqrt(r), 2*sqrt(r), 2*sqrt(r)];
        rectangle('Position', para, 'Curvature', [1 1]);
        para = [ex-Rs, ey-Rs, 2*Rs, 2*Rs];
        rectangle('Position', para, 'Curvature', [1 1]);
        hold on;
        for i=1:3
            para = [-i*Rc, -i*Rc, 2*i*Rc, 2*i*Rc];
            rectangle('Position', para, 'Curvature', [1 1],'EdgeColor','r');
            hold on;
        end
        xlim([-(sqrt(r)+10) sqrt(r)+10])
        ylim([-(sqrt(r)+10) sqrt(r)+10])
        axis equal
        scatter(0,0,'r');% SN
        for i=1:n
            if N(i).type==0
                scatter(N(i).x,N(i).y,'g');
                c = num2str(i);
                text(N(i).x,N(i).y,c);
            else
                if N(i).ANc~=1
                    scatter(N(i).x,N(i).y,'b');
                end
            end
            if N(i).E<=0
                scatter(N(i).x,N(i).y,'k');
            elseif N(i).ANc==1
                scatter(N(i).x,N(i).y,'y');            
            end
        end
    end
    colorselect=0;
    color={'k','r','c','b'};
    colorlen=size(color,2);
    linetypeselect=0;
    linetype={'-','--','-.'};
    linetypelen=size(linetype,2);
    for i=1:n % 初始化
        if N(i).type~=0
            N(i).type=-1;
        end
        N(i).path=[]; % 重置path
        N(i).INpath=[];
        N(i).INnpklen=0;
        N(i).INrp=[];
        N(i).INn=[];
        N(i).IN=0;
        N(i).bfspath=[];
        N(i).flagpath=1;
    end
    for i=EN
        pklen=packlen;
        curnode=i;
        while(curnode~=-1) % -1代表SN
            nbs=N(curnode).nb;
            Wrp=[];
            foundIN=0; % 判断有没有找到现有的IN
            if N(curnode).steps~=1
                for j=nbs %结合能量和跳数，找出邻居中其他路由节点最少的路由节点的路径
                    rpnbtmp=0;
                    if N(j).E>0&&N(j).type~=0&&N(j).type~=2&&N(j).ANc~=1&&size(find(N(i).path==j),2)==0 % 监督节点\事件节点不作为路由；不能选择已经在路径的节点
                        weight=a*N(j).steps+b*Eo/N(j).E;%+0.1*rpnbtmp/size(N(j).nb,2);
                        Wrp=[Wrp,weight];
                    else
                        Wrp=[Wrp,inf];
                    end
                end
                [v,idx]=min(Wrp);
                RP=nbs(idx);
                if N(RP).type~=0
                    N(RP).type=1;
                end
                if v<inf
                    N(RP).type=1;
                    N(i).path=[N(i).path,RP];
                    colorstr=color{mod(colorselect,colorlen)+1}; % 选择颜色
                    linetypestr=linetype{mod(linetypeselect,linetypelen)+1}; % 选择线型
                    if N(RP).steps==1
                        if mod(ep,s)==0
                            plot([N(curnode).x,N(RP).x],[N(curnode).y,N(RP).y],[linetypestr,colorstr]);
                            plot([N(RP).x,0],[N(RP).y,0],[linetypestr,colorstr]);
                            text(N(curnode).x,N(curnode).y,num2str(curnode));
                            text(N(RP).x,N(RP).y,num2str(RP));
                            hold on;
                        end
                        curnode=-1;
                    else
                        if mod(ep,s)==0
                            plot([N(curnode).x,N(RP).x],[N(curnode).y,N(RP).y],[linetypestr,colorstr]);
                            text(N(curnode).x,N(curnode).y,num2str(curnode));
                            hold on;
                        end
                        curnode=RP;
                    end
                end
            else % 事件节点本身就是一跳，直接连接SN(这时的path还是保持[])
                colorstr=color{mod(colorselect,colorlen)+1}; % 选择颜色
                linetypestr=linetype{mod(linetypeselect,linetypelen)+1}; % 选择线型
                if mod(ep,s)==0
                    plot([N(curnode).x,0],[N(curnode).y,0],[linetypestr,colorstr]);
                    hold on;
                end
                curnode=-1;
            end
        end
        N(i).path=[i,N(i).path,n+1];
        preIN=i;
        preRP=i;
        for rp=N(i).path
            rpnbs=N(rp).nb;
            preINnbs=N(preIN).nb;
            Win=[];
            INcds=[];
            if rp==n+1 % 到达SN
                if sum(find(rpnbs==preIN))~=0
                    % 如果上一个IN已经是一跳节点了,说明由上个IN直接传给SN
                    % 这时INpath里显示SN的编号
                    N(i).INpath=[N(i).INpath,n+1];
                    if mod(ep,s)==0
                        plot([N(preIN).x,0],[N(preIN).y,0],[linetypestr,'m']);
                        hold on;
                    end
                    break;
                else
                    % 如果上一个IN还不是一跳节点,说明需要再找一个IN
                    % 这时INpath里面显示IN的编号，并且其INn为SN
                    targetnodes=intersect(rpnbs,N(preRP).nb);
                end
            else
                targetnodes=intersect(rpnbs,N(preRP).nb); %1.先找RP的公共邻居
            end
            if rp==i % 开始，EN不需要找IN
                % targetnodes=intersect(rpnbs,N(N(i).path(2)).nb);
                IN=i;
                N(i).INpath=[N(i).INpath,IN];
                N(IN).INn=[N(IN).INn;i,i];
                N(rp).IN=i;
            else
                [idxin,Win,INcds]=findbestIN(targetnodes,Win,INcds,N,Eo,i);
                IN=INcds(idxin); % 此处的两个RP间的IN默认存在，假如没有我也不知道怎么办
                if N(IN).ANc==1
                    IN;
                end
                if sum(N(preIN).nb==IN)~=0 % 如果IN可以跟preIN通信
                    if mod(ep,s)==0
                        plot([N(IN).x,N(preIN).x],[N(IN).y,N(preIN).y],[linetypestr,'m']);
                        plot([N(IN).x,N(rp).x],[N(IN).y,N(rp).y],['-','g']);
                        text(N(IN).x,N(IN).y,num2str(IN));
                        hold on;
                    end
                    N(i).INpath=[N(i).INpath,IN];
                    N(IN).type=2;
                    N(IN).INn=[N(IN).INn;i,rp];
                    N(rp).IN=IN;
                    preIN=N(rp).IN;
                    preRP=rp;
                else % 如果IN不能跟preIN通信(运行bfs)
                    S=[];
                    t=[];
                    for p=1:n
                        for q=p:n
                            if dis(p,q)~=inf&&(N(p).ANc==0&&N(q).ANc==0)&&(N(p).E>0&&N(q).E>0)
                                if (N(p).type==-1||N(p).type==2)&&(N(q).type==-1||N(q).type==2)
                                    S=[S,q];
                                    t=[t,p];
                                end
                            end
                        end
                    end
                    G = graph(S,t);
                    [v,E] = bfsearch(G,preIN,'edgetonew');
                    vc1=v(:,1);
                    vc2=v(:,2);
                    Win=[];
                    INcds=[];
                    targetnodes=[];
                    for p=intersect(rpnbs,N(preRP).nb)
                        if (N(p).type==-1||N(p).type==2)&&sum(find(S==p))>0
                            targetnodes=[targetnodes,p];
                            curnode=p;
                            N(p).bfspath=[N(p).bfspath,p];
                            while(curnode~=preIN)
                                for q=1:size(v,1)
                                    if curnode==vc2(q)
                                        curnode=vc1(q);
                                        N(p).bfspath=[curnode,N(p).bfspath];
                                        break;
                                    end
                                end
                            end
                        end
                    end
                    [idxin,Win,INcds]=findbestIN(targetnodes,Win,INcds,N,Eo,i);
                    idx=1;
                    if sum(idxin)==0
                        N(i).flagpath=0;
                        break; % 如果bfs还找不到，可能是网络产生分割，从preIN找不到路径了;那么就直接跳出，flagpath置0;
                    end
                    IN=INcds(idxin);
                    if N(IN).ANc==1
                        IN;
                    end
                    bfspath=N(IN).bfspath;
                    for p=bfspath
                        if p~=preIN&&p~=IN
                            N(preIN).INrp=[N(preIN).INrp;bfspath(idx+1),p];
                        end
                        if p~=IN
                            if mod(ep,s)==0
                                plot([N(bfspath(idx)).x,N(bfspath(idx+1)).x],[N(bfspath(idx)).y,N(bfspath(idx+1)).y],[linetypestr,'m']);
                                text(N(bfspath(idx)).x,N(bfspath(idx)).y,num2str(bfspath(idx)));
                                hold on;
                            end
                        end
                        if p~=i
                            N(p).type=2;
                        end
                        idx=idx+1;
                    end
                    if mod(ep,s)==0
                        plot([N(IN).x,N(rp).x],[N(IN).y,N(rp).y],['-','g']);
                        text(N(IN).x,N(IN).y,num2str(IN));
                    end
                    N(i).INpath=[N(i).INpath,IN];
                    N(IN).type=2;
                    N(IN).INn=[N(IN).INn;i,rp];
                    N(rp).IN=IN;
                    preIN=IN;
                    preRP=rp;
                end
            end
        end
        if N(i).flagpath
            % 传输数据包
            for rp=N(i).path
                if rp~=n+1 % 不是SN
                    % rp行为
                    if rp~=i % 不是起始的事件节点
                        N(rp).E=N(rp).E-pklen*Er; % 接受
                    end
                    pklen=pklen*N(rp).rpp; % 选择性转发长度
                    N(rp).E=N(rp).E-pklen*Et; % 发送给IN
                    % IN行为
                    N(N(rp).IN).E=N(N(rp).IN).E-pklen*Er; % IN接收来自rp的长度
                    N(N(rp).IN).INnpklen=pklen; % IN记录选择性转发长度
                    % rp行为
                    N(rp).E=N(rp).E-pklen*Et; % 选择性转发给下一个
                end
            end
            % IN结算转发率
            for in=N(i).INpath
                if in~=N(i).INpath(end) % 不是SN或者中继到SN的IN
                    if in~=N(i).INpath(1) % 不是起始的事件节点的IN
                        N(in).E=N(in).E-prepklen*Er; % 接受
                    else
                        prepklen=packlen; % 初始化前跳包长
                    end
                    pklen=N(in).INnpklen; % 读取监督的RP的转发包长
                    % 更新监督的RP或EN的信誉度
                    k=0;
                    for j=N(in).INn(:,1)
                        k=k+1;
                        if j==i
                            v=N(in).INn(k,2);
                            if sum(find(N(i).path==v))~=0
                                N(v).r=N(v).r+prepklen;
                                N(v).t=N(v).t+pklen;
                                N(v).credit=(N(v).t)/(N(v).r); % 更新信誉度
                            end
                        end
                    end
                    prepklen=N(in).INnpklen*N(in).rpp; % 更新监督到的包长，作为下一个RP实际接受到的包长
                    % IN之间也存在不完全转发，所以必须要用IN的rpp更新
                    N(in).E=N(in).E-prepklen*Et; % 发送给下一个IN
                else
                    N(in).E=N(in).E-prepklen*Er; % 接受
                    N(in).E=N(in).E-N(in).rpp*prepklen*Et; % 发送给SN
                    % 如果这时IN是SN，由于E是inf，上面的操作没有啥影响
                    % 如果这时IN是中继到SN的IN，单纯用rpp处理prepklen，然后转发就行了
                end
            end
        end
        colorselect=colorselect+1;
        linetypeselect=linetypeselect+1;
    end
    for i=EN
        if N(i).type~=0
            N(i).type=0; % EN可能作为rp,先复原
        end
    end
    for i=1:n
        if N(i).credit<0.5
            N(i).ANc=1;
            if N(i).type==0
                newEN=EN(EN~=i);
                EN=newEN;
            end
        end
    end
    ep
    if sum(EN)==0
        break;
    end
end
%% 计算误检率
numANc=0;
numFANc=0;
for i=1:n
    if N(i).ANc
        numANc=numANc+1;
        if ~N(i).AN
            numFANc=numFANc+1;
        end
    end
end
error=numFANc/numANc;
error


