clc;
clear all;
close all;
%%
r=40000/pi;
rch=r/2;% 簇首成簇半径
%%
N=load('N.mat');
N=N.N;
a=0.3;
b=1-a;
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
    N(i).credit=1; % 信誉度
    N(i).nb=[]; % Rc范围内邻居
    N(i).nbhf=[]; % Rc/2范围内邻居
    N(i).d=0; % 到SN的距离
    N(i).path=[]; % 到SN的路径
    N(i).INpath=[]; % IN路径
    N(i).r=1; % 实际接收
    N(i).t=1; % 实际发送
    N(i).ANc=0; % 申明为恶意的节点
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
for ep=1:200
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
    end
    for i=EN
        pklen=packlen;
        curnode=i;
        while(curnode~=-1) % -1代表SN
            nbs=N(curnode).nb;
            Wrp=[];
            foundIN=0; % 判断有没有找到现有的IN
            if N(curnode).steps~=1
                for j=nbs %1.找出路由节点的路径
                    if N(j).E>0&&N(j).type~=2&&N(j).ANc~=1&&size(find(N(i).path==j),2)==0 % 监督节点不作为路由,事件节点可能作路由；不能选择已经在路径的节点
                        weight=a*N(j).steps+b*Eo/N(j).E;
                        Wrp=[Wrp,weight];
                    else
                        Wrp=[Wrp,inf];
                    end
                end
                [v,idx]=min(Wrp);
                idxmins=find(Wrp==v); % 对符合最小权值的路由节点进一步筛选
                nextstepnbmax=inf;
                if size(idxmins,2)>1
                    for j=idxmins % 找出下一跳邻居最多的最小权值路由节点
                        nextstepnbtmp=0;
                        for k=N(nbs(j)).nb
                            if N(k).ANc~=1
                                if N(k).steps==N(i).steps-1
                                    nextstepnbtmp=nextstepnbtmp+1;
                                end
                                if nextstepnbtmp>nextstepnbmax
                                    idx=j;
                                    nextstepnbmax=nextstepnbtmp;
                                    v=Wrp(idx);
                                end
                            end
                        end
                    end
                end
                RP=nbs(idx);
                if N(RP).type~=0
                    N(RP).type=1;
                end
                if v<inf
                    %N(RP).E=N(RP).E-packlen*Er;
                    N(RP).type=1;
                    %N(curnode).E=N(curnode).E-packlen*Et;
                    N(i).path=[N(i).path,RP];
                    colorstr=color{mod(colorselect,colorlen)+1}; % 选择颜色
                    linetypestr=linetype{mod(linetypeselect,linetypelen)+1}; % 选择线型
                    if N(RP).steps==1
                        if mod(ep,s)==0
                            plot([N(curnode).x,N(RP).x],[N(curnode).y,N(RP).y],[linetypestr,colorstr]);
                            plot([N(RP).x,0],[N(RP).y,0],[linetypestr,colorstr]);
                            text(N(curnode).x,N(curnode).y,num2str(curnode));
                            text(N(RP).x,N(RP).y,num2str(RP));
                            %N(RP).E=N(RP).E-packlen*Et;
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
            Win=[];
            INcds=[];
            if rp==i % 开始
                targetnodes=intersect(rpnbs,N(N(i).path(2)).nb);
            elseif rp==n+1 % 到达SN
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
                    targetnodes=intersect(rpnbs,N(preIN).nb);
                end
            else
                targetnodes=intersect(rpnbs,N(preIN).nb);
            end
            [idxin,Win,INcds]=findbestIN(targetnodes,Win,INcds,N,rp,Eo,i);
            IN=INcds(idxin);
            if sum(IN)~=0
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
            else
                for j=N(preIN).nb % 从邻居里面寻找rp(有时候还是会找不到路径)
                    if N(j).type==-1&&N(j).ANc~=1
                        targetnodes=intersect(rpnbs,N(j).nb);
                        Win=[];
                        INcds=[];
                        if sum(targetnodes)~=0
                            [idxin,Win,INcds]=findbestIN(targetnodes,Win,INcds,N,rp,Eo,i);
                            IN=INcds(idxin);
                            if sum(IN)~=0
                                if mod(ep,s)==0
                                    plot([N(j).x,N(preIN).x],[N(j).y,N(preIN).y],[linetypestr,'m']);
                                    plot([N(IN).x,N(j).x],[N(IN).y,N(j).y],[linetypestr,'m']);
                                    plot([N(IN).x,N(rp).x],[N(IN).y,N(rp).y],['-','g']);
                                    text(N(IN).x,N(IN).y,num2str(IN));
                                    hold on;
                                end
                                N(i).INpath=[N(i).INpath,IN];
                                N(IN).type=2;
                                N(IN).INn=[N(IN).INn;i,rp];
                                N(rp).IN=IN;
                                N(preIN).INrp=[N(preIN).INrp;IN,j];
                                preIN=N(rp).IN;
                                preRP=rp;
                                break;
                            end
                        end
                    end
                end
%                 if sum(IN)==0
%                     
%                 end
            end
        end
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
        colorselect=colorselect+1;
        linetypeselect=linetypeselect+1;
    end
    for i=1:n
        if N(i).credit<0.8
            N(i).ANc=1;
            if N(i).type==0
                newEN=EN(EN~=i);
                EN=newEN;
            end
        end
        if N(i).type==0&&N(i).E<=0
            newEN=EN(EN~=i);
            EN=newEN;
        end
    end
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


