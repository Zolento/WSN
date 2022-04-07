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
Rc=20; % 通信半径
Eo=1; % 初始能量
Et=0.0003; % 发送能量
Er=0.0001; % 接收能量
dis=zeros(n); % 距离矩阵
EN=[];
er=sqrt(r*rand(1));
etheta=2*pi*rand(1);
ex=er*cos(etheta);
ey=er*sin(etheta);
%% 初始化1 需要提前固定的字段
for i=1:n 
    N(i).type=-1; % 普通节点 -1 事件节点(EN 0) 中继节点（RP 1） 监督节点（IN 2)
    N(i).steps=0; % 跳数
    N(i).E=Eo; % 初始能量
    N(i).credit=1; % 信誉度
    N(i).nb=[]; % Rc范围内邻居
    N(i).nbhf=[]; % Rc/2范围内邻居
    N(i).d=0; % 到SN的距离
    N(i).path=[]; % 到SN的路径
end
%% 初始化2
for i=1:n
    N(i).d=sqrt((N(i).x)^2+(N(i).y)^2);
    for j=i:n
        dis(i,j) = sqrt((N(i).x-N(j).x)^2+(N(i).y-N(j).y)^2);
        if dis(i,j)<=Rc
            N(i).nb = [N(i).nb,j];
            N(j).nb = [N(j).nb,i];
        end
        if 0<dis(i,j)&&dis(i,j)<=Rc/2
            N(i).nbhf = [N(i).nbhf,j];
            N(j).nbhf = [N(j).nbhf,i];
        end
    end
    if (N(i).x-ex)^2+(N(i).y-ey)^2<=Rc^2
        N(i).type=0;
        EN=[EN,i];
    end
    N(i).steps=ceil(N(i).d/Rc);
end
dis=dis+dis';
dis(dis>Rc|dis==0)=inf;
%% 寻找每一个节点的下一跳节点(倒序)
for ep=1:1
    figure(ep);
    hold on;
    para = [-sqrt(r), -sqrt(r), 2*sqrt(r), 2*sqrt(r)];
    rectangle('Position', para, 'Curvature', [1 1]);
    para = [ex-Rc, ey-Rc, 2*Rc, 2*Rc];
    rectangle('Position', para, 'Curvature', [1 1]);
    for i=1:5
        para = [-i*Rc, -i*Rc, 2*i*Rc, 2*i*Rc];
        rectangle('Position', para, 'Curvature', [1 1],'EdgeColor','r');
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
            scatter(N(i).x,N(i).y,'b');
        end
    end
    colorselect=0;
    color={'k','g','r','m','c','b'};
    colorlen=size(color,2);
    linetypeselect=0;
    linetype={'-','--','-.'};
    linetypelen=size(linetype,2);
    for i=EN
        N(i).path=[]; % 重置path
        curnode=i;
        while(curnode~=-1) % -1代表SN
            W=[];
            nbs=N(curnode).nb;
            if N(curnode).steps~=1
                for j=nbs
                    if N(j).E>0%&&N(j).type~=0
                        weight=a*N(j).steps+b*Eo/N(j).E;
                        W=[W,weight];
                    else
                        W=[W,inf];
                    end
                end
                [v,idx]=min(W);
                idxmins=find(W==v);
                nextstepnbmax=inf;
                for j=idxmins
                    nextstepnbtmp=0;
                    for k=N(nbs(j)).nb
                        if N(k).steps==N(i).steps-1
                            nextstepnbtmp=nextstepnbtmp+1;
                        end
                        if nextstepnbtmp>nextstepnbmax
                            idx=j;
                            nextstepnbmax=nextstepnbtmp;
                            v=W(idx);
                        end
                    end
                end    
                if v<inf
                    N(nbs(idx)).E=N(nbs(idx)).E-Er;
                    N(nbs(idx)).type=1;
                    N(curnode).E=N(curnode).E-Et;
                    N(i).path=[N(i).path,nbs(idx)];
                    colorstr=color{mod(colorselect,colorlen)+1}; % 选择颜色
                    linetypestr=linetype{mod(linetypeselect,linetypelen)+1}; % 选择线型
                    if N(nbs(idx)).steps==1
                        plot([N(curnode).x,N(nbs(idx)).x],[N(curnode).y,N(nbs(idx)).y],[linetypestr,colorstr]);
                        plot([N(nbs(idx)).x,0],[N(nbs(idx)).y,0],[linetypestr,colorstr]);
                        text(N(curnode).x,N(curnode).y,num2str(curnode));
                        text(N(nbs(idx)).x,N(nbs(idx)).y,num2str(nbs(idx)));
                        N(nbs(idx)).E=N(nbs(idx)).E-Et;
                        hold on;
                        curnode=-1;
                    else
                        plot([N(curnode).x,N(nbs(idx)).x],[N(curnode).y,N(nbs(idx)).y],[linetypestr,colorstr]);
                        text(N(curnode).x,N(curnode).y,num2str(curnode));
                        hold on;
                        curnode=nbs(idx);
                    end
                end
            else % 事件节点本身就是一跳，直接连接SN(这时的path保持[])
                colorstr=color{mod(colorselect,colorlen)+1}; % 选择颜色
                linetypestr=linetype{mod(linetypeselect,linetypelen)+1}; % 选择线型
                plot([N(curnode).x,0],[N(curnode).y,0],[linetypestr,colorstr]);
                hold on;
            end
        end
        colorselect=colorselect+1;
        linetypeselect=linetypeselect+1;
    end
end
%%



