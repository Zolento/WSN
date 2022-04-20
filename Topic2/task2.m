clc;
clear all;
close all;
%%
r=40000/pi;
rch=r/2;% 簇首成簇半径
para = [-sqrt(r), -sqrt(r), 2*sqrt(r), 2*sqrt(r)];
rectangle('Position', para, 'Curvature', [1 1]);
hold on;
xlim([-(sqrt(r)+10) sqrt(r)+10])
ylim([-(sqrt(r)+10) sqrt(r)+10])
axis equal
%%
N=load('N.mat');
N=N.N;
p=0.08; % 簇首占比
n=500;
Rc=30; % 通信半径
Eo=1; % 初始能量
dis=zeros(n); % 距离矩阵
EN=[];
inRc=[]; % 在SN的Rc半径范围内的节点集合
inRc_notEN=[]; % 在SN的Rc半径范围内但不是EN的节点集合
er=sqrt(r*rand(1));
etheta=2*pi*rand(1);
% ex=er*cos(etheta);
% ey=er*sin(etheta);
ex = 0;
ey = 0;
scatter(0,0,'r');% SN
hold on;
para = [-Rc, -Rc, 2*Rc, 2*Rc];
rectangle('Position', para, 'Curvature', [1 1],'EdgeColor','r');
hold on;
%% 初始化1 需要提前固定的字段
for i=1:n 
    N(i).type=-1;% 普通节点 -1 簇成员节点（MN 0） 簇首节点（CH 1） 监督节点（IN 2)
    N(i).EN=0; % 是否为事件节点
    N(i).E=Eo; % 初始能量
    N(i).credit=1; % 信誉度
    N(i).nb=[]; % Rc范围内邻居
    N(i).nbhf=[]; % Rc/2范围内邻居
    N(i).MN=[]; % 簇成员节点
    N(i).cluster=0; % 簇头
    N(i).inrange=0; % 能否与SN直接通信(即在SN的Rc半径范围内)
    N(i).d=0; % 到SN的距离
end
%% 初始化2
for i=1:n
    N(i).d=sqrt((N(i).x)^2+(N(i).y)^2);
    for j=i:n
        dis(i,j) = sqrt((N(i).x-N(j).x)^2+(N(i).y-N(j).y)^2);
        if Rc/2<dis(i,j)&&dis(i,j)<=Rc
            N(i).nb = [N(i).nb,j];
            N(j).nb = [N(j).nb,i];
        end
        if 0<dis(i,j)&&dis(i,j)<=Rc/2
            N(i).nbhf = [N(i).nbhf,j];
            N(j).nbhf = [N(j).nbhf,i];
        end
    end
%     if (N(i).x-ex)^2+(N(i).y-ey)^2<=50^2
    if (N(i).x-0)^2+(N(i).y-0)^2<=sqrt(r)^2
        if N(i).d<=50
            N(i).inrange=1;
            inRc=[inRc,i];
        end
        N(i).EN=1;
        N(i).rleft=0;
        EN=[EN,i];
        scatter(N(i).x,N(i).y,'g');
        hold on;
    else
        if N(i).d<=50
            N(i).inrange=1;
            inRc=[inRc,i];
            inRc_notEN=[inRc_notEN,i];
        end
        scatter(N(i).x,N(i).y,'b');
        hold on;
    end
end
dis=dis+dis';
dis(dis>Rc|dis==0)=inf;
%% 簇首选举，用这种方式确定锚点
found=zeros(1,size(EN,2));
target=EN;
chsall=[]; % 存放选举出的所有簇首
r=1;
while sum(found)<size(EN,2) % 对事件节点进行分簇
    t=(p/(1-p*(mod(r,1/p))));% 阈值
    tleft=mod(r,round(1/p));
    numCLheads=0;
    chs=[];
    target_tmp=[];
    chtarget=[];
    for t=target
        if found(EN==t)==0
            target_tmp=[target_tmp,t];
        end
    end
    for i=target_tmp % step1 选出候选簇头
        if N(i).rleft>0
           N(i).rleft=N(i).rleft-1;
        end
        if N(i).E>0 && N(i).rleft==0
            tmp=rand;
            if tmp<t*N(i).E/Eo
                N(i).rleft=round(1/p)-tleft;% 更新到下次选举轮数
                chtarget = [chtarget,i];% 本轮候选簇头集
                numCLheads = numCLheads+1;
            end
        end
    end
    if numCLheads>0
        conflict = zeros(numCLheads); 
        del = [];
        for i=1:numCLheads% step2 小于的Rc/2范围 簇头竞争
            for j=i:numCLheads
                if size(find(N(chtarget(i)).nbhf==chtarget(j)),2)~= 0
                    conflict(i,j)=1;
                end
            end
        end
        while sum(conflict(:))~=0
            [row,col]=find(conflict==1);
            if N(chtarget(row(1))).E>=N(chtarget(col(1))).E
                conflict(:,col(1))=0;
                del = [del,col(1)];
                N(chtarget(col(1))).type=-1;
            else
                conflict(row(1),:)=0;
                del = [del,row(1)];
                N(chtarget(row(1))).type=-1;
            end
        end
        for i=1:numCLheads % 划归簇成员节点
            if sum(i==del)==0
                N(chtarget(i)).type=1;
                chs = [chs,chtarget(i)];
                found(EN==chtarget(i))=1;
                scatter(N(chtarget(i)).x,N(chtarget(i)).y,'m');
                hold on;
                for j=N(chtarget(i)).nbhf
                    if N(j).type==-1&&N(j).EN==1% Rc/2内的普通节点一律声明为候选节点
                        N(j).type=0;
                        N(j).cluster=chtarget(i);
                        N(chtarget(i)).MN=[N(chtarget(i)).MN,j];
                        found(EN==j)=1;
                        plot([N(chtarget(i)).x,N(j).x],[N(chtarget(i)).y,N(j).y],'m');
                        hold on;
                    end
                end
            end
        end
        chsall=[chsall,chs];
    end
    r=r+1;
end
%% 均衡簇节点不够的簇首，为其多分配几个节点；选出每个簇的IN节点
% flag=0;
% while(~flag) 
%     flag=1;
%     for i=chsall
%         while(size(N(i).MN,2)<1)
%             nbhf_sorted=N(i).nbhf;
%             length=size(N(i).nbhf,2);
%             for p=1:length-1
%                 for q=1:length-p
%                     if dis(i,nbhf_sorted(q))>dis(i,nbhf_sorted(q+1))
%                         tmp=nbhf_sorted(q);
%                         nbhf_sorted(q)=nbhf_sorted(q+1);
%                         nbhf_sorted(q+1)=tmp;
%                     end
%                 end
%             end
%             
%             for j=nbhf_sorted
%                 if N(j).EN==1&&N(j).type==0&&size(N(N(j).cluster).MN,2)>3&&size(N(i).MN,2)<3
%                     N(i).MN=[N(i).MN,j];
%                     tmp=[];
%                     for k=N(N(j).cluster).MN
%                         if k~=j
%                             tmp=[tmp,k];
%                         end
%                     end
%                     N(N(j).cluster).MN=tmp;
%                     N(j).cluster=i;
%                     plot([N(i).x,N(j).x],[N(i).y,N(j).y],'c');
%                     hold on;
%                     flag=0; % 发生了重分配，flag置0
%                 end
%             end
%         end
%     end
% end
%% 计算簇首邻接矩阵
% chsalllen=size(chsall,2);
% dis_2=zeros(chsalllen+1);
% mind_2=inf;
% numinRc=0;
% innerRP=[];
% T=2;
% for i=1:chsalllen
%     if mind_2(1)>sqrt(N(chsall(i)).x^2+N(chsall(i)).y^2)
%        mind_2=[sqrt(N(chsall(i)).x^2+N(chsall(i)).y^2),i];
%     end
%     if N(chsall(i)).inrange==1
%         dis_2(i,chsalllen+1)=N(chsall(i)).d^2;
%         numinRc=numinRc+1;
%     end
%     for j=i:chsalllen
%         if size(find(N(chsall(i)).nb==chsall(j)),2)~=0
%             dis_2(i,j)=dis(chsall(i),chsall(j))^2;
% %             plot([N(chsall(i)).x,N(chsall(j)).x],[N(chsall(i)).y,N(chsall(j)).y],'k');
% %             hold on;
%         end
%     end
% end
% dis_2=dis_2+dis_2';
% dis_2(dis_2==0|dis_2>Rc^2)=inf;
% if 0<=numinRc&&numinRc<=T % 如果在一跳范围内的簇首个数小于一定的值
%  % 把最近的簇首的内圈范围内的点都用作中继(实际上内圈最远的情况也差不多略大于Rc)
%     root=mind_2(2);
%     flag=0;
%     for i=inRc_notEN
%         if dis(i,chsall(root))<N(chsall(root)).d&&N(i).d<N(chsall(root)).d
%             innerRP=[innerRP,i];
%         end
%     end
%     targetnode=[chsall,innerRP]; % 添加内圈节点
%     targetnodelen=size(targetnode,2);
%     dis_2=zeros(targetnodelen+1); % 重新计算距离矩阵
%     for i=1:targetnodelen
%         if N(targetnode(i)).inrange==1
%             dis_2(i,targetnodelen+1)=N(targetnode(i)).d^2;
%             numinRc=numinRc+1;
%         end
%         for j=i:targetnodelen
%             dis_2(i,j)=dis(targetnode(i),targetnode(j))^2;
%         end
%     end
%     dis_2=dis_2+dis_2';
%     dis_2(dis_2==0|dis_2>Rc^2)=inf;
%     found=zeros(1,targetnodelen+1);
%     found(targetnodelen+1)=1;
%     d_2=zeros(1,targetnodelen+1);
%     d_2(1:targetnodelen)=inf;
%     pre_node=-1*ones(1,targetnodelen+1);
%     path=-1*ones(1,targetnodelen+1);
%     root=targetnodelen+1; % root重新指向SN
% else % 如果在一跳范围内的簇首个数足够多
%     flag=1;
%     found=zeros(1,chsalllen+1);
%     found(chsalllen+1)=1;
%     d_2=zeros(1,chsalllen+1);
%     d_2(1:chsalllen)=inf;
%     pre_node=-1*ones(1,chsalllen+1);
%     path=-1*ones(1,chsalllen+1);
%     root=chsalllen+1; % root直接指向SN
% end
% if ~flag
%     pathnodes=targetnode;
% else
%     pathnodes=chsall;
% end
% % 画出最近的簇首的内圈范围
% SNsr=sqrt(N(chsall(mind_2(2))).x^2+N(chsall(mind_2(2))).y^2);
% para = [-SNsr, -SNsr, 2*SNsr, 2*SNsr];
% rectangle('Position', para, 'Curvature', [1 1]);
% hold on;
% % 画出参与路由的节点序号
% for i=pathnodes
%     c = num2str(i);
%     text(N(i).x,N(i).y,c);
%     hold on;
% end
% %% 寻找簇首到SN的路径(d^2最小)
% SNs=root;
% % 求路径 run dijkstra
% while sum(found)<size(found,2)  %看是否所有的点都标记为P标号
%     target=find(found==0); 
%     dtmp=d_2;
%     d_2(target)=min(d_2(target),d_2(SNs)+dis_2(SNs,target));%对还没有存入found的点更新最短通路
%     minp_nodeidx=find(d_2(target)==min(d_2(target)));
%     updated=find(d_2(target)~=dtmp(target));%最短路径被更新过的节点
%     if size(updated,2)~=0
%         pre_node(target(updated))=SNs;
%     end
% %     for i=minp_nodeidx
%     SNs=target(minp_nodeidx(1));   %将当前符合最短通路的且离SN最近的节点设为源点
%     found(SNs)=1;%已找到最小路径 found(i)=1
%     path(SNs)=pre_node(SNs);
%     k=path(SNs); % 画出所有路径
%     p=SNs;
%     if k<=size(pathnodes,2)
%         plot([N(pathnodes(p)).x,N(pathnodes(k)).x],[N(pathnodes(p)).y,N(pathnodes(k)).y],'k');
%         hold on;
%     else
%         plot([N(pathnodes(p)).x,0],[N(pathnodes(p)).y,0],'k');
%         hold on;
%     end
% end
% path(path==size(path,2))=-1; % 将指向SN的置为-1
% %% 输出路径 v2;找前一跳是RP的簇首
% path_traceback=cell(1,chsalllen);
% E=zeros(1,chsalllen);
% Ed=zeros(1,chsalllen);
% N_RP=[];
% for i=1:chsalllen
%     curnode=i;
%     Ed(i)=N(pathnodes(i)).d^2;
%     while curnode~=-1
%         path_traceback{i}=[pathnodes(curnode),path_traceback{i}];
%         curnode=path(curnode);
%     end
%     pathtmp=path_traceback{i};
%     if size(pathtmp,2)>2&&size(find(innerRP==pathtmp(end-1)),2)~=0
%         N_RP=[N_RP,i];
%     end
%     for j=1:size(pathtmp,2)
%         if j==1
%             E(i)=E(i)+N(pathtmp(j)).x^2+N(pathtmp(j)).y^2;
%             plot([0,N(pathtmp(j)).x],[0,N(pathtmp(j)).y],'k');
%             hold on;
%         else
%             E(i)=E(i)+(N(pathtmp(j)).x-N(pathtmp(j-1)).x)^2+(N(pathtmp(j)).y-N(pathtmp(j-1)).y)^2;
%             plot([N(pathtmp(j)).x,N(pathtmp(j-1)).x],[N(pathtmp(j)).y,N(pathtmp(j-1)).y],'k');
%             hold on;
%         end
%     end
% end
%%






