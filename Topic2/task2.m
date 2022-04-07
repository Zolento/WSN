clc;
clear all;
close all;
%%
r=40000/pi;
rch=r/2;% ���׳ɴذ뾶
para = [-sqrt(r), -sqrt(r), 2*sqrt(r), 2*sqrt(r)];
rectangle('Position', para, 'Curvature', [1 1]);
hold on;
xlim([-(sqrt(r)+10) sqrt(r)+10])
ylim([-(sqrt(r)+10) sqrt(r)+10])
axis equal
%%
N=load('N.mat');
N=N.N;
p=0.08; % ����ռ��
n=500;
Rc=50; % ͨ�Ű뾶
Eo=1; % ��ʼ����
dis=zeros(n); % �������
EN=[];
inRc=[]; % ��SN��Rc�뾶��Χ�ڵĽڵ㼯��
inRc_notEN=[]; % ��SN��Rc�뾶��Χ�ڵ�����EN�Ľڵ㼯��
er=sqrt(r*rand(1));
etheta=2*pi*rand(1);
ex=er*cos(etheta);
ey=er*sin(etheta);
scatter(0,0,'r');% SN
hold on;
para = [-Rc, -Rc, 2*Rc, 2*Rc];
rectangle('Position', para, 'Curvature', [1 1],'EdgeColor','r');
hold on;
%% ��ʼ��1 ��Ҫ��ǰ�̶����ֶ�
for i=1:n 
    N(i).type=-1;% ��ͨ�ڵ� -1 �س�Ա�ڵ㣨MN 0�� ���׽ڵ㣨CH 1�� �ල�ڵ㣨IN 2)
    N(i).EN=0; % �Ƿ�Ϊ�¼��ڵ�
    N(i).E=Eo; % ��ʼ����
    N(i).credit=1; % ������
    N(i).nb=[]; % Rc��Χ���ھ�
    N(i).nbhf=[]; % Rc/2��Χ���ھ�
    N(i).MN=[]; % �س�Ա�ڵ�
    N(i).cluster=0; % ��ͷ
    N(i).inrange=0; % �ܷ���SNֱ��ͨ��(����SN��Rc�뾶��Χ��)
    N(i).d=0; % ��SN�ľ���
end
%% ��ʼ��2
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
    if (N(i).x-ex)^2+(N(i).y-ey)^2<=50^2
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
%% ����ѡ��
found=zeros(1,size(EN,2));
target=EN;
chsall=[]; % ���ѡ�ٳ������д���
r=1;
while sum(found)<size(EN,2) % ���¼��ڵ���зִ�
    t=(p/(1-p*(mod(r,1/p))));% ��ֵ
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
    for i=target_tmp % step1 ѡ����ѡ��ͷ
        if N(i).rleft>0
           N(i).rleft=N(i).rleft-1;
        end
        if N(i).E>0 && N(i).rleft==0
            tmp=rand;
            if tmp<t*N(i).E/Eo
                N(i).rleft=round(1/p)-tleft;% ���µ��´�ѡ������
                chtarget = [chtarget,i];% ���ֺ�ѡ��ͷ��
                numCLheads = numCLheads+1;
            end
        end
    end
    if numCLheads>0
        conflict = zeros(numCLheads); 
        del = [];
        for i=1:numCLheads% step2 С�ڵ�Rc/2��Χ ��ͷ����
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
        for i=1:numCLheads % ����س�Ա�ڵ�
            if sum(i==del)==0
                N(chtarget(i)).type=1;
                chs = [chs,chtarget(i)];
                found(EN==chtarget(i))=1;
                scatter(N(chtarget(i)).x,N(chtarget(i)).y,'m');
                hold on;
                for j=N(chtarget(i)).nbhf
                    if N(j).type==-1&&N(j).EN==1% Rc/2�ڵ���ͨ�ڵ�һ������Ϊ��ѡ�ڵ�
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
%% ����ؽڵ㲻���Ĵ��ף�Ϊ�����伸���ڵ㣻ѡ��ÿ���ص�IN�ڵ�
flag=0;
while(~flag) 
    flag=1;
    for i=chsall
        while(size(N(i).MN,2)<3)
            nbhf_sorted=N(i).nbhf;
            length=size(N(i).nbhf,2);
            for p=1:length-1
                for q=1:length-p
                    if dis(i,nbhf_sorted(q))>dis(i,nbhf_sorted(q+1))
                        tmp=nbhf_sorted(q);
                        nbhf_sorted(q)=nbhf_sorted(q+1);
                        nbhf_sorted(q+1)=tmp;
                    end
                end
            end
            
            for j=nbhf_sorted
                if N(j).EN==1&&N(j).type==0&&size(N(N(j).cluster).MN,2)>3&&size(N(i).MN,2)<3
                    N(i).MN=[N(i).MN,j];
                    tmp=[];
                    for k=N(N(j).cluster).MN
                        if k~=j
                            tmp=[tmp,k];
                        end
                    end
                    N(N(j).cluster).MN=tmp;
                    N(j).cluster=i;
                    plot([N(i).x,N(j).x],[N(i).y,N(j).y],'c');
                    hold on;
                    flag=0; % �������ط��䣬flag��0
                end
            end
        end
    end
end
for i=chsall
    
end
%% ��������ڽӾ���
chsalllen=size(chsall,2);
dis_2=zeros(chsalllen+1);
mind_2=inf;
numinRc=0;
innerRP=[];
T=2;
for i=1:chsalllen
    if mind_2(1)>sqrt(N(chsall(i)).x^2+N(chsall(i)).y^2)
       mind_2=[sqrt(N(chsall(i)).x^2+N(chsall(i)).y^2),i];
    end
    if N(chsall(i)).inrange==1
        dis_2(i,chsalllen+1)=N(chsall(i)).d^2;
        numinRc=numinRc+1;
    end
    for j=i:chsalllen
        if size(find(N(chsall(i)).nb==chsall(j)),2)~=0
            dis_2(i,j)=dis(chsall(i),chsall(j))^2;
%             plot([N(chsall(i)).x,N(chsall(j)).x],[N(chsall(i)).y,N(chsall(j)).y],'k');
%             hold on;
        end
    end
end
dis_2=dis_2+dis_2';
dis_2(dis_2==0|dis_2>Rc^2)=inf;
if 0<=numinRc&&numinRc<=T % �����һ����Χ�ڵĴ��׸���С��һ����ֵ
 % ������Ĵ��׵���Ȧ��Χ�ڵĵ㶼�����м�(ʵ������Ȧ��Զ�����Ҳ����Դ���Rc)
    root=mind_2(2);
    flag=0;
    for i=inRc_notEN
        if dis(i,chsall(root))<N(chsall(root)).d&&N(i).d<N(chsall(root)).d
            innerRP=[innerRP,i];
        end
    end
    targetnode=[chsall,innerRP]; % �����Ȧ�ڵ�
    targetnodelen=size(targetnode,2);
    dis_2=zeros(targetnodelen+1); % ���¼���������
    for i=1:targetnodelen
        if N(targetnode(i)).inrange==1
            dis_2(i,targetnodelen+1)=N(targetnode(i)).d^2;
            numinRc=numinRc+1;
        end
        for j=i:targetnodelen
            dis_2(i,j)=dis(targetnode(i),targetnode(j))^2;
        end
    end
    dis_2=dis_2+dis_2';
    dis_2(dis_2==0|dis_2>Rc^2)=inf;
    found=zeros(1,targetnodelen+1);
    found(targetnodelen+1)=1;
    d_2=zeros(1,targetnodelen+1);
    d_2(1:targetnodelen)=inf;
    pre_node=-1*ones(1,targetnodelen+1);
    path=-1*ones(1,targetnodelen+1);
    root=targetnodelen+1; % root����ָ��SN
else % �����һ����Χ�ڵĴ��׸����㹻��
    flag=1;
    found=zeros(1,chsalllen+1);
    found(chsalllen+1)=1;
    d_2=zeros(1,chsalllen+1);
    d_2(1:chsalllen)=inf;
    pre_node=-1*ones(1,chsalllen+1);
    path=-1*ones(1,chsalllen+1);
    root=chsalllen+1; % rootֱ��ָ��SN
end
if ~flag
    pathnodes=targetnode;
else
    pathnodes=chsall;
end
% ��������Ĵ��׵���Ȧ��Χ
SNsr=sqrt(N(chsall(mind_2(2))).x^2+N(chsall(mind_2(2))).y^2);
para = [-SNsr, -SNsr, 2*SNsr, 2*SNsr];
rectangle('Position', para, 'Curvature', [1 1]);
hold on;
% ��������·�ɵĽڵ����
for i=pathnodes
    c = num2str(i);
    text(N(i).x,N(i).y,c);
    hold on;
end
%% Ѱ�Ҵ��׵�SN��·��(d^2��С)
SNs=root;
% ��·�� run dijkstra
while sum(found)<size(found,2)  %���Ƿ����еĵ㶼���ΪP���
    target=find(found==0); 
    dtmp=d_2;
    d_2(target)=min(d_2(target),d_2(SNs)+dis_2(SNs,target));%�Ի�û�д���found�ĵ�������ͨ·
    minp_nodeidx=find(d_2(target)==min(d_2(target)));
    updated=find(d_2(target)~=dtmp(target));%���·�������¹��Ľڵ�
    if size(updated,2)~=0
        pre_node(target(updated))=SNs;
    end
%     for i=minp_nodeidx
    SNs=target(minp_nodeidx(1));   %����ǰ�������ͨ·������SN����Ľڵ���ΪԴ��
    found(SNs)=1;%���ҵ���С·�� found(i)=1
    path(SNs)=pre_node(SNs);
    k=path(SNs); % ��������·��
    p=SNs;
    if k<=size(pathnodes,2)
        plot([N(pathnodes(p)).x,N(pathnodes(k)).x],[N(pathnodes(p)).y,N(pathnodes(k)).y],'k');
        hold on;
    else
        plot([N(pathnodes(p)).x,0],[N(pathnodes(p)).y,0],'k');
        hold on;
    end
end
path(path==size(path,2))=-1; % ��ָ��SN����Ϊ-1
%% ���·�� v2;��ǰһ����RP�Ĵ���
path_traceback=cell(1,chsalllen);
E=zeros(1,chsalllen);
Ed=zeros(1,chsalllen);
N_RP=[];
for i=1:chsalllen
    curnode=i;
    Ed(i)=N(pathnodes(i)).d^2;
    while curnode~=-1
        path_traceback{i}=[pathnodes(curnode),path_traceback{i}];
        curnode=path(curnode);
    end
    pathtmp=path_traceback{i};
    if size(pathtmp,2)>2&&size(find(innerRP==pathtmp(end-1)),2)~=0
        N_RP=[N_RP,i];
    end
    for j=1:size(pathtmp,2)
        if j==1
            E(i)=E(i)+N(pathtmp(j)).x^2+N(pathtmp(j)).y^2;
            plot([0,N(pathtmp(j)).x],[0,N(pathtmp(j)).y],'k');
            hold on;
        else
            E(i)=E(i)+(N(pathtmp(j)).x-N(pathtmp(j-1)).x)^2+(N(pathtmp(j)).y-N(pathtmp(j-1)).y)^2;
            plot([N(pathtmp(j)).x,N(pathtmp(j-1)).x],[N(pathtmp(j)).y,N(pathtmp(j-1)).y],'k');
            hold on;
        end
    end
end
%% Ѱ��RP�ļල�ڵ�
for i=N_RP
    for j=path_traceback{i}
        
    end
end






