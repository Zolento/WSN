clc;
clear all;
close all;
%%
r=40000/pi;
rch=r/2;% ���׳ɴذ뾶
%%
N=load('N.mat');
N=N.N;
a=0.3;
b=1-a;
p=0.08; % ����ռ��
n=500;
Rc=30; % ͨ�Ű뾶
Rs=10; % ��Ӧ�뾶
Eo=1; % ��ʼ����
Et=0.0003; % ��������
Er=0.0001; % ��������
dis=zeros(n); % �������
EN=[];
onestep=[];
er=sqrt(r*rand(1));
etheta=2*pi*rand(1);
ex=er*cos(etheta);
ey=er*sin(etheta);
%% ��ʼ��1 ��Ҫ��ǰ�̶����ֶ�
for i=1:n 
    N(i).type=-1; % ��ͨ�ڵ� -1 �¼��ڵ�(EN 0) �м̽ڵ㣨RP 1�� �ල�ڵ㣨IN 2)
    N(i).steps=0; % ����
    N(i).E=Eo; % ��ʼ����
    N(i).IN=[]; % �ýڵ������ļල�ڵ�(������һ������)
    N(i).INn=0; % ����Ǽල�ڵ㣬�ýڵ����ල�Ľڵ� 
    N(i).credit=1; % ������
    N(i).nb=[]; % Rc��Χ���ھ�
    N(i).nbhf=[]; % Rc/2��Χ���ھ�
    N(i).d=0; % ��SN�ľ���
    N(i).path=[]; % ��SN��·��
    N(i).INpath=[]; % IN·��
end
%% ��ʼ��2
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
N(n+1).type=2;
N(n+1).steps=0;
dis=dis+dis';
dis(dis>Rc|dis==0)=inf;
%% Ѱ��ÿһ���ڵ����һ���ڵ�(����)
for ep=1:1
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
            N(i).type=-1;
            scatter(N(i).x,N(i).y,'b');
        end
        
    end
    colorselect=0;
    color={'k','r','c','b'};
    colorlen=size(color,2);
    linetypeselect=0;
    linetype={'-','--','-.'};
    linetypelen=size(linetype,2);
    for i=EN
        N(i).path=[]; % ����path
        
        curnode=i;
        while(curnode~=-1) % -1����SN
            nbs=N(curnode).nb;
            Wrp=[];
%             INs=[];
            foundIN=0; % �ж���û���ҵ����е�IN
            if N(curnode).steps~=1
                for j=nbs %1.�ҳ�·�ɽڵ��·��
                    if N(j).E>0&&N(j).type~=2 % �ල�ڵ㲻��Ϊ·��,�¼��ڵ������·��
                        weight=a*N(j).steps+b*Eo/N(j).E;
                        Wrp=[Wrp,weight];
                    else
                        Wrp=[Wrp,inf];
                    end
%                     if N(j).type==2
%                         foundIN=1;
%                         INs=[INs,j];
%                     end
                end
                [v,idx]=min(Wrp);
                idxmins=find(Wrp==v); % �Է�����СȨֵ��·�ɽڵ��һ��ɸѡ
                nextstepnbmax=inf;
                if size(idxmins,2)>1
                    for j=idxmins % �ҳ���һ���ھ�������СȨֵ·�ɽڵ�
                        nextstepnbtmp=0;
                        for k=N(nbs(j)).nb
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
                N(idx).type=1;  
                if v<inf
                    N(nbs(idx)).E=N(nbs(idx)).E-Er;
                    N(nbs(idx)).type=1;
                    N(curnode).E=N(curnode).E-Et;
                    N(i).path=[N(i).path,nbs(idx)];
                    colorstr=color{mod(colorselect,colorlen)+1}; % ѡ����ɫ
                    linetypestr=linetype{mod(linetypeselect,linetypelen)+1}; % ѡ������
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
            else % �¼��ڵ㱾�����һ����ֱ������SN(��ʱ��path���Ǳ���[])
                colorstr=color{mod(colorselect,colorlen)+1}; % ѡ����ɫ
                linetypestr=linetype{mod(linetypeselect,linetypelen)+1}; % ѡ������
                plot([N(curnode).x,0],[N(curnode).y,0],[linetypestr,colorstr]);
                hold on;
                curnode=-1;
            end
        end
        N(i).path=[N(i).path,n+1];
        preIN=i;
        preRP=i;
        allpath=[i,N(i).path];
        for rp=allpath
            rpnbs=N(rp).nb;
            Win=[];
            INcds=[];
            if rp==i % ��ʼ
                targetnodes=intersect(rpnbs,N(allpath(2)).nb);
            elseif rp==n+1 % ����SN
                if sum(find(rpnbs==preIN))~=0 % �����һ��IN�Ѿ���һ���ڵ���
                    N(i).INpath=[N(i).INpath,n+1];
                    plot([N(preIN).x,0],[N(preIN).y,0],[linetypestr,'m']);
                    hold on;
                    break;
                else
                    targetnodes=intersect(rpnbs,N(preIN).nb);
                end
            else
                targetnodes=intersect(rpnbs,N(preIN).nb);
            end
%             if sum(targetnodes)==0
%                 
%             end
            [idxin,Win,INcds]=findbestIN(targetnodes,Win,INcds,N,rp);
            plot([N(INcds(idxin)).x,N(preIN).x],[N(INcds(idxin)).y,N(preIN).y],[linetypestr,'m']);
            plot([N(INcds(idxin)).x,N(rp).x],[N(INcds(idxin)).y,N(rp).y],['-','g']);
            text(N(INcds(idxin)).x,N(INcds(idxin)).y,num2str(INcds(idxin)));
            hold on;
            N(i).INpath=[N(i).INpath,INcds(idxin)];
            N(INcds(idxin)).type=2;
            N(rp).IN=INcds(idxin);
            preIN=N(rp).IN;
            preRP=rp;
        end
        colorselect=colorselect+1;
        linetypeselect=linetypeselect+1;
    end
end
%% 



