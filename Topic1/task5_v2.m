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
a(a==0|a>225)=inf;% ����ͳ���������ֵ�Ľڵ���Ϊinf
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

rdata = zeros(1,n);% ÿ�������յ�����Ϣ����
sent = zeros(1,n);% �����ھӷ��͹����ݵı�ʶ
prenode = zeros(1,n);% ǰ��ڵ㣬��ֹ�ش�(ͬ���ڵ�����ɴ��ڻش�)
allflag = 0;% ֹͣ��־�����нڵ㶼�յ�����Ϣ����1
cnt = 0;% ��¼����
nbcom = 0;% ͨ�ŵ��ھ���
SN = 3;
ST = SN;% һ��ʼ������Ϣ�Ľڵ�


while ~allflag
    pause(0.1);
    if ~sent(ST)
        %rdata(nb{ST}) = rdata(nb{ST})+1;
        alltarget = ST;
        sent(ST) = 1;
    else
        alltarget = find(rdata~=0&sent~=1);
        sent(alltarget) = 1;% ��־Ϊ�ѷ���
    end
    for i = 1:size(alltarget,2)
        spread_target = nb{alltarget(i)};% ��ɢĿ��;�ҳ��ھ�
        spread_target = spread_target(spread_target~=prenode(alltarget(i)));% ����ɢĿ����ȥ����ǰĿ���ǰ��ڵ�
        spread_target = spread_target(rdata(spread_target)==0);% ����ɢĿ����ȥ�����յ���Ϣ�Ľڵ�
        prenode(spread_target(rdata(spread_target)==0)) = alltarget(i);% ������ɢĿ���ǰ��ڵ�
        rdata(spread_target)=rdata(spread_target)+1;% ��������
        nbcom = nbcom+size(spread_target,2);% ���㷢����
        for j = 1:size(spread_target,2)
            plot([x(alltarget(i)),x(spread_target(j))],[y(alltarget(i)),y(spread_target(j))],'m');%����target�㵽�����ھӵ�·��
            hold on;
        end
    end
    cnt = cnt+1;
    allflag = size(find(rdata==0),2)==1;% ���һ���յ��Ľڵ㲻��ش����ھ���
end
avgdata = sum(rdata)/(n-1);
avgnbcom = nbcom/n;