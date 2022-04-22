clear all;
close all;
load 'rppsummary.mat' rppsummary
load 'NwithResults.mat' N
%% 记录mexep轮后的rpp
% 就是rppsummary的最后一列
n=500;
maxep=400;
sr=200;
lastrpp=zeros(1,n);
for ep=1:maxep
    if mod(ep,sr)==0
        rd=floor(ep/sr);
        % 就是rppsummary的最后一列
        rppsummaryslice=rppsummary(:,(rd-1)*sr+1:rd*sr);
        varsummary=zeros(1,n);
        X=zeros(n,3);
        for i=1:n
            lastrpp(i)=N(i).credit;
        end
%         figure(maxep+4*(rd-1)+1)
%         title('记录的累计转发率分布')
%         for i=1:n
%             if ~N(i).del
%                 if N(i).AN
%                     scatter(i,lastrpp(i),'r')
%                     hold on
%                 else
%                     scatter(i,lastrpp(i),'b')
%                     hold on
%                 end
%                 if N(i).isabnormal
%                     scatter(i,lastrpp(i),'k','+')
%                     hold on
%                 end
%             end
%         end
%         figure(maxep+4*(rd-1)+2)
%         title('理想的转发率平均值分布')
%         for i=1:n
%             if ~N(i).del
%                 if N(i).AN
%                     scatter(i,N(i).rppavg,'r')
%                     hold on
%                 else
%                     scatter(i,N(i).rppavg,'b')
%                     hold on
%                 end
%                 if N(i).isabnormal
%                     scatter(i,N(i).rppavg,'k','+')
%                     hold on
%                 end
%             end
%         end
        % dbscan聚类
        for i=1:n
            tmp=[];
            if ~N(i).del
                pretmp=rppsummaryslice(i,1);
                for j=1:sr
                    if rppsummaryslice(i,j)<1&&rppsummaryslice(i,j)~=pretmp
                        tmp=[tmp,rppsummaryslice(i,j)];
                        pretmp=rppsummaryslice(i,j);
                    end
                end
                if N(i).ispath
                    varsummary(i)=var([1,tmp]);
                else
                    varsummary(i)=0;
                end
                X(i,1)=lastrpp(i);
                X(i,2)=varsummary(i);
                X(i,3)=var([1,tmp])/mean([1,tmp]);
            else
                X(i,:)=0;
            end
        end
        X(:,2)=X(:,2)/max(X(:,2));
        X(:,3)=X(:,3)/max(X(:,3));
        figure(maxep+4*(rd-1)+3)
        title('聚类结果')
%         xlabel('normed VAR')
%         ylabel('CFR')
        idx = dbscan(X,0.1,5);
        for i=1:n
            if ~N(i).del
                if idx(i)==-1
                    scatter3(X(i,1),X(i,2),X(i,3),'r');
                    hold on;
                elseif idx(i)==2
                    scatter3(X(i,1),X(i,2),X(i,3),'b');
                    hold on;
                else
                    scatter3(X(i,1),X(i,2),X(i,3),'g');
                    hold on;
                end
            end
        end
        figure(maxep+4*(rd-1)+4)
        title('实际情况')
%         xlabel('normed VAR')
%         ylabel('CFR')
        for i=1:n
            if ~N(i).del
                if N(i).AN
                    scatter3(X(i,1),X(i,2),X(i,3),'r')
                    hold on
                else
                    scatter3(X(i,1),X(i,2),X(i,3),'b')
                    hold on
                end
                if N(i).isabnormal
                    scatter3(X(i,1),X(i,2),X(i,3),'k','+');
                    hold on;
                end
                if ~N(i).ispath
                    scatter3(X(i,1),X(i,2),X(i,3),'k','*');
                    hold on;
                end
            end
        end
        % 计算误检率
        numANc=0;
        numFANc=0;
        for i=1:n
            N(i).ANc=idx(i);
            if N(i).ANc==-1&&N(i).ispath&&~N(i).del
                N(i).del=1;
                numANc=numANc+1;
                if ~N(i).AN
                    numFANc=numFANc+1;
                end
            end
        end
        error=numFANc/numANc;
        fprintf('第%d个200轮检测正确率:%f\n',rd,1-error);
        % fprintf('第%d个200轮漏检率:%f',rd,1-error);
    end
end