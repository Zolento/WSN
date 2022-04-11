function [idxin,Win,INcds]=findbestIN(targetnodes,Win,INcds,N,curnode,Eo,i)
    for j=targetnodes
        flag=1; % 规定:不能复用同一条INpath上的IN；距离如果太远，也舍去
        if N(j).type==2
            for t=N(j).INn(:,1)
                if t==i % 该IN已经在同一条INpath中
                    flag=0;
                end
            end
%             if (N(j).x-N(curnode).x)^2+(N(j).y-N(curnode).y)^2>10^2
%                 flag=0;
%             end
        end
        if (N(j).type==-1||(N(j).type==2&&flag))&&N(j).E>0&&N(i).ANc~=1
            d=(N(j).x-N(curnode).x)^2+(N(j).y-N(curnode).y)^2;
            Win=[Win,0.3*N(j).credit+0.3*N(j).E/Eo+0.4*30/d];
            INcds=[INcds,j];
        end
    end
    [vin,idxin]=max(Win);
%     idxinmaxs=find(Win==vin); % 对符合要求的候选监督节点进一步筛选
%     W=inf;
%     if size(idxinmaxs,2)>1
%         for j=idxinmaxs % 如果信誉度一样，考虑离rp距离
%             d=(N(INcds(j)).x-N(curnode).x)^2+(N(INcds(j)).y-N(curnode).y)^2;
%             Wtmp=d;
%             % Wtmp=0.7*d/30+0.3*Win(j);
%             if Wtmp<W
%                 W=Wtmp;
%                 idxin=j;
%             end
%         end
%     end
end