function [idxin,Win,INcds]=findbestIN(targetnodes,Win,INcds,N,curnode)
for j=targetnodes
    if N(j).type==-1||N(j).type==2
        Win=[Win,N(j).credit];
        INcds=[INcds,j];
    end
end
[vin,idxin]=max(Win);
idxinmaxs=find(Win==vin); % 对符合要求的候选监督节点进一步筛选
% nextstepinmin=inf;
if size(idxinmaxs,2)>1
%     for j=idxinmaxs % 找出跳数最小的的候选监督节点
%         nextstepintmp=N(INcds(j)).steps;
%     %     for k=N(targetnodes(j)).nb
%     %         if N(k).steps==N(rp).steps-1
%     %             nextstepintmp=nextstepintmp+1;
%     %         end
%         if nextstepintmp<nextstepinmin
%             idxin=j;
%             nextstepinmin=nextstepintmp;
%         end
%     end
%     if nextstepinmin>=N(curnode).steps % 如果选出的节点跳数没有减小，可能周围的可用路径已经很少了
    nextstepinmin=inf;
    for j=idxinmaxs % 找出离当前rp最近的候选监督节点
        nextstepintmp=(N(INcds(j)).x-N(curnode).x)^2+(N(INcds(j)).y-N(curnode).y)^2;
        if nextstepintmp<nextstepinmin
            idxin=j;
            nextstepinmin=nextstepintmp;
        end
    end
%     end
end
end