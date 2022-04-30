function [AUC_results,idx]=runiforest(ADLabels,Data,ADprop)
% ADLabels: groundtruth labels, PosLabel = 1& NegLabel =0;
% 正样本：异常节点；负样本：正常节点
rounds = 10; % rounds of repeat
% parameters for iForest
NumTree = 100; % number of isolation trees
NumSub = size(Data, 1); % subsample size
NumDim = size(Data, 2); % do not perform dimension sampling 
auc = zeros(rounds, 1);
mtime = zeros(rounds, 2);
rseed = zeros(rounds, 1);

for r = 1:rounds
    % disp(['rounds ', num2str(r), ':']);
    
    rseed(r) = sum(100 * clock);
    Forest = IsolationForest(Data, NumTree, NumSub, NumDim, rseed(r));
    mtime(r, 1) = Forest.ElapseTime;
%    [Mass, mtime(r, 2)] = IsolationEstimation(Data, Forest);
    [Mass, ~] = IsolationEstimation(Data, Forest);
    Score = - mean(Mass, 2);
    [auc(r), idx] = Measure_AUC(Score, ADLabels,ADprop);
    % disp(['auc = ', num2str(auc(r)), '.']);
%    [~,~,~,AUClog(r)] = perfcurve(logical(ADLabels),Score,'true');    
end

% myresults = [mean(auc), var(auc), mean(mtime(:, 1)), mean(mtime(:, 2))] 
AUC_results = [mean(auc), std(auc)]; % average AUC over 10 trials
end