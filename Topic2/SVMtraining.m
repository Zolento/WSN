clear;
close all;
load 'Data.mat' Data
tmp1=Data;
load 'Data1.mat' Data
Data=[Data;tmp1];
load 'ADLabels.mat' ADLabels
tmp1=ADLabels;
load 'ADLabels1.mat' ADLabels
ADLabels=[ADLabels;tmp1];
figure(1);
newData=[];
newADLabels=[];
for i=1:size(Data,1)
    if ~ADLabels(i)
        newData=[newData;Data(i,:)];
        newADLabels=[newADLabels;ADLabels(i)];
    else
        newData=[newData;Data(i,:)];
        newADLabels=[newADLabels;ADLabels(i)];
    end
end
for i=1:size(newData,1)
    if newADLabels(i)
        text(newData(i,1),newData(i,2),num2str(i));
        h1=scatter(newData(i,1),newData(i,2),'r');
        hold on;
    else
        h2=scatter(newData(i,1),newData(i,2),'b');
        hold on;
    end
end
xlabel('CFR');
ylabel('σ/μ');
legend([h1,h2],'恶意节点','正常节点')
title('实际情况');
SVMModel=fitcsvm(newData,newADLabels,'KernelFunction','gaussian','KernelScale','auto');
classOrder = SVMModel.ClassNames;
%%
figure(2);
label = predict(SVMModel,newData);
for i=1:size(newData,1)
    if label(i)
        h1=scatter(newData(i,1),newData(i,2),'r');
        hold on;
    else
        h2=scatter(newData(i,1),newData(i,2),'b');
        hold on;
    end
end