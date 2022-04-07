clc;
clear all;
close all;

data = load("lab1data.txt");
x = data(:,1);
y = data(:,2);
scatter(x,y);
hold on;

n = size(data);
n = n(1);
min = x(1)^2+y(1)^2;
SN = 1;
for i=2:n
    if min>x(i)^2+y(i)^2
        min = x(i)^2+y(i)^2;
        SN = i;
    end
end

scatter(x(SN),y(SN), 'r');
    