clc;
clear all;
close all;

r=40000/pi;
para = [-sqrt(r), -sqrt(r), 2*sqrt(r), 2*sqrt(r)];
rectangle('Position', para, 'Curvature', [1 1]);
hold on;
xlim([-(sqrt(r)+10) sqrt(r)+10])
ylim([-(sqrt(r)+10) sqrt(r)+10])
axis equal

N=struct;
er=sqrt(r*rand(1));
etheta=2*pi*rand(1);
ex=er*cos(etheta);
ey=er*sin(etheta);
scatter(0,0,'r');
hold on;

for i=1:500
    N(i).r=sqrt(r*rand(1));
    N(i).theta=2*pi*rand(1);
    N(i).x=N(i).r*cos(N(i).theta);
    N(i).y=N(i).r*sin(N(i).theta);
    scatter(N(i).x,N(i).y,'b');
    hold on;
end




