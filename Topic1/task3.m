clc;
clear all;
close all;

C = ones(1,3);
B = 500*1e3;
pn = 1.38e-23*300*B;%Pn = KTB 1.38e-23 300
dmin = 222.1734;
f=[430,900,2400];
for i=1:3
    pt(i) = getpt(f(i),dmin, -100);%单位dbm
    C = B*log2(1+0.001*10^(-100/10)/pn);%bps 
    C1 = C*1.1;
    pr = 10*log10((2^(C1/B)-1)*pn/0.001);
    new_pt(i) = getpt(f(i),dmin,pr);
end
bar(f,0.001*10.^(new_pt/10),'r');
hold on;
bar(f,0.001*10.^(pt/10),'b');
xlabel('f/Mhz');
ylabel('pt/mW')
hold on;
