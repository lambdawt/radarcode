function [s_echo_1]=gaofang(f0,B,fs,s_echo_1)
n=50;
Wn=[(f0-B/2)/(fs/2),(f0+B/2)/(fs/2)];         % 通带范围10M～30M  
b=fir1(n,Wn);         % n阶带通滤波器系数
% s_echo=[s_echo,zeros(1,25)];
s_echo_1=filter(b,1,s_echo_1); 