function [s_echo_1]=gaofang(f0,B,fs,s_echo_1)
n=50;
Wn=[(f0-B/2)/(fs/2),(f0+B/2)/(fs/2)];         % ͨ����Χ10M��30M  
b=fir1(n,Wn);         % n�״�ͨ�˲���ϵ��
% s_echo=[s_echo,zeros(1,25)];
s_echo_1=filter(b,1,s_echo_1); 