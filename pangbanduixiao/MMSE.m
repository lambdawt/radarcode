%基于MMSE准则的自适应旁瓣对消  
clc;
clear all;close all
N=16;%主天线阵元个数
M=2;%辅助天线阵列
d_lamda=0.5;%阵列元间隔比与波长
theta0=0;%波束指向
thetaj=[20,40];%干扰方向
theta=-60:1:60;%方向图扫描角度范围
JNR=40;%干燥比
Ns=32;%干扰采样数
nj=length(thetaj);
epsilon=0.00001;
j=sqrt(-1);
Vs=exp(j*2*pi*d_lamda*(0:(N-1))'*sin(theta*pi/180));
Vs0=exp(j*2*pi*d_lamda*(0:(N-1))'*sin(theta0*pi/180));
Vsj=exp(j*2*pi*d_lamda*(0:(N-1))'*sin(thetaj*pi/180));
AJ=10^(JNR/20)*0.707*(randn(nj,Ns)+j*randn(nj,Ns));
noise=0.707*(randn(N,Ns)+j*randn(N,Ns));
Xs=Vsj*AJ+noise;%天线主干扰信号
Xj=Xs(1:M,: );
D=Vs0'*Xs;
R11=Xj*Xj'/Ns;
r01=Xj*D'/Ns;
W=R11\r01;
pattern1=abs(Vs0'*Vs-W'*Vs(1:M,: ))+epsilon;
pattern1=20*log10(pattern1/max(pattern1));
pattern0=abs(Vs0'*Vs)+epsilon;
pattern0=20*log10(pattern0/max(pattern0));
plot(theta,pattern0,'r--',theta,pattern1);
grid on;
xlabel('azimuth/\circ');
ylabel('gain/db');