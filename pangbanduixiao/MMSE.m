%����MMSE׼�������Ӧ�԰����  
clc;
clear all;close all
N=16;%��������Ԫ����
M=2;%������������
d_lamda=0.5;%����Ԫ������벨��
theta0=0;%����ָ��
thetaj=[20,40];%���ŷ���
theta=-60:1:60;%����ͼɨ��Ƕȷ�Χ
JNR=40;%�����
Ns=32;%���Ų�����
nj=length(thetaj);
epsilon=0.00001;
j=sqrt(-1);
Vs=exp(j*2*pi*d_lamda*(0:(N-1))'*sin(theta*pi/180));
Vs0=exp(j*2*pi*d_lamda*(0:(N-1))'*sin(theta0*pi/180));
Vsj=exp(j*2*pi*d_lamda*(0:(N-1))'*sin(thetaj*pi/180));
AJ=10^(JNR/20)*0.707*(randn(nj,Ns)+j*randn(nj,Ns));
noise=0.707*(randn(N,Ns)+j*randn(N,Ns));
Xs=Vsj*AJ+noise;%�����������ź�
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