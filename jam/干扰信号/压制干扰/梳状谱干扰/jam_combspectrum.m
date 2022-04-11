function [noise_szp,t_noise] = jam_combspectrum(fs,Bj,Ns,fj,frame,Prj,Tr)
%========================================================================%
%功能：产生梳状谱干扰信号                                              %
%输入：
%   fs          采样频率   
%   Bj          噪声信号带宽
%   Ns           梳状谱个数
%   fj          噪声信号中心频率
%   frame       脉冲数
%   Prj         干扰功率
%   Tr          脉冲重复周期
%输出：
%   noise_szp    梳状谱噪声信号
%========================================================================%
% % 预设输入参数
% fs=160e6;
% Bj=20e6;
% Ns=3;
% frame=64;
% Prj=10;
% Tr=40e-6;
% fj=[1e7,4e7,7e7];
t=0:1/fs:frame*Tr-1/fs;N=length(t);
Bn=Bj;
sj_all=zeros(1,N);
u_all=sqrt(Prj)*randn(1,Ns*N);

for k=1:Ns
    u=u_all(N*(k-1)+1:N*k);
    
%频域截取
Fu=fftshift(fft(u));
n1=floor((fs-Bn)*N/2/fs);n2=floor((fs+Bn)*N/2/fs);
Fuj=[zeros(1,n1-1),Fu(n1:n2),zeros(1,N-n2)];

%从频域到时域
uj=ifft(ifftshift(Fuj));

%噪声调制
phi=2*pi*rand(1,1);
sj=uj.*exp(1i*(2*pi*fj(k)*t+phi));

sj_all=sj_all+sj;
end
 
noise_szp=sj_all;
t_noise=t;
end

