function [noise_sp] = shepinzaosheng(fs,Bj,fj,frame,Prj,Tr)
%========================================================================%
%功能：产生射频噪声干扰信号                                              %
%输入：
%   fs          采样频率   
%   Bj          噪声信号带宽
%   fj          噪声信号中心频率
%   frame       脉冲数
%   Vj          干扰幅度
%   Tr          脉冲重复周期
%输出：
%   noise_sp    射频噪声信号
%========================================================================%
t=0:1/fs:frame*Tr-1/fs;N=length(t);
vRand=randn(1,2*N);
u=vRand(1:N)+1i*vRand(N+1:2*N);
Bn=Bj;

%频域截取
Fu=fftshift(fft(u));
n1=floor((fs-Bn)*N/2/fs);n2=floor((fs+Bn)*N/2/fs);
Fuj=[zeros(1,n1-1),Fu(n1:n2),zeros(1,N-n2)];
%从频域到时域
uj=ifft(ifftshift(Fuj));
Vj=dBmtoVwithoutr(Prj);
uj=PowerWeighted(Vj,uj);



%噪声调制
phi=2*pi*rand(1,1);
sj=uj.*exp(1i*(2*pi*fj*t+phi));


noise_sp=sj;
end

