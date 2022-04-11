function [noise_tp] =zaoshengtiaopin(fs,Kfm,Pn,Prj,Bn,fj,frame,Tr)
%========================================================================%
%功能：产生噪声调频干扰信号                                              %
%输入：
%   fs          采样频率  
%   Kfm         调频斜率
%   Prj          噪声调频信号功率
%   Bj          噪声调频信号的等效干扰带宽
%   Bn          调制噪声带宽
%   fj          噪声信号中心频率
%   frame       脉冲数
%   Tr          脉冲重复周期
%输出：
%   noise_tp    噪声调频信号
%========================================================================%
t=0:1/fs:frame*Tr-1/fs;
N=length(t);
%产生调制噪声（窄带高斯白噪声）
vGsRand=randn(1,2*N);
u=vGsRand(1:N)+1i*vGsRand(N+1:2*N);
Pn0=sqrt(fs*Pn/Bn);
u=Pn0*u;
%%%频域截取
Fu=fftshift(fft(u));
n1=floor((fs-Bn)*N/2/fs);n2=floor((fs+Bn)*N/2/fs);%
Fuj=[zeros(1,n1-1),Fu(n1:n2),zeros(1,N-n2)];
uj=ifft(ifftshift(Fuj));

phi=2*pi*rand(1,1);
%%%调制噪声积分
Uj=zeros(1,N);
for i=1:N
    if i==1
        Uj(i)=real(uj(i))/fs;
    else
    Uj(i)=Uj(i-1)+real(uj(i))/fs; 
    end
end
sj=exp(1i*(2*pi*fj*t+2*pi*Kfm*Uj+phi));
Vj=dBmtoV(Prj);
sj=PowerWeighted(Vj,sj);

 noise_tp=sj;
end

