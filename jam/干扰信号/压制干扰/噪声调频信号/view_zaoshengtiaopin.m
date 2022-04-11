function [noise_tp] =view_zaoshengtiaopin(fs,Kfm,Pn,Prj,Bn,fj,frame,Tr)
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
% % 预设输入参数
% fs=160e6;%HZ
% Bj=20e6;%HZ
% fj=10e6;%HZ
% frame=64;
% Prj=40;%dBmW
% Tr=40e-6;%s
% Kfm=4e6;
% Bn=10e6;%HZ
% Pn=(Bj/2/(2.5*Kfm))^2;

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
figure(1);
subplot(311);
plot(t,real(uj));grid;xlabel('时间(s)');ylabel('幅度(V)');
title('调制噪声时域波形');
subplot(312);
[y,x]=ksdensity(real(uj));bar(x,y);grid;xlabel('幅度(V)');ylabel('概率密度');
title('调制噪声概率密度函数');
subplot(313);
window=hamming(100);noverlap=20;n=N;
yrj=pwelch(uj,window,noverlap,n,fs);
Pdb=10*log10(yrj);
frj=0:fs/length(Pdb):(length(Pdb)-1)*(fs/length(Pdb));
 plot(frj-fs/2,fftshift(Pdb));grid;xlabel('频率(Hz)');ylabel('功率谱(dB/Hz)');title('调制噪声功率谱密度');
 
% waitbar(.67,h,'干扰产生中，请稍等！');
% pause(1)

%频率调制
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
figure(2);
subplot(311);
plot(t,real(sj));grid;xlabel('时间(s)');ylabel('幅度(V)');
title('噪声调频信号时域波形');
subplot(312);
[y,x]=ksdensity(real(sj));bar(x,y);grid;xlabel('幅度(V)');ylabel('概率密度');
title('噪声调频信号概率密度函数');
subplot(313);
window=hamming(100);noverlap=20;n=N;
yrj=pwelch(sj,window,noverlap,n,fs);
Pdb=10*log10(yrj);
frj=0:fs/length(Pdb):(length(Pdb)-1)*(fs/length(Pdb));
 plot(frj-fs/2,fftshift(Pdb));grid;xlabel('频率(hz)');ylabel('功率谱(db/hz)');title('噪声调频功率谱密度');
 noise_tp=sj;
%  close(h);
end

