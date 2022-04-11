function [noise_tf] =view_zaoshengtiaofu(fs,Bj,fj,Prj,Tr,frame)
%========================================================================%
%功能：产生噪声调幅干扰信号                                              %
%输入：
%   fs          采样频率  
%   Bj          干扰带宽
%   fj          噪声信号中心频率
%   Prj         干扰功率
%   Tr          脉冲重复周期
%   frame       脉冲数
%输出：
%   noise_tf    噪声调频信号
%========================================================================%
% % 预设输入参数
% fs=160e6;
% Bj=20e6;
% fj=10e6;
% frame=64;
% Prj=40;%dBmW
% Tr=40e-6;
h = waitbar(0,'干扰产生中，请稍等！');
pause(.5)

t=0:1/fs:frame*Tr-1/fs;
N=length(t);
Bn=Bj;
vGsRandnum=randn(1,2*N);
u=vGsRandnum(1:N)+1i*vGsRandnum(N+1:2*N);
Prj=dBmtoW(Prj);
Pn=Prj*2/19;
Pn0=sqrt(fs*Pn/Bj);
u=Pn0*u;
%%%频域截取
Fu=fftshift(fft(u));
n1=floor((fs-Bn)*N/2/fs);n2=floor((fs+Bn)*N/2/fs);%
Fuj=[zeros(1,n1-1),Fu(n1:n2),zeros(1,N-n2)];
uj=ifft(ifftshift(Fuj));

waitbar(.33,h,'干扰产生中，请稍等！');
pause(1)

figure(1);
subplot(311);
plot(t,real(uj));grid;xlabel('时间(s)');ylabel('幅度(V)');
title('调制噪声时域波形');

subplot(312);
[y,x]=ksdensity(real(uj));bar(x,y);grid;xlabel('幅度(V)');ylabel('概率密度');title('调制噪声概率密度函数');

subplot(313);
window=hamming(100);noverlap=20;n=N;
yrj=pwelch(uj,window,noverlap,n,fs);
Pdb=10*log10(yrj);
frj=0:fs/length(Pdb):(length(Pdb)-1)*(fs/length(Pdb));
 plot(frj-fs/2,fftshift(Pdb));grid;xlabel('频率(Hz)');ylabel('功率谱(dB/Hz)');title('调制噪声功率谱密度');
waitbar(.67,h,'干扰产生中，请稍等！');
pause(1)
 
%幅度调制
 phi=2*pi*rand(1,1);
 U0=sqrt(Prj-Pn/2);
 sj=(U0+real(uj)).*exp(1i*(2*pi*fj*t+phi));
 figure(2);
subplot(311);
plot(t,real(sj));grid;xlabel('时间(s)');ylabel('幅度(V)');
title('噪声调幅信号时域波形');
subplot(312);
[y,x]=ksdensity(real(sj));bar(x,y);grid;xlabel('幅度(V)');ylabel('概率密度');
title('噪声调幅信号概率密度函数');
subplot(313);
window=hamming(100);noverlap=20;n=N;
yrj=pwelch(sj,window,noverlap,n,fs);
Pdb=10*log10(yrj);
frj=0:fs/length(Pdb):(length(Pdb)-1)*(fs/length(Pdb));
 plot(frj-fs/2,fftshift(Pdb));grid;xlabel('频率(hz)');ylabel('功率谱(db/hz)');title('噪声调幅信号功率谱密度');
 
 noise_tf=sj;
 close(h);
end

