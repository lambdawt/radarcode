function [noise_sp] = view_shepinzaosheng(fs,Bj,fj,frame,Prj,tr)
%========================================================================%
%功能：产生射频噪声干扰信号                                              %
%输入：
%   fs          采样频率   
%   Bj          噪声信号带宽
%   fj          噪声信号中心频率
%   frame       脉冲数
%   Vj          干扰幅度
%   tr          脉冲重复周期
%输出：
%   noise_sp    射频噪声信号
%========================================================================%
h = waitbar(0, '干扰产生中，请稍等！');
pause(.5)
t=0:1/fs:frame*tr-1/fs;N=length(t);
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
waitbar(.33,h, '干扰产生中，请稍等！');
pause(.1)

figure(1);
subplot(311);
plot(t,real(u));grid;xlabel('时间(s)');ylabel('幅度(V)');
title('调制噪声时域波形');

subplot(312);
% [y,x]=hist(real(u),100);y=y./N./mean(diff(x));bar(x,y,1);grid;xlabel('幅度(V)');ylabel('概率密度');
% title('调制噪声概率密度函数');
[y,x]=ksdensity(real(u));bar(x,y);grid;xlabel('幅度(V)');ylabel('概率密度');title('调制噪声概率密度函数');
%功率谱

window=hamming(100);noverlap=20;n=length(t);
yrj=pwelch(u,window,noverlap,n,fs);
Pdb=10*log10(yrj);
frj=0:fs/length(Pdb):(length(Pdb)-1)*(fs/length(Pdb));
subplot(313)
plot(frj-fs/2,fftshift(Pdb));grid;xlabel('频率(Hz)');ylabel('功率谱(dBW/Hz)');title('调制噪声功率谱密度');
waitbar(.67,h, '干扰产生中，请稍等！');
pause(1)

%噪声调制
phi=2*pi*rand(1,1);
sj=uj.*exp(1i*(2*pi*fj*t+phi));
figure(2);
subplot(311);
plot(t,real(sj));grid;xlabel('时间(s)');ylabel('幅度(V)');
title('射频噪声信号时域波形');
subplot(312);
% [y,x]=hist(real(sj),100);y=y./N./mean(diff(x));bar(x,y,1);grid;xlabel('幅度(V)');ylabel('概率密度');
% title('射频噪声信号概率密度函数');
[y,x]=ksdensity(real(sj));bar(x,y);grid;xlabel('幅度(V)');ylabel('概率密度');title('射频噪声信号概率密度函数');
subplot(313);
window=hamming(100);noverlap=20;n=N;
yrj=pwelch(sj,window,noverlap,n,fs);
Pdb=10*log10(yrj);
frj=0:fs/length(Pdb):(length(Pdb)-1)*(fs/length(Pdb));
 plot(frj-fs/2,fftshift(Pdb));grid;xlabel('频率(Hz)');ylabel('功率谱(dBW/Hz)');title('射频噪声信号功率谱密度');

noise_sp=sj;
close(h);
end

