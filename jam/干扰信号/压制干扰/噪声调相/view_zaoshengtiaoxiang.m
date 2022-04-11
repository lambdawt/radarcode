function [noise_tx] =view_zaoshengtiaoxiang(fs,Bj,fj,Prj,Tr,frame)
%========================================================================%
%功能：产生噪声调相干扰信号                                              %
%输入：
%   fs          采样频率  
%   Bj          干扰带宽
%   fj          噪声信号中心频率
%   Prj         干扰功率（单位dBm）
%   Tr          脉冲重复周期
%   frame       脉冲数
%输出：
%   noise_tx    噪声调相信号
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
%调制噪声产生
u=2*pi*rand(1,N)-pi;
waitbar(.33,h,'干扰产生中，请稍等！');
pause(1)

figure(1);
subplot(311);
plot(t,real(u));grid;xlabel('时间(s)');ylabel('幅度(V)');
title('调制噪声时域波形');

subplot(312);
[y,x]=ksdensity(u);bar(x,y);grid;xlabel('幅度(V)');ylabel('概率密度');
title('调制噪声概率密度函数');

subplot(313);
window=hamming(100);noverlap=20;n=N;
yrj=pwelch(u,window,noverlap,n,fs);
Pdb=10*log10(yrj);
frj=0:fs/length(Pdb):(length(Pdb)-1)*(fs/length(Pdb));
plot(frj-fs/2,fftshift(Pdb));grid;xlabel('频率(Hz)');ylabel('功率谱(dB/Hz)');title('调制噪声功率谱密度');
 
waitbar(.67,h,'干扰产生中，请稍等！');
pause(1)
 %噪声调相
 phi=2*pi*rand(1,1);
 sj=zeros(1,N);
 for i=1:N
 sj(i)=exp(1i*(2*pi*fj*t(i)+u(floor(t(i)*Bj)+1)+phi));
 end
 Vj=dBmtoV(Prj);
 sj=PowerWeighted(Vj,sj);
figure(2);
subplot(311);
plot(t,real(sj));grid;xlabel('时间(s)');ylabel('幅度(V)');
title('噪声调相信号时域波形');

subplot(312);
[y,x]=ksdensity(real(sj));bar(x,y);grid;xlabel('幅度(V)');ylabel('概率密度');
title('噪声调相信号概率密度函数');

subplot(313);
window=hamming(100);noverlap=20;n=N;
yrj=pwelch(sj,window,noverlap,n,fs);
Pdb=10*log10(yrj);
frj=0:fs/length(Pdb):(length(Pdb)-1)*(fs/length(Pdb));
 plot(frj-fs/2,fftshift(Pdb));grid;xlabel('频率(hz)');ylabel('功率谱(db/hz)');title('噪声调相功率谱密度');

 noise_tx=sj;
  close(h);
end

