function [noise_tx] =zaoshengtiaoxiang(fs,Bj,fj,Prj,Tr,frame)
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
t=0:1/fs:frame*Tr-1/fs;
N=length(t);
%调制噪声产生
u=2*pi*rand(1,N)-pi;


 %噪声调相
 phi=2*pi*rand(1,1);
 sj=zeros(1,N);
 for i=1:N
 sj(i)=exp(1i*(2*pi*fj*t(i)+u(floor(t(i)*Bj)+1)+phi));
 end
 Vj=dBmtoV(Prj);
 sj=PowerWeighted(Vj,sj);

 noise_tx=sj;
end

