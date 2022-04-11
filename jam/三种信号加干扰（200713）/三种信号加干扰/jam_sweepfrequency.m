function [ sig_noise,t_noise ] = jam_sweepfrequency( fs,Bj,fj,frame,Prj,Tr,T_fr,Time_begin,K_sweep )
%========================================================================%
%功能：产生扫频干扰信号                                              %
%输入：
%   fs          采样频率   
%   Bj          噪声信号带宽
%   fj          噪声信号中心频率
%   frame       脉冲数
%   Prj         干扰功率
%   Tr          脉冲重复周期
%   T_fr        扫频周期    （最好在时间轴长度范围内，这样可以达到最大频率）
%   Time_begin  干扰开始时间（即当前时刻，信号的起始时间，可以为任意值）
%   K_sweep     中心频率变化斜率（K_sweep*T_fr应该与fj在一个数量级上，否则变化不明显）
%输出：
%   sig_noise    扫频噪声信号
%   t_noise      输出信号的时间轴
%========================================================================%

t=0:1/fs:frame*Tr-1/fs;N=length(t);
Bn=Bj;
u=sqrt(Prj)*randn(1,N);
%频域截取
Fu=fftshift(fft(u));
n1=floor((fs-Bn)*N/2/fs);n2=floor((fs+Bn)*N/2/fs);
Fuj=[zeros(1,n1-1),Fu(n1:n2),zeros(1,N-n2)];
%从频域到时域
uj=ifft(ifftshift(Fuj));

%噪声调制
F_fr=1/T_fr;
t_begin=Time_begin-floor(Time_begin*F_fr)*T_fr;      %t_begin 是当前时刻在一个扫频周期中的位置
t_f=t_begin+(0:N-1)/fs;                              %从扫频周期中的t_begin位置起，后面时间长度N的范围都对中心频率进行线性变化
tj=Time_begin+t;
sj=uj.*exp(1i*(2*pi*(fj+K_sweep*mod(t_f,T_fr)).*tj));%用mod是由于时间轴是开始时间是从t_begin开始而且在加上长度N*ts的时间有可能超出扫频周期，此时要把超出的部分折回从扫频周期的0时刻开始。
                                                     %扫频的最小值是fj，最大值是:（时间超过一个扫频周期）fj+K_sweep*T_fr；（时间小于一个扫频周期）;
sig_noise=sj;
t_noise=tj;
end

