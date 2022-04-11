function [s_noise]=zaosheng(frame,N,An,B,fs,ts)
%%  输入 
%   frame   仿真脉冲数
%   N       一个脉冲重复周期采样点数长度
%   An      噪声强度
%   B       带宽
%   fs      采样频率
%%  输出
%   s_noise  噪声信号
%%  功能  产生噪声信号
%%
s_noise=wgn(1,frame*N,An);
%hold on,plot(0:ts:(frame*N-1)*ts,abs(s_noise),'g');
%figure,plot(0:fs/length(abs(fft(s_noise))):fs-fs/length(abs(fft(s_noise))),abs(fft(s_noise)));
window=chebwin(51,40);
[b,a]=fir1(50,2*B/fs,window);
s_noise=filter(b,a,s_noise);