function [M,match_filter_fft]=maiyaxishu(f0,fs,y,tr,ts,N)
%%  输入
%   f0    中心频率
%   fs    采样频率
%   y     发射信号
%   tr    脉冲重复周期
%   ts    采样周期
%   N     一个脉冲重复周期采样点数长度
%%  输出
%   M     脉冲采样点数（一个PRI内）
%   match_filter_fft     匹配滤波器系数
%%  功能   产生脉冲压缩系数（匹配滤波器频率响应）
%%


n=0:N-1;%调频信号时间序列的长度
local_oscillator_i=cos(n*(-f0)/fs*2*pi);
local_oscillator_q=sin(n*(-f0)/fs*2*pi);
fbb_i=local_oscillator_i.*y;
fbb_q=local_oscillator_q.*y;
% window=chebwin(51,40);%切比雪夫窗
% [b,a]=fir1(50,2*B/fs,window);%低通滤波器，截止频率为2*B/fs
% fbb_i=[fbb_i,zeros(1,25)];
% fbb_q=[fbb_q,zeros(1,25)];
% fbb_i=filter(b,a,fbb_i);%一维数字滤波器，a为分母，b为分子
% fbb_q=filter(b,a,fbb_q);
% fbb_i=fbb_i(26:end);
% fbb_q=fbb_q(26:end);
fbb=fbb_i+1j*fbb_q;
% figure,plot(ft,fbb_i),xlabel('t(单位：秒)'),ylabel('y(单位：伏)'),title('解调后的I路信号');
% figure,plot(ft,fbb_q),xlabel('t(单位：秒)'),ylabel('y(单位：伏)'),title('解调后的Q路信号');
% figure,plot(0:fs/length(abs(fft(fbb))):fs-fs/length(abs(fft(fbb))),abs((fft(fbb)))),xlabel('f(单位：Hz)'),ylabel('y'),title('解调后信号的频谱');
%%%%%%%%%%%%%%%%%%%产生理想脉压系数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=fix(tr/ts);%2^nextpow2(tr/ts);%M为接近PRI长度的2的K次方
match_filter=fliplr(conj(fbb));%conj为求共轭，flipr为求矩阵翻转压缩后幅度增加sqrt(D)倍，FFT要求得真实幅度，应当乘以2除以n
% win=hamming(O)';
% match_filter=match_filter.*win;
match_filter_fft_1=fft(match_filter,M);
match_filter_fft=match_filter_fft_1.';
% figure,plot(real(match_filter_fft)),title('脉冲压缩系数(实部）');
% figure,plot(imag(match_filter_fft)),title('脉冲压缩系数(虚部）');