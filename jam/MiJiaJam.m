close all;
clear all;
clc; 
%%信号参数
A=1;                  %发射信号的振幅
Phi0=0;               %发射信号的随机初相
T=20e-6;              %信号时宽
B=10e6;               %信号带宽
F0=0e6;               %中频频率，即载频频率
Fs=40e6;        %采样频率

%%%%% 信号的参数设置 %%%%%
K=B/T;         %调频斜率
Ts=1/Fs;       %采样周期
N=T/Ts;        %采样点数

%%%%% 生成线性调频发射信号 %%%%%
t=linspace(-T/2,T/2,N);
freq=linspace(-Fs/2,Fs/2,N);
sn=A*exp(1j*(2*pi*F0*t+pi*K*t.^2+Phi0));

%%%%% 生成密集多假目标 %%%%%
L=2e-6; % 假目标间隔
M=20;  % 假目标个数
k=T/L; % 分段数目
duanNum=k+M-1;  % 重构后总段数
T_g=duanNum*L;  % 重构后时间长度
gn=zeros(1,T_g*Fs);
gn(1,1:N)=sn;
for i=1:M-1
    cur = gn(1,i*L*Fs+1:i*L*Fs+N);
    gn(1,i*L*Fs+1:i*L*Fs+N)=cur+A*exp(1j*(2*pi*F0*(t)+pi*K*(t).^2+Phi0));
end

%%%%% 对密集多假目标干扰信号进行脉压 %%%%%
S_rec = gn;
S_ref = A*exp(1j*(2*pi*F0*t+pi*K*t.^2+Phi0));
temp = conv(S_rec,conj(fliplr(S_ref)));
S_pc = temp(length(S_ref):end);
                                     
hfig = figure(1);
subplot(3,1,1);
plot((t+T/2)*1e6,real(sn));
title('\fontsize{10.5}\fontname{宋体}线性调频信号实部');
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10.5);
grid on;
axis tight;
xlabel('\fontsize{10.5}\fontname{宋体}时间\fontname{Times New Roman}/us', 'FontSize', 10.5);
ylabel('\fontsize{10.5}\fontname{宋体}幅度\fontname{Times New Roman}/V', 'FontSize', 10.5)

subplot(3,1,2);
plot(linspace(0,T_g,T_g*Fs)*1e6,abs(gn));
title('\fontsize{10.5}\fontname{宋体}密集多假目标干扰信号包络');
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10.5);
grid on;
axis tight;
xlabel('\fontsize{10.5}\fontname{宋体}时间\fontname{Times New Roman}/us', 'FontSize', 10.5);
ylabel('\fontsize{10.5}\fontname{宋体}幅度\fontname{Times New Roman}/V', 'FontSize', 10.5);

subplot(3,1,3);
plot(linspace(0,T_g,T_g*Fs)*1e6,abs(S_pc));
title('\fontsize{10.5}\fontname{宋体}密集多假目标干扰信号脉压结果');
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10.5);
grid on;
axis tight;
xlabel('\fontsize{10.5}\fontname{宋体}时间\fontname{Times New Roman}/us', 'FontSize', 10.5);
ylabel('\fontsize{10.5}\fontname{宋体}幅度\fontname{Times New Roman}/V', 'FontSize', 10.5);

% 设置图片输出格式
figWidth = 14;
figHeight = 12.9;  % 单栏建议14cm,8.6cm，黄金比例为0.618
set(hfig, 'PaperUnits', 'centimeters');  %centimeters %inches
set(hfig, 'PaperPosition', [0 0 figWidth figHeight]);
print(hfig, ['密集多假目标干扰.', 'tif'], '-r600', '-dtiff');   
% print(hfig, ['密集多假目标干扰.', 'jpg'], '-r600', '-djpeg');  