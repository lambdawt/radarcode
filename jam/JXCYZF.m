% 2021-05-14 
% 灵巧噪声干扰

% 2022-02-18 
% 间歇采样转发 直接 转发干扰 
% 仿真分析；

% 2022-02-19
% 间歇采样转发 重复 转发干扰
% 仿真分析；

close all;
clear;
clc;
%% 
Tp = 20e-6;                  % 脉宽；
B = 10e6;                   % 带宽；
gama = B/Tp;



C = 3e8;                    % 光速；  
fs = 4*B;                  % 采样率； 
R_wave_gate = 10e3;         % 波门中心；
T_wave_gate = Tp*2;         % 波门时宽；
nrn = round(fs*T_wave_gate);% 距离向采样点数；
Tstart = 2*R_wave_gate/C - nrn/2/fs;    % 距离向采样起始时刻；
Tend = 2*R_wave_gate/C + nrn/2/fs - 1/fs;   % 距离向采样终止时刻；
tnrn = Tstart:1/fs:Tend;    % 距离向采样时刻；


%%-------目标参数
tar_num = 1;                % 目标个数；
% R_tar = [10e3 10e3+500 10e3-100];
R_tar = [8e3 10e3+500 10e3-100];
R_tar = R_tar(1:tar_num);   % 目标径向距离（初始时刻）；
V_tar = [50 100 150 200];
V_tar = V_tar(1:tar_num);   % 目标径向速度（初始时刻）；目标运动状态：假设匀速直线运动；

tao_tar = 2*R_tar/C;        % 第一个脉冲，目标的时延；

%% ----目标回波
S_tar = zeros(1,nrn);
for num = 1:tar_num
    win = (tnrn-tao_tar(1,num)-Tp/2)>-Tp/2 & (tnrn-tao_tar(1,num)-Tp/2)<Tp/2;
    S_tar = S_tar + win.*exp(1j*pi*gama*(tnrn-tao_tar(1,num)-Tp/2).^2);
end
%------参考信号（雷达发射信号）；
S_ref = exp(1j*pi*gama*(-Tp/2:1/fs:(Tp/2-1/fs)).^2);  


%% ------间歇采样干扰生成（直接转发or重复转发，通过repeat_num这个参数来控制）
% 重复转发干扰，参数1；
jam_amp = 1;
sampl_time = 2e-6;     % 采样时长；
sampl_period = 6e-6;   % 采样周期；
ISRJ_signal_1 = ISRJ(S_tar,jam_amp,sampl_time,sampl_period,fs);

% figure,subplot(2,1,1),plot(real(S_tar));title('雷达发射信号');
% subplot(2,1,2),plot(real(ISRJ_signal_1));title('间歇采样重复转发干扰信号');

% 论文绘图
[c, r] = size(S_tar);
x = linspace(0, r/fs, r) * 1e3;% ms
F = linspace(-fs/2, fs/2, r);% Hz
hfig = figure(1);
subplot(3, 1, 1);
plot(x*1e3, real(S_tar));
title('\fontsize{10.5}\fontname{宋体}雷达发射信号时域波形');
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10.5);
grid on;
axis tight;
xlabel('\fontname{宋体}时间\fontname{Times New Roman}/us', 'FontSize', 10.5);
ylabel('\fontname{宋体}幅度\fontname{Times New Roman}/V', 'FontSize', 10.5);

subplot(3, 1, 2);
plot(x*1e3, real(ISRJ_signal_1));
title('\fontsize{10.5}\fontname{宋体}间歇采样转发干扰信号时域波形');
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10.5);
grid on;
axis tight;
xlabel('\fontname{宋体}时间\fontname{Times New Roman}/us', 'FontSize', 10.5);
ylabel('\fontname{宋体}幅度\fontname{Times New Roman}/V', 'FontSize', 10.5);

subplot(3, 1, 3);
plot(F/1e6, abs(fftshift(fft(ISRJ_signal_1))));
% plot(F/1e6, db(abs(fftshift(fft(ISRJ_signal_1))) / max(abs(fftshift(fft(ISRJ_signal_1))))));
title('\fontsize{10.5}\fontname{宋体}间歇采样转发干扰信号频谱');  
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10.5);
grid on;
% axis tight;
xlabel('\fontname{宋体}频率\fontname{Times New Roman}/MHz', 'FontSize', 10.5);
ylabel('\fontname{宋体}幅度\fontname{Times New Roman}/V', 'FontSize', 10.5);


% 设置图片输出格式
% figWidth = 14;
% figHeight = 12.9;  % 单栏建议14cm,8.6cm，黄金比例为0.618
% set(hfig, 'PaperUnits', 'centimeters');  %centimeters %inches
% set(hfig, 'PaperPosition', [0 0 figWidth figHeight]);
% print(hfig, ['间歇采样转发干扰.', 'tif'], '-r600', '-dtiff');   