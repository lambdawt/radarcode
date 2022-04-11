close all;
clear all;
clc; 
Tp = 20e-6;                  % 脉宽；
B = 10e6;                   % 带宽；
gama = B/Tp;                % 调频斜率
C = 3e8;                    % 光速；  
fs = 4*B;                  % 采样率； 
R_wave_gate = 10e3;         % 波门中心；
T_wave_gate = Tp*5;         % 波门时宽；
nrn = round(fs*T_wave_gate);% 距离向采样点数；
Tstart = 2*R_wave_gate/C - nrn/2/fs;    % 距离向采样起始时刻；
Tend = 2*R_wave_gate/C + nrn/2/fs - 1/fs;   % 距离向采样终止时刻；
tnrn = Tstart:1/fs:Tend;    % 距离向采样时刻；


%%-------目标参数
tar_num = 2;                % 目标个数；
R_tar = [10e3 12e3];
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

figure;
plot(tnrn*1e6, abs(S_tar));
xlabel('\fontsize{10.5}\fontname{宋体}时间\fontname{Times New Roman}/us', 'FontSize', 10.5);
ylabel('\fontsize{10.5}\fontname{宋体}幅度\fontname{Times New Roman}/V', 'FontSize', 10.5);


%%%%% 生成密集多假目标 %%%%%
% L=2e-6; % 假目标间隔
% M=20;  % 假目标个数
% k=T/L; % 分段数目
% duanNum=k+M-1;  % 重构后总段数
% T_g=duanNum*L;  % 重构后时间长度
% gn=zeros(1,T_g*Fs);
% gn(1,1:N)=sn;
% for i=1:M-1
%     cur = gn(1,i*L*Fs+1:i*L*Fs+N);
%     gn(1,i*L*Fs+1:i*L*Fs+N)=cur+A*exp(1j*(2*pi*F0*(t)+pi*K*(t).^2+Phi0));
% end
% 
% %%%%% 对密集多假目标干扰信号进行脉压 %%%%%
% S_rec = gn;
% S_ref = A*exp(1j*(2*pi*F0*t+pi*K*t.^2+Phi0));
% temp = conv(S_rec,conj(fliplr(S_ref)));
% S_pc = temp(length(S_ref):end);
%                                      
% hfig = figure(1);
% subplot(3,1,1);
% plot((t+T/2)*1e6,real(sn));
% title('\fontsize{10.5}\fontname{宋体}线性调频信号实部');
% set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10.5);
% grid on;
% axis tight;
% xlabel('\fontsize{10.5}\fontname{宋体}时间\fontname{Times New Roman}/us', 'FontSize', 10.5);
% ylabel('\fontsize{10.5}\fontname{宋体}幅度\fontname{Times New Roman}/V', 'FontSize', 10.5)
% 
% subplot(3,1,2);
% plot(linspace(0,T_g,T_g*Fs)*1e6,abs(gn));
% title('\fontsize{10.5}\fontname{宋体}密集多假目标干扰信号包络');
% set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10.5);
% grid on;
% axis tight;
% xlabel('\fontsize{10.5}\fontname{宋体}时间\fontname{Times New Roman}/us', 'FontSize', 10.5);
% ylabel('\fontsize{10.5}\fontname{宋体}幅度\fontname{Times New Roman}/V', 'FontSize', 10.5);
% 
% subplot(3,1,3);
% plot(linspace(0,T_g,T_g*Fs)*1e6,abs(S_pc));
% title('\fontsize{10.5}\fontname{宋体}密集多假目标干扰信号脉压结果');
% set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10.5);
% grid on;
% axis tight;
% xlabel('\fontsize{10.5}\fontname{宋体}时间\fontname{Times New Roman}/us', 'FontSize', 10.5);
% ylabel('\fontsize{10.5}\fontname{宋体}幅度\fontname{Times New Roman}/V', 'FontSize', 10.5);