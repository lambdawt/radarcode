% 相干转发
%均匀线阵
clc;
clear;
close all;

%% 常数
c = 3e8;          % 光速（m/s）
k = 1.38*1e-23;     %玻尔兹曼常数（J/K）
To = 290;               % 标准室温 17℃（K 开尔文）.

%% 雷达系统参数
M=30;
fo = 1.2*1e9;         % 载频(Hz)
Pt = 15*1e3;          % 发射机峰值功率(w）
Gt = 33;              % 发射机增益（dB）
B = 10*1e6;            % 接收机带宽
fs = 2*B;             %采样频率
fr = 2000;            % 脉冲重复频率(Hz)
Tr = 1/fr;            % 脉冲重复周期 （秒）
Tp = 2e-5;              % 脉宽 （秒）
lambda = c/fo;           % 波长 （m）
d = lambda/2;            % 阵元间距 （m）
Gele = 5 ;              %阵元接收增益（dB）
F = 12;                 %接收机噪声系数（dB）
Ts = To*(10^(F/10)-1);      %接收机温度 (K)
Nn = k*Ts;                 % 接收机噪声功率谱密度 （W/Hz）
sigmaN2 = 1;               % 噪声功率（W）
dR  = c/2/B;               % 距离增量 （m）    

%**************************原始信号形式
K = B/Tp;  %线性调频斜率
t0 = 0:1/fs:Tp;
s0 = exp(1i*pi*K*t0.^2);
t = [t0,Tp+1/fs:1/fs:Tr-1/fs];
s = exp(1i*pi*K*t.^2).*rectpuls(t-Tp/2,Tp);%单个线性调频脉冲信号
Ns = length(t);            %一周期采样点数

figure(1)
subplot(2,1,1)
plot(t,real(s));
subplot(2,1,2)
plot(t(1:Tp/Tr*Ns),real(s(1:Tp/Tr*Ns)))
%****************************目标回波
%% 目标信号
R_t = 35e3; 
vt = 35;
fdt = 2*vt/lambda;
tao = 2*R_t/c;
K = B/Tp;  %线性调频斜率
sr = exp(1i*pi*K*(t-tao).^2).*rectpuls(t-Tp/2-tao,Tp);

figure(2)
subplot(2,1,1)
plot(t,real(sr));
subplot(2,1,2)
plot(t(tao/Tr*Ns:(tao+Tp)/Tr*Ns),real(sr(tao/Tr*Ns:(tao+Tp)/Tr*Ns)));

%匹配滤波器
s0 = exp(1i*pi*K*t0.^2);
h = conj(fliplr(s0));
H = fft(h,Ns);

S0 = repmat(sr,1,M);
SigAll = S0.*exp(1i*2*pi*fdt*(0:Ns*M-1)/fs);
SigAll= reshape(SigAll,Ns,M);

%% 间歇采样转发干扰

TTs = 4e-6;       %采样信号pri
tau = TTs/2;      %采样信号Tp
tau1 = TTs;
pt0 = square(2*pi*(t+tau/2)/tau1,50); %生成周期方波信号
pt= (pt0+1)/2;       %生成采样信号

figure(11);
plot(t,real(pt));

s11 = sr.*pt;
s11=[zeros(1,ceil(k*tau*fs)),s11(1:Ns-ceil(k*tau*fs))];

% s11 = s11*1;
S1 = repmat(s11,1,M);
JamAll = S1.*exp(1i*2*pi*fdt*(0:Ns*M-1)/fs);
JamAll= reshape(JamAll,Ns,M);


figure(3)
subplot(2,1,1)
plot(t,real(sr),'k');
xlabel('时间（s）','FontSize',19,'FontName','宋体','Color','k');
ylabel('幅度（V）','FontSize',19,'FontName','宋体','Color','k');
subplot(2,1,2)
plot(t,real(JamAll(:,1)),'k');
xlabel('时间（s）','FontSize',19,'FontName','宋体','Color','k');
ylabel('幅度（V）','FontSize',19,'FontName','宋体','Color','k');
% hold on
% plot(t,real(SigAll(:,1)),'r:')
% ylim([-5 5]);
% figure(2)
% Yt=fftshift(fft(SigAll(:,1)));
% f=linspace(-fs/2,fs/2,length(t));
% plot(f,abs(Yt))


Jam_pp = zeros(Ns,M);
Sig_pp = zeros(Ns,M);
for m=1:M
    Jam_fft = fft(JamAll(:,m).');
    Sig_fft = fft(SigAll(:,m).');
    Jam_pp(:,m) = ifft(H.*Jam_fft);
    Sig_pp(:,m) = ifft(H.*Sig_fft);
end

% save('jxcy_zr','Jam_pp','Sig_pp');

mtd1 = fftshift(fft(Jam_pp,16,2),2);
mtd1 = mtd1/max(max(abs(mtd1)));
mtd2 = fftshift(fft(Sig_pp,16,2),2);
mtd2 = mtd2/max(max(abs(mtd2)));
figure(4);
mesh(abs(mtd1));
hold on
mesh(abs(mtd2));
shading interp
colormap jet
axis tight
xlabel('多普勒通道','FontSize',19,'FontName','宋体','Color','k');
ylabel('距离门','FontSize',19,'FontName','宋体','Color','k');
zlabel('归一化幅度','FontSize',19,'FontName','宋体','Color','k');

All = [Sig_pp(:,1);Jam_pp(:,1)];
figure(5)
plot(t,abs(Sig_pp(:,1))/max(abs(All)),'r--');
hold on
plot(t,abs(Jam_pp(:,1)/max(abs(All))),'k');
% axis([1.8e-4 2.5e-4 0 1200])
% legend('目标回波','间歇采样干扰','FontSize',19)
xlabel('时间（s）','FontSize',19,'FontName','宋体','Color','k');
ylabel('幅度（V）','FontSize',19,'FontName','宋体','Color','k')










