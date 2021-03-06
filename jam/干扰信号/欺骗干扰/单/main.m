load('data_LFMParameter.mat')
load('data_target0Parameter.mat')
load('data_DeceptionJammingParameter1.mat')
temp1=3;%1为速度、2为距离、3为联合
c=3e8;%光速
R1=3000;%假目标距离
v1=300;%假目标速度
congmubiao=1;%假目标个数
sigma=2;%雷达目标截面积
% R=2500;%目标距离
% v=250;%假目标速度
% sigma=10;%目标雷达截面积
% tr_LFM=40e-6;%雷达信号重复周期
% fz_LFM=1e10;%雷达载频
% tau_LFM=1e-6;%雷达信号脉宽
% Pt_LFM=4.6e5;%雷达峰值功率
% Gt_LFM=30;%雷达发射天线增益
% Gr_LFM=30;%雷达接受天线增益
% L_LFM=3;%雷达接收机损失
% B1_LFM=20e6;%雷达信号调制带宽
% fs_LFM=160e9;%采样频率
% f0_LFM=20e6;%雷达信号中心频率
% F_LFM=4;%雷达接收机噪声系数
% B_LFM=20e6;%雷达接收机带宽
% Te_LFM=290;%雷达接收机等效噪声温度
% frame_LFM=32;%脉冲积累个数
% Pj=9.5;%干扰机发射功率
% Gj=30;%干扰机发射天线增益
% Gjr=10;%干扰机发射天线增益
fr=1/tr_LFM;%脉冲重复频
lamta=c/fz_LFM;%波长
tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%一个脉冲重复周期采样序列
N=length(tm);%一个脉冲重复周期采样点数长度
An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%噪声强度
ts=1/fs_LFM;
k=B1_LFM/tau_LFM;   
Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %#ok<*NODEF> %目标回波信号功率
A=sqrt(Prs);%回波信号幅度
f_doppler=2*v/lamta;%真目标多普勒频率%线性调频信号调制系数
Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L_LFM/10)));
Aj=sqrt(Prj);
f_doppler1=2*v1/lamta;
[y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
[~,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y,tr_LFM,ts,N);
[s_echo_2,echo]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
[s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
[s_ft,echo3]=LFMDeceptionJamming(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao,temp1);
s_echo_1=s_echo_2+s_noise+s_ft;%%%%干扰信号
figure(1),plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('雷达接收信号');
figure(2),plot((0:fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_LFM-fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('雷达接收信号的频谱');
