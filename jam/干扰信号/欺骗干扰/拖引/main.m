load('data_LFMParameter.mat')
load('data_target0Parameter.mat')
load('data_DeceptionJammingParameter.mat')
temp1=2;%1为速度、2为距离、3为联合
c=3e8;
lamta=c/fz_LFM;
Rj=2e3;
Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L_LFM/10)));
Aj=sqrt(Prj);
tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%一个脉冲重复周期采样序列
N=length(tm);%一个脉冲重复周期采样点数长度
f_doppler=2*v/lamta;
f_doppler1=f_doppler;
k=B1_LFM/tau_LFM;%线性调频信号调制系数
ts=1/fs_LFM;%采样间隔
%生成发射信号
Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L/10))); %目标回波信号功率
A=sqrt(Prs);
Te=1143;%温度
F=4;%噪声系数
An=10*log10((1.382e-23)*Te*B_LFM*10^(F/10));%噪声强度
[y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k); 
[~,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y,tr_LFM,ts,N);
%生成回波信号
[s_echo_2,echo]=LFMtuoyinhuiboxinhao(y,R,tr_LFM,ts,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
%figure,plot(0:ts:(N-1)*ts,real(s_echo_2(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('回波信号');
%figure,plot((0:fs/length(abs(fft(abs(fftshift(echo(1,:)))))):fs-fs/length(abs(fft(abs(fftshift(echo(1,:))))))),abs((fft(echo(1,:))))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('回波信号的频谱');

%生成干扰信号
[s_ft,echo3]=LFMtuoyinganrao(R,frame_LFM,tf,Aj,tm,tau_LFM,f0_LFM,B1_LFM,N,f_doppler1,f_doppler,vf,lamta,k,ts,c,fs_LFM,temp1);

%生成噪声
[s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
%目标回波信号、假目标信号、噪声叠加在一起送入接收机
s_echo_1=s_echo_2+s_noise+s_ft;
figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('雷达接收信号');
figure,plot((0:fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_LFM-fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('雷达接收信号的频谱');
