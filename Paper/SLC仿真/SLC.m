%===========================================================================================%
%                       自适应旁瓣相消(线阵）   ----ASLC                                    %
%===========================================================================================%
%===================================================================================%
%                     仿真参数初始化 (线性调频信号＋点频干扰)                       %
%===================================================================================%
close all;clear all;clc;
C=3.0e8; %光速(m/s)
Fs=40.0e6;%采样频率
F0=30e6;   %中心频率
TimeWidth=25e-6; %时宽   
BandWidth=5e6; %带宽
number=fix(Fs*TimeWidth); %数据采样个数 大小影响计算精度
N=3;   %主天线阵元数  
auN=2; %辅助天线个数
%=================================加入干扰信息==========================================%
fd1=31e6;    % 干扰1频率(HZ)
fd2=29e6;
if rem(number,2)~=0  % if number is not pow of 2,then number+1 
   number=number+1;
end    
Interfer1=zeros(1,number);
Interfer2=zeros(1,number);
signal = zeros(1,number);
for i=-fix(number/2):fix(number/2)-1
 Interfer1(i+fix(number/2)+1)=exp(j*2*pi*fd1*i/Fs);  % 干扰1的采样值(模值为1)
 Interfer2(i+fix(number/2)+1)=exp(j*2*pi*fd2*i/Fs);  % 干扰2的采样值 
 signal(i+fix(number/2)+1)=exp(j*(2*pi*F0*i/Fs+pi*(BandWidth/TimeWidth)*(i/Fs)^2));  % 信号的采样值
end

figure(1);
subplot(2,1,1);
plot(real(signal));grid on;zoom on;
xlabel('采样点数'),ylabel('幅度'),title('实部');
subplot(2,1,2);
plot(imag(signal));grid on;zoom on;
xlabel('采样点数'),ylabel('幅度'),title('虚部');
% signal_fft=fft(signal+Interfer1+Interfer2);
% figure(2);
% plot(Fs/number*(0:number-1),abs(signal_fft));grid on;
% xlabel('频率(单位：Hz）'),ylabel('幅度'),title('线性调频信号频谱');
%====================================================================================%
%                                   产生系统噪声信号                                 %
%====================================================================================%
Systemnoiseamp=1;
SLCSystemNoise=Systemnoiseamp*randn(auN,number); %  mean value is 0, 辅助天线噪声
DBFSystemNoise=Systemnoiseamp*randn(N,number);  %主天线噪声
%====================================================================================%
%                                 形成通道数据                                       %
%=====================================================================================%
INR1=20;   %dB 干噪比 功率
INR2=20;   %dB
SNR=10;    %dB信噪比
% f=200e6; %载波中心频率
wavelength=0.1;%载波信号波长:10cm
%===================== 干扰1方向 ===============%
direction1x_interfere=41*pi/180; 
%==================== 干扰2方向 ================%
direction2x_interfere=-30*pi/180;
%====================目标信号方向===============%
directionx_signal=10*pi/180; %方位角
dx=wavelength/2; %阵元间距
dass=[-1*dx 0*dx]; %辅助天线的坐标    表示
% d=[1*dx 2*dx 3*dx 4*dx 5*dx 6*dx 7*dx 8*dx 9*dx 10*dx 11*dx 12*dx];% 13*dx 14*dx...
  %15*dx 16*dx 17*dx 18*dx 19*dx 20*dx 21*dx];% 主天线的坐标表示
d=dx:dx:N*dx;
k1=sin(direction1x_interfere);%干扰1
steer_interfer1=exp(j*(2*pi/wavelength*k1*d));%主天线中的干扰1导向矢量:传播方向矢量与坐标相乘 
auxsteer_interfer1=exp(j*(2*pi/wavelength*k1*dass));%辅助中的干扰1导向矢量  
k2=sin(direction2x_interfere);%干扰2的传播方向矢量
steer_interfer2=exp(j*(2*pi/wavelength*k2*d));%主天线中的干扰2导向矢量  
auxsteer_interfer2=exp(j*(2*pi/wavelength*k2*dass));%辅助天线中的干扰2导向矢量
k3=sin(directionx_signal);%信号的传播方向矢量
steer_signal=exp(j*(2*pi/wavelength*k3*d)); %主天线信号的导向矢量
auxsteer_signal=exp(j*(2*pi/wavelength*k3*dass)); %辅助天线中的信号导向矢量
%
interfer1=10^(INR1/10)*Interfer1; %幅度用20*log10(x)
interfer2=10^(INR2/10)*Interfer2; %功率用10*log10(x)
signal=10^(SNR/10)*signal;
% signal_fft=fft(signal+interfer1+interfer2);
% figure(2);
% plot(Fs/number*(0:number-1),abs(signal_fft));grid on;
% xlabel('频率(单位：Hz）'),ylabel('幅度'),title('线性调频信号频谱');
%==========================================================================%
%              interference  echo   
%==========================================================================%
input_main=steer_interfer1.'*interfer1 + steer_interfer2.'*interfer2 +...
           DBFSystemNoise; % 主天线接收到的信号(干扰1+干扰2+噪声）
input_aux=auxsteer_interfer1.'*interfer1 + auxsteer_interfer2.'*interfer2 +...
         SLCSystemNoise; % 辅助天线接收到的信号     
input_main1=steer_interfer1.'*interfer1 + steer_interfer2.'*interfer2 +...
           steer_signal.'*signal+DBFSystemNoise; % 主天线接收到的信号(干扰1+干扰2+信号+噪声）
input_aux1=auxsteer_interfer1.'*interfer1 + auxsteer_interfer2.'*interfer2 +...
           auxsteer_signal.'*signal+SLCSystemNoise;%     
%========================SMI===============================================
Rxx=input_main*input_main'/number;
SMI_W=inv(Rxx)*steer_signal.';%权矢量（线性约束最小方差(LCMV)准则） 
SMI_OUT=SMI_W'*input_main;
normal_W=(taywin(N,-20,8)'.*steer_signal).';
% normal_W=steer_signal.';
input_main_new=steer_signal.'*signal;
signal_new1=(normal_W)'*input_main_new;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %==============================================================================%
% %                       ASLC                                                  %
% %==============================================================================%
% SLL=-25;
% m=3;
% tay=taywin(N,SLL,m);
% main= (steer_signal.'.*tay)'*input_main; %主阵列的输出
main= (normal_W)'*input_main; %主阵列的输出
[size_x,size_y]=size(input_aux); %返回辅助天线的个数和采样点数
r_xd=zeros(1,size_y); R=zeros(size_x,size_x);
% R=input_aux*input_aux'/number; %辅助天线自相关矩阵
R=zeros(auN,auN);
r_xd=zeros(auN,1);
SLCW=zeros(auN,1);
M=50;%:20:1000;%快拍数
% NN=number/M;%每个快拍的点数
% ii=1:length(M)
for i=1:M;
    R=input_aux(:,i)*input_aux(:,i)'+R;%%%%%%%%%%是几个辅助通道的值一块计算其自相关函数不是辅助天线各自的值计算吗？？？
    r_xd=input_aux(:,i)*main(:,i)'+r_xd;
end
    R=R/M;
    r_xd=r_xd/M;
% r_xd=input_aux*main'/number; %主天线与辅助天线的互相关向量
SLCW=inv(R)*r_xd; %辅助天线
% r_xd=input_aux*main'/number; %主天线与辅助天线的互相关向量
% SLCW=inv(R)*r_xd; %计算旁瓣相消最佳权(这里的期望信号d就是主阵的输出) 
out_aux=SLCW'*input_aux1; %辅助天线输出
main=(normal_W)'*input_main1; %主天线输出   
% main2=(steer_signal.')'*(steer_signal.'*signal+DBFSystemNoise);
out_main=main-out_aux;  %对消后输出

out_aux1=SLCW'*input_aux;
main1=(normal_W)'*input_main;
out_main1=main1-out_aux1;

p_a=norm(main1).^2/number;%相消前平均功率
Es0=10*log10(p_a);
p_b=norm(out_main1).^2/number;  %相消后平均功率
Es1=10*log10(p_b);
CR1=10*log10(p_a/p_b)%对消比

% figure(100)
% plot(real(main1));grid on;hold on;
% plot(real(aux),'k');hold on;
% plot(real(out_main),'r');
% legend('对消前主天线','辅助天线','对消后主天线');

figure(2);
plot(Fs/number*(0:number-1),abs(fft(main)),'r');grid on;hold on;
plot(Fs/number*(0:number-1),(abs(fft(out_main)))-850*ones(1,length(out_main)));grid on;
xlabel('频率（Hz）');
ylabel('幅度');
title('相消前后信号频谱');
legend('相消前','相消后');


%========================计算主通道的信干噪比========================================
S_main=norm((normal_W)'*(steer_signal.'*signal)).^2/number;
I_main=norm((normal_W)'*(steer_interfer1.'*interfer1 + steer_interfer2.'*interfer2)).^2/number;
N_main=norm((normal_W)'*DBFSystemNoise).^2/number;
SNR_main=10*log10(S_main/N_main);
INR_main=10*log10(I_main/N_main);
%========================计算主通道的信干噪比=======================================

%===============================================================================%
%                    波束扫描  Forming beam pattern                             %
%===============================================================================%


row=0; %行
for loop_v=-90*pi/180:0.2*pi/180:90*pi/180       
     row=row+1
     k= sin(loop_v); %传播方向矢量
     steer_inter=exp(j*(2*pi/wavelength*k*d)).'; %扫描导向矢量,searching variable
%      SMI_W=SMI_W.*hamming(length(SMI_W));
     P_SMI(row)=abs((SMI_W)'*steer_inter); % SMI DBF(方向图:见阵列信号处理3)
     P_normal(row)=((normal_W)'*steer_inter);
%    P_normal(row)=abs((normal_W.*tay)'*steer_inter);  % Usual DBF 
end
u=[-90*pi/180:0.2*pi/180:90*pi/180]*180/pi;
%%==========================SMI主天线方向图====================================%
%%%%%%%%%%%%%%%%%%%%%%%%%%%SMI主天线方向图是什么？？？%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3);
P_SMI_1= P_SMI./(max(P_SMI)); %归一化
P_SMI_1=20*log10(P_SMI_1); %dB
plot(u,P_SMI_1);grid on;zoom on; %画归一化后的SMI向图：指向信号
axis([-100 100 -60 0]);
xlabel('角度');ylabel('增益(dB)');
title('SMI方向图');
%%%%%%%%%%%%%%%%%%%%%%%%%%%SMI主天线方向图是什么？？？%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%========================对消前的方向图=======================================%
figure(4);
% P_normal_1= P_normal./(max( P_normal));%归一化
P_normal_1=abs( P_normal);
% P_normal_1=20*log10( P_normal_1./(max( P_normal_1))); %dB
P_normal_1=20*log10( P_normal_1); %dB
plot(u, P_normal_1);grid on;zoom on; %画对消前的方向图
% axis([-100 100 -60 0]); %限定坐标
xlabel('角度');ylabel('增益(dB)');
title('对消前方向图');
% legend('3度方向图');
% %============================================================================%
% %                     Figure  SLC                                            %
% %============================================================================%
figure(5);
subplot(2,2,1);plot(real(main));grid on;zoom on; %对消前的输出
xlabel('采样点数');ylabel('幅度');title('相消前输出信号的实部');
subplot(2,2,3);plot(real(out_main));grid on;zoom on; %对消后的输出
xlabel('采样点数');ylabel('幅度');title('相消后输出信号的实部');
subplot(2,2,2);plot(imag(main));grid on;zoom on; %对消前的输出
xlabel('采样点数');ylabel('幅度');title('相消前输出信号的虚部');
subplot(2,2,4);plot(imag(out_main));grid on;zoom on; %对消后的输出
xlabel('采样点数');ylabel('幅度');title('相消后输出信号的虚部');
 
figure(8);
subplot(2,1,1);plot(20*log10(abs(main)));grid on;zoom on; %对消前的输出
xlabel('采样点数');ylabel('幅度(dB)');title('相消前输出信号的幅度');
subplot(2,1,2);plot(20*log10(abs(out_main)));grid on;zoom on; %对消后的输出
 xlabel('采样点数'); ylabel('幅度(dB)');title('相消后输出信号的幅度');

%  subplot(3,1,3);plot(abs(SMI_OUT));grid on;zoom on; %SMI输出
% ylabel('幅度');title('SMI输出');
%==============================对消比============================================%
% CR=max(abs(main))./max(abs(out_main));
% CR=20*log10(CR)
p_a=norm(main).^2/number;%相消前平均功率
Es0=10*log10(p_a)
p_b=norm(out_main).^2/number;  %相消后平均功率
Es1=10*log10(p_b)
CR=10*log10(p_a/p_b)%对消比
%%=============================辅助天线方向图====================================%
row=0; %行
for loop_v=-90*pi/180:0.2*pi/180:90*pi/180       
     row=row+1;
     k= sin(loop_v); %传播方向矢量
    steer_inter=exp(j*(2*pi/wavelength*k*dass)).'; %扫描导向矢量,searching variable
     P_aux(row)=(SLCW'*steer_inter);%辅助天线方向图
  end
%   P_aux_1=P_aux./max(P_normal);
  P_aux_1=abs(P_aux);
  P_aux_1=20*log10(P_aux_1);
figure(6);
plot(u,P_aux_1,'r');grid on;zoom on;hold on;
plot(u,P_normal_1);
% axis([-100 100 -60 0]);
xlabel('角度');ylabel('增益(dB)');
legend('对消前方向图','辅助天线方向图')
title('辅助天线方向图');
%%=============================相减后的方向图====================================%
P =abs(P_normal-P_aux);
% P_1=P./(max(P_normal)); % 归一化
% P_1=P./(max(P)); % 归一化
% P_2=P;
P_1=20*log10(P); %dB 幅度用20*log10(x)
% P_2=20*log10(P_2);
figure(7);
plot(u,P_1,'r');grid on;zoom on;hold on;
% plot(u,P_aux_1,'g');grid on;zoom on;hold on;
% plot(u,P_normal_1);
xlabel('角度');ylabel('增益(dB)');
title('对消后方向图');
% text(30,-3,'干扰方向为35和-25度');
% legend('相消后','辅助天线','相消前');
%=================================END=============================================%
%  abc=P_1(35)-P_normal_1(35);
h=conj(fliplr(signal)); 
% %----对匹配滤波器系数进行抽取
% d1=fix(number/a);
% for m=1:d1
%     hh(m)=h(a*(m-1)+1);%----对滤波器进行抽取
% end
% h1=hh./sqrt(sum(abs(hh).^2)); %匹配滤波器系数归一化
h1=h./sqrt(sum(abs(h).^2));%匹配滤波器系数归一化
Match_out=conv(h1,out_main);
figure(14)
plot(20*log10(abs(Match_out)/max(abs(Match_out))));grid on;
xlabel('点数');
ylabel('幅度');
title('匹配滤波后结果');

t=linspace(-TimeWidth/2,TimeWidth/2,length(h1));
K1=0.08;n=2;
W=K1+(1-K1)*(cos(pi*t/TimeWidth)).^n;%----海明窗
% W=0.42+0.5*cos(2*pi*t/TimeWidth)+0.08*cos(4*pi*t/TimeWidth);%----Blackman窗
Wt=W.*h1;
Wt=Wt./sqrt(sum(abs(Wt).^2));
Window_out=conv(Wt,out_main);
figure(15)
plot(20*log10(abs(Window_out)/max(abs(Window_out))));grid on;
xlabel('点数');
ylabel('幅度');
title('加Hamming窗后结果');

RMS=rms1(Window_out);
q=20*log10(max(abs(Match_out))/max(abs(Window_out)));