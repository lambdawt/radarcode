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
k=BandWidth/TimeWidth;
N1=32;   %主天线阵元数 
N2=32;
auN1=1; %辅助天线个数
auN2=2;
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
t=linspace(-TimeWidth/2,TimeWidth/2,number);
signal1=exp(j*(2*pi*F0*t+k*pi*t.^2)); 
% Interfer1=exp(j*2*pi*fd1*t);
% Interfer2=exp(j*2*pi*fd2*t);

figure(101)
subplot(3,1,1)
plot(real(signal));grid on;
title('线性调频信号实部');
subplot(3,1,2)
plot(imag(signal));grid on;
title('线性调频信号虚部');
subplot(3,1,3)
plot(abs(signal));grid on;
title('线性调频信号包络');

% Interfer1=randn(1,number);
% Interfer2=randn(1,number);
%====================================================================================%
%                                   产生系统噪声信号                                 %
%====================================================================================%
Systemnoiseamp=1;
SLCSystemNoise=Systemnoiseamp*randn(auN1*auN2,number); %  mean value is 0, 辅助天线噪声
DBFSystemNoise=Systemnoiseamp*randn(N1*N2,number);  %主天线噪声
%====================================================================================%
%                                 形成通道数据                                       %
%=====================================================================================%
INR1=40;   %dB 干噪比 功率
INR2=40;   %dB
SNR=-10;    %dB信噪比
wavelength=0.1;%载波信号波长:10cm
%===================== 干扰1方向 ===============%
direction11x_interfere=0*pi/180;%干扰1方位角 
direction12x_interfere=10*pi/180;%干扰1俯仰角
%==================== 干扰2方向 ================%
direction21x_interfere=20*pi/180;%干扰2方位角
direction22x_interfere=10*pi/180;%干扰2俯仰角
%====================目标信号方向===============%Fs
direction1x_signal=10*pi/180; %期望信号方位角
direction2x_signal=10*pi/180; %期望信号俯仰角
dx=wavelength/2; %行阵元间距
dy=wavelength/2; %列阵元间距

ganraofw=[direction11x_interfere direction21x_interfere];
ganraofy=[direction12x_interfere direction22x_interfere];
figure(2)
plot(ganraofy*180/pi,ganraofw*180/pi,'ro');
hold on;grid on;
plot(direction2x_signal*180/pi,direction1x_signal*180/pi,'b*');
axis([0 60 -60 60 ]);
title('目标角度与干扰角度的分布示意');
xlabel('俯仰角/度');
ylabel('方位角/度');

theta11=sin(direction11x_interfere)*cos(direction12x_interfere);
% theta12=cos(direction11x_interfere)*cos(direction12x_interfere);
theta12=sin(direction12x_interfere);
theta21=sin(direction21x_interfere)*cos(direction22x_interfere);
% theta22=cos(direction21x_interfere)*cos(direction22x_interfere);
theta22=sin(direction22x_interfere);
theta1=sin(direction1x_signal)*cos(direction2x_signal);
% theta2=cos(direction1x_signal)*cos(direction2x_signal);
theta2=sin(direction2x_signal);


X=[0:(N1-1)]*dx;
Y=[0:(N2-1)]*dy;
X2=kron(ones(1,N2),X);
Y2=kron(Y,ones(1,N1));
X1=[3 6]*dx;
Y1=[-1 -1]*dy;
figure(1);
plot(X2,Y2,'.');axis equal;grid on;
title('天线阵');xlabel('距离（m）');ylabel('距离（m）');hold on;
plot(X1,Y1,'r*');
legend('主天线','辅助天线');

steer_interfer11=exp(j*(2*pi/wavelength*X*theta11));%主天线中的干扰1导向矢量:传播方向矢量与坐标相乘 
steer_interfer13=exp(j*(2*pi/wavelength*Y*theta12));
steer_interfer21=exp(j*(2*pi/wavelength*X*theta21));%主天线中的干扰2导向矢量
steer_interfer23=exp(j*(2*pi/wavelength*Y*theta22));
steer_signal1=exp(j*(2*pi/wavelength*X*theta1)); %主天线信号的导向矢量  
steer_signal3=exp(j*(2*pi/wavelength*Y*theta2));
steer_interfer12=kron(steer_interfer11,steer_interfer13);%先取行后取列
steer_interfer22=kron(steer_interfer21,steer_interfer23);
steer_signal2=kron(steer_signal1,steer_signal3);

auxsteer_interfer11=exp(j*(2*pi/wavelength*X1*theta11));%辅助中的干扰1导向矢量 
auxsteer_interfer13=exp(j*(2*pi/wavelength*Y1*theta12));
auxsteer_interfer21=exp(j*(2*pi/wavelength*X1*theta21));%辅助天线中的干扰2导向矢量
auxsteer_interfer23=exp(j*(2*pi/wavelength*Y1*theta22));
auxsteer_signal1=exp(j*(2*pi/wavelength*X1*theta1)); %辅助天线中的信号导向矢量
auxsteer_signal3=exp(j*(2*pi/wavelength*Y1*theta2));

aux_Gain=1;
auxsteer_interfer12=aux_Gain*auxsteer_interfer11.*auxsteer_interfer13;
auxsteer_interfer22=aux_Gain*auxsteer_interfer21.*auxsteer_interfer23;
auxsteer_signal2=aux_Gain*auxsteer_signal1.*auxsteer_signal3;
% auxsteer_interfer12=kron(auxsteer_interfer11,auxsteer_interfer13);
% auxsteer_interfer22=kron(auxsteer_interfer21,auxsteer_interfer23);
% auxsteer_signal2=kron(auxsteer_signal1,auxsteer_signal3);
%
interfer1=10^(INR1/20)*Interfer1; %幅度用20*log10(x)
interfer2=10^(INR2/20)*Interfer2; %功率用10*log10(x)
signal=10^(SNR/20)*signal;
% signal_fft=fft(signal+interfer1+interfer2);
% figure(2);
% plot(Fs/number*(0:number-1),abs(signal_fft));grid on;
% xlabel('频率(单位：Hz）'),ylabel('幅度'),title('线性调频信号频谱');
%==========================================================================%
%              interference  echo   
%==========================================================================%
input_main=steer_interfer12.'*interfer1 + steer_interfer22.'*interfer2 +...
           DBFSystemNoise; % 主天线接收到的信号(干扰1+干扰2+噪声）
input_main1=steer_interfer12.'*interfer1 + steer_interfer22.'*interfer2 +...
           steer_signal2.'*signal+DBFSystemNoise; % 主天线接收到的信号(干扰1+干扰2+信号+噪声）
input_aux=auxsteer_interfer12.'*interfer1 + auxsteer_interfer22.'*interfer2 +...
           SLCSystemNoise; % 辅助天线接收到的信号(干扰1+干扰2+噪声）   
input_aux1=auxsteer_interfer12.'*interfer1 + auxsteer_interfer22.'*interfer2 +...
           auxsteer_signal2.'*signal+SLCSystemNoise;% 辅助天线接收到的信号(干扰1+干扰2+信号+噪声）
steer_signal1=taywin(N1,-30,8)'.*exp(j*(2*pi/wavelength*X*theta1)); %主天线信号的导向矢量  
steer_signal3=taywin(N2,-30,8)'.*exp(j*(2*pi/wavelength*Y*theta2));
normal_W=kron(steer_signal1,steer_signal3);
% normal_W=steer_signal2.';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%========================计算主通道的信干噪比========================================
% S1_main=norm(signal).^2/number;
% SNR1_main=10*log10(S1_main/N1_main);
S_main=norm((normal_W.')'*(steer_signal2.'*signal)).^2/number;
I_main=norm((normal_W.')'*(steer_interfer12.'*interfer1 + steer_interfer22.'*interfer2)).^2/number;
N_main=norm((normal_W.')'*DBFSystemNoise).^2/number;
SNR_main=10*log10(S_main/N_main);
INR_main=10*log10(I_main/N_main)
%========================计算主通道的信干噪比=======================================

% %==============================================================================%
% %                       ASLC                                                  %
% %==============================================================================%
main= (normal_W.')'*input_main; %主阵列的输出
main1=(normal_W.')'*input_main1; %主天线输出

%%%%%%%%%%%%%%%主辅通道数据的量化%%%%%%%%%%%%%%%
Maximum1=[max(abs(real(main))),max(abs(imag(main))),max(abs(real(main1))),max(abs(imag(main1))),max(abs(real(input_aux(1,:)))),max(abs(real(input_aux(2,:)))),max(abs(imag(input_aux(1,:)))),max(abs(imag(input_aux(2,:)))),max(abs(real(input_aux1(1,:)))),max(abs(real(input_aux1(2,:)))),max(abs(imag(input_aux1(1,:)))),max(abs(imag(input_aux1(2,:))))];
Maximum=max(Maximum1);
main=round(main./Maximum*(2^15-1));
main1=round(main1./Maximum*(2^15-1));
input_aux=round(input_aux./Maximum*(2^15-1));
input_aux1=round(input_aux1./Maximum*(2^15-1));
%%%%%%%%%%%%%%%主辅通道数据的量化%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%主辅通道数据的的直接取整%%%%%%%%%%%%%
% main=round(main);
% main1=round(main1);
% input_aux=round(input_aux);
% input_aux1=round(input_aux1);
% %%%%%%%%%%%%%%主辅通道数据的的直接取整%%%%%%%%%%%%%

M=50;%:20:1000;%快拍数
% % NN=number/M;%每个快拍的点数
% % ii=1:length(M)
R=zeros(auN1*auN2,auN1*auN2);
r_xd=zeros(auN1*auN2,1);

Xi=main(1,1:M);
Yi=input_aux(:,1:M);
R=Yi*Yi'/M;
r_xd=Yi*Xi'/M;
% SLCW=[50+j*31;30-j*50];
SLCW=inv(R)*r_xd;%辅助天线
out_aux=SLCW'*input_aux1; %辅助天线输出
out_main=main1-out_aux;

out_aux1=SLCW'*input_aux;
out_main1=main-out_aux1;

number=length(out_main);
figure(4)
plot(Fs/number*(0:number-1),abs(fft(out_main)));grid on;hold on;
plot(Fs/number*(0:number-1),abs(fft(main1)),'r');
xlabel('频率(单位：Hz）'),ylabel('幅度')
title('相消前后信号频谱');
legend('相消后','相消前');
figure(105)
subplot(2,1,1),plot(real((main1)));grid on;hold on;xlabel('点数');ylabel('幅度');title('相消前输出信号的实部');
subplot(2,1,2),plot(real((out_main)));xlabel('点数'),ylabel('幅度');grid on;title('相消后输出信号的实部');


%===============================================================================%
%                    波束扫描  Forming beam pattern                             %
%===============================================================================%
k1=1;
for thta=-60:1:60
    k2=1;
    for phi=0:0.5:60
        th_a1=exp(j*2*pi/wavelength*X*sin(thta*pi/180)*cos(phi*pi/180));
%         th_p1=exp(j*2*pi/wavelength*Y*cos(thta*pi/180)*cos(phi*pi/180));
        th_p1=exp(j*2*pi/wavelength*Y*sin(phi*pi/180));
        th_all1=kron(th_a1,th_p1).';
        th_a2=exp(j*2*pi/wavelength*X1*sin(thta*pi/180)*cos(phi*pi/180));
%         th_p2=exp(j*2*pi/wavelength*Y1*cos(thta*pi/180)*cos(phi*pi/180));
        th_p2=exp(j*2*pi/wavelength*Y1*sin(phi*pi/180));
        th_all2=aux_Gain*(th_a2.*th_p2).';
%         th_all2=kron(th_a2,th_p2).';
%         yy440_cancel(k1,k2)=abs((normal_W')*th_all1-SLCW'*th_all2);
        yy440_qian(k1,k2)=abs(((normal_W.')')*th_all1);
        yy440_fuzhu(k1,k2)=abs(SLCW'*th_all2);
        yy440_cancel(k1,k2)=abs((((normal_W.')')*th_all1)-(SLCW'*th_all2));
        k2=k2+1;
    end
    k1=k1+1;
end
% y_cancel=20*log10((yy440_cancel)./max(max(yy440_cancel)));
% y_qian=20*log10((yy440_qian)./max(max(yy440_qian)));
% y_fuzhu=20*log10((yy440_fuzhu)./max(max(yy440_fuzhu)));
y_cancel=20*log10(yy440_cancel);
y_qian=20*log10(yy440_qian);
y_fuzhu=20*log10(yy440_fuzhu);
p1=0:0.5:60;
a1=-60:1:60;
figure(5)
mesh(p1,a1,y_qian);
title('对消前波束图');
xlabel('俯仰角/度');ylabel('方位角/度');zlabel('幅度/dB');
% axis([0 60 -60 60 -40 40]);
figure(100)
mesh(p1,a1,y_fuzhu);
title('辅助天线方向图');
xlabel('俯仰角/度');ylabel('方位角/度');zlabel('幅度/dB');
figure(6)
plot(p1,y_qian(1+direction1x_signal*180/pi+60,:));
xlabel('俯仰角/度');ylabel('幅度/dB');
grid on;hold on;
figure(7)
plot(a1,y_qian(:,1+direction2x_signal*180/pi/0.5));
xlabel('方位角/度');ylabel('幅度/dB');
grid on;hold on;


%% 论文绘图使用
hfig=figure(8);
mesh(p1,a1,y_cancel);
view([90,0]);
title('\fontname{宋体}自适应旁瓣对消后天线方向图', 'FontSize',10.5);
xlabel('\fontname{宋体}俯仰角\fontname{Times new roman}/°');
ylabel('\fontname{宋体}方位角\fontname{Times new roman}/°');
zlabel('\fontname{宋体}幅度\fontname{Times new roman}/dB');
set(gca,'FontName','Times New Roman','FontSize',10.5);
grid on;
axis tight;
% 关于图片大小，应该通过计算得到宽度和高度。
% 例如A4纸的宽度为21cm，如果左右页边距各为2.5cm，正文的宽度就是16cm。
% 如果页面是单栏，设为宽度14cm，高度8.6cm比较合适（黄金比例0.618）。
% 双栏的话，宽7cm，高4.3cm即可。
% figWidth = 12;  % 设置图片宽度
% figHeight = 7.416;  % 设置图片高度
% set(hfig,'PaperUnits','centimeters'); % 图片尺寸所用单位，inches—英寸，centimeters—厘米
% set(hfig,'PaperPosition',[0 0 figWidth figHeight]);
% fileout = 'SLC处理后天线方向图.'; % 输出图片的文件名
% print(hfig,[fileout,'tif'],'-r600','-dtiff'); % 设置图片格式、分辨率
%%


% axis([0 60 -60 60 -40 40]);
figure(9)
plot(p1,y_cancel(1+direction1x_signal*180/pi+60,:),'r');
xlabel('俯仰角/度');ylabel('幅度/dB');
title('对消后主通道俯仰维波束')
grid on
figure(10)
plot(a1,y_cancel(:,1+direction2x_signal*180/pi/0.5),'r');
xlabel('方位角/度');ylabel('幅度/dB');
title('对消后主通道方位维波束')
grid on
% =====================================================================%
figure(11);
subplot(2,1,1);plot(real(main));grid on;zoom on; %对消前的输出
ylabel('幅度');title('对消前输入');
subplot(2,1,2);plot(real(out_main));grid on;zoom on; %对消后的输出
 xlabel('采样点数'); ylabel('幅度');title('对消后输出');
%  subplot(3,1,3);plot(abs(SMI_OUT));grid on;zoom on; %SMI输出
% ylabel('幅度');title('SMI输出');





