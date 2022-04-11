%===========================================================================================%
%                       自适应旁瓣相消(线阵）   ----ASLC                                    %
close all;clear all;clc;
N1=30;   %主天线阵元数
N2=30;
wavelength=0.1;%载波信号波长:10cm

%====================目标信号方向===============%Fs
direction1x_signal=0*pi/180; %期望信号方位角
direction2x_signal=20*pi/180; %期望信号俯仰角
dx=wavelength/2; %行阵元间距
dy=wavelength/2; %列阵元间距

theta1=sin(direction1x_signal)*cos(direction2x_signal);
theta2=sin(direction2x_signal);

X=[0:(N1-1)]*dx;
Y=[0:(N2-1)]*dy;
X2=kron(ones(1,N2),X);
Y2=kron(Y,ones(1,N1));

% figure
% plot(X2,Y2,'.');axis equal;grid on;
% title('天线阵');xlabel('距离（m）');ylabel('距离（m）');

steer_signal1=exp(j*(2*pi/wavelength*X*theta1)); %主天线信号的导向矢量  
steer_signal3=exp(j*(2*pi/wavelength*Y*theta2));
steer_signal2=kron(steer_signal1,steer_signal3);
normal_W=steer_signal2;

%===============================================================================%
%                    波束扫描  Forming beam pattern                             %
%===============================================================================%
k1=1;
for thta=-90:0.5:90
    k2=1;
    for phi=0:0.3:72
        th_a1=exp(j*2*pi/wavelength*X*sin(thta*pi/180)*cos(phi*pi/180));
        th_p1=exp(j*2*pi/wavelength*Y*sin(phi*pi/180));
        th_all1=kron(th_a1,th_p1).';
        yy440_qian(k1,k2)=abs(((normal_W.')')*th_all1);
        k2=k2+1;
    end
    k1=k1+1;
end
y_qian=20*log10(yy440_qian/max(max(yy440_qian)));
p1=-90:0.5:90;
a1=0:0.3:72;
figure
mesh(a1,p1,y_qian);
title('矩形面阵方向图');
xlabel('俯仰角');ylabel('方位角');zlabel('幅度(dB)');
axis([0 72 -100 100 -150 20]);

% figure
% plot(a1,y_qian(1+direction1x_signal*180/pi+90,:));
% xlabel('俯仰角/度');ylabel('幅度/dB');
% grid on;hold on;
% 
% figure
% plot(p1,y_qian(:,1+direction2x_signal*180/pi/0.5));
% xlabel('方位角/度');ylabel('幅度/dB');
% grid on;hold on;