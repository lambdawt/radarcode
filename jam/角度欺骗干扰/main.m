close all
clear 
clc
%雷达参数设定
T = 10e-6;                                                    % 脉冲持续时间
B = 15e7;                                                     % 脉冲带宽
ratio = 5;                                                    % 过采样率
Fs = ratio*B;                                                 % 采样频率
dt = 1/Fs;                                                    % 采样间隔
Nr = ceil(T/dt);                                              % 采样点数
t0 = ((0:Nr-1)-Nr/2)/Nr*T;                                    % 基本时间轴
c=3e8;                                               
                                                          %天线间距
fc=1e8;                                                       %雷达工作频率
lambda=c/fc;                                                  %工作波长
d=lambda/2;
theta_3db=0.5*lambda/d;                                       %天线波束宽度

%%放置目标
target_position=[10000,70000];                          
real_angle=atand((target_position(1))/(target_position(2)));
disp('真实角度：');disp(real_angle)
phase_deta=2*pi*d*sin(real_angle/180*pi)/lambda;
ft=Antenna_pattern(real_angle*pi/180,theta_3db);              %主瓣增益影响
%%目标回波
S1=ft*exp(1i*2*pi*fc*t0);
S2=ft*exp(1i*(2*pi*fc*t0+phase_deta));

%设置环境  信噪比20dB
S1=awgn(S1, 20);                     
S2=awgn(S2, 20);

%放置相干两点干扰源
fj=1e8;                                               %干扰源频率
J1=[60000,30000];
J2=[-20000,50000];
theta_J1=atand(J1(1)/J1(2));
theta_J2=atand(J2(1)/J2(2));
ft_J1=Antenna_pattern(theta_J1*pi/180,theta_3db);     %考虑主瓣增益
ft_J2=Antenna_pattern(theta_J2*pi/180,theta_3db);
if fj~=fc
    S11=S1+ft_J1*exp(1i*(2*pi*fj*t0+theta_J1/180*pi))+ft_J2*exp(1i*(2*pi*fj*t0+theta_J2/180*pi));
    S12=S2+ft_J2*exp(1i*(2*pi*fj*t0+theta_J2/180*pi+2*pi*fj/c*d*sin(theta_J2)))+ft_J1*exp(1i*(2*pi*fj*t0+theta_J1/180*pi+2*pi*fj/c*d*sin(theta_J1)));
else
    S11=S1+ft_J1*exp(1i*(2*pi*fj*t0+theta_J1/180*pi))+ft_J2*exp(1i*(2*pi*fj*t0+theta_J2/180*pi));;
    S12=S11*exp(1i*pi);
end
%%相位和差法测角
theta_measure=phasesd( S11 , S12 );
disp('测量角度：/°');disp(real(theta_measure));
disp('误差绝对值/°');disp(abs(theta_measure-real_angle));

%---------------------------------------------------------------------------
% theta_k=real_angle*pi/180;                   %相对等场强方向的波束倾斜角
% theta=-3.5*theta_3db:0.01:3.5*theta_3db;
% f1=exp(-1.3863*(theta-theta_k-theta_J1*pi/180).^2/theta_3db^2);
% f2=exp(-1.3863*(theta+theta_k+theta_J1*pi/180).^2/theta_3db^2);
% St1=f1+f2;
% St2=f1-f2;
% tt=St2./St1;
% figure(1)
% plot(theta,f1);hold on
% plot(theta,f2);hold off
% figure(2)
% plot(-350:350,tt);
% xlim([-150 150]);
% grid on;xlabel('theta/rad')
% figure(3)
% plot(theta,St1);title('和波束')
% figure(4)
% plot(theta,St2);title('差波束')


























% theta_k=theta_3db/3; %相对等场强方向的波束倾斜角
% theta=-2*theta_3db:0.2:2*theta_3db;
% f1=exp(-1.3863*(theta-theta_k).^2/theta_3db^2);
% f2=exp(-1.3863*(theta+theta_k).^2/theta_3db^2);
% sum1=f1+f2;
% deta1=f1-f2;
% figure(1)
% plot(real(sum1))
% figure(2)
% plot(real(deta1))









