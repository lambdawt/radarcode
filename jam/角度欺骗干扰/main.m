close all
clear 
clc
%�״�����趨
T = 10e-6;                                                    % �������ʱ��
B = 15e7;                                                     % �������
ratio = 5;                                                    % ��������
Fs = ratio*B;                                                 % ����Ƶ��
dt = 1/Fs;                                                    % �������
Nr = ceil(T/dt);                                              % ��������
t0 = ((0:Nr-1)-Nr/2)/Nr*T;                                    % ����ʱ����
c=3e8;                                               
                                                          %���߼��
fc=1e8;                                                       %�״﹤��Ƶ��
lambda=c/fc;                                                  %��������
d=lambda/2;
theta_3db=0.5*lambda/d;                                       %���߲������

%%����Ŀ��
target_position=[10000,70000];                          
real_angle=atand((target_position(1))/(target_position(2)));
disp('��ʵ�Ƕȣ�');disp(real_angle)
phase_deta=2*pi*d*sin(real_angle/180*pi)/lambda;
ft=Antenna_pattern(real_angle*pi/180,theta_3db);              %��������Ӱ��
%%Ŀ��ز�
S1=ft*exp(1i*2*pi*fc*t0);
S2=ft*exp(1i*(2*pi*fc*t0+phase_deta));

%���û���  �����20dB
S1=awgn(S1, 20);                     
S2=awgn(S2, 20);

%��������������Դ
fj=1e8;                                               %����ԴƵ��
J1=[60000,30000];
J2=[-20000,50000];
theta_J1=atand(J1(1)/J1(2));
theta_J2=atand(J2(1)/J2(2));
ft_J1=Antenna_pattern(theta_J1*pi/180,theta_3db);     %������������
ft_J2=Antenna_pattern(theta_J2*pi/180,theta_3db);
if fj~=fc
    S11=S1+ft_J1*exp(1i*(2*pi*fj*t0+theta_J1/180*pi))+ft_J2*exp(1i*(2*pi*fj*t0+theta_J2/180*pi));
    S12=S2+ft_J2*exp(1i*(2*pi*fj*t0+theta_J2/180*pi+2*pi*fj/c*d*sin(theta_J2)))+ft_J1*exp(1i*(2*pi*fj*t0+theta_J1/180*pi+2*pi*fj/c*d*sin(theta_J1)));
else
    S11=S1+ft_J1*exp(1i*(2*pi*fj*t0+theta_J1/180*pi))+ft_J2*exp(1i*(2*pi*fj*t0+theta_J2/180*pi));;
    S12=S11*exp(1i*pi);
end
%%��λ�Ͳ���
theta_measure=phasesd( S11 , S12 );
disp('�����Ƕȣ�/��');disp(real(theta_measure));
disp('������ֵ/��');disp(abs(theta_measure-real_angle));

%---------------------------------------------------------------------------
% theta_k=real_angle*pi/180;                   %��Եȳ�ǿ����Ĳ�����б��
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
% plot(theta,St1);title('�Ͳ���')
% figure(4)
% plot(theta,St2);title('���')


























% theta_k=theta_3db/3; %��Եȳ�ǿ����Ĳ�����б��
% theta=-2*theta_3db:0.2:2*theta_3db;
% f1=exp(-1.3863*(theta-theta_k).^2/theta_3db^2);
% f2=exp(-1.3863*(theta+theta_k).^2/theta_3db^2);
% sum1=f1+f2;
% deta1=f1-f2;
% figure(1)
% plot(real(sum1))
% figure(2)
% plot(real(deta1))









