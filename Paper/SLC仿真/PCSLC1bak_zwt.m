%===========================================================================================%
%                       ����Ӧ�԰�����(����   ----ASLC                                    %
%===========================================================================================%
%===================================================================================%
%                     ���������ʼ�� (���Ե�Ƶ�źţ���Ƶ����)                       %
%===================================================================================%
close all;clear all;clc;
C=3.0e8; %����(m/s)
Fs=40.0e6;%����Ƶ��
F0=30e6;   %����Ƶ��
TimeWidth=25e-6; %ʱ�� 
BandWidth=5e6; %����
number=fix(Fs*TimeWidth); %���ݲ������� ��СӰ����㾫��
k=BandWidth/TimeWidth;
N1=32;   %��������Ԫ�� 
N2=32;
auN1=1; %�������߸���
auN2=2;
%=================================���������Ϣ==========================================%
fd1=31e6;    % ����1Ƶ��(HZ)
fd2=29e6;
if rem(number,2)~=0  % if number is not pow of 2,then number+1 
   number=number+1;
end    
Interfer1=zeros(1,number);
Interfer2=zeros(1,number);
signal = zeros(1,number);
for i=-fix(number/2):fix(number/2)-1
    Interfer1(i+fix(number/2)+1)=exp(j*2*pi*fd1*i/Fs);  % ����1�Ĳ���ֵ(ģֵΪ1)
    Interfer2(i+fix(number/2)+1)=exp(j*2*pi*fd2*i/Fs);  % ����2�Ĳ���ֵ 
    signal(i+fix(number/2)+1)=exp(j*(2*pi*F0*i/Fs+pi*(BandWidth/TimeWidth)*(i/Fs)^2));  % �źŵĲ���ֵ
end
t=linspace(-TimeWidth/2,TimeWidth/2,number);
signal1=exp(j*(2*pi*F0*t+k*pi*t.^2)); 
% Interfer1=exp(j*2*pi*fd1*t);
% Interfer2=exp(j*2*pi*fd2*t);

figure(101)
subplot(3,1,1)
plot(real(signal));grid on;
title('���Ե�Ƶ�ź�ʵ��');
subplot(3,1,2)
plot(imag(signal));grid on;
title('���Ե�Ƶ�ź��鲿');
subplot(3,1,3)
plot(abs(signal));grid on;
title('���Ե�Ƶ�źŰ���');

% Interfer1=randn(1,number);
% Interfer2=randn(1,number);
%====================================================================================%
%                                   ����ϵͳ�����ź�                                 %
%====================================================================================%
Systemnoiseamp=1;
SLCSystemNoise=Systemnoiseamp*randn(auN1*auN2,number); %  mean value is 0, ������������
DBFSystemNoise=Systemnoiseamp*randn(N1*N2,number);  %����������
%====================================================================================%
%                                 �γ�ͨ������                                       %
%=====================================================================================%
INR1=40;   %dB ����� ����
INR2=40;   %dB
SNR=-10;    %dB�����
wavelength=0.1;%�ز��źŲ���:10cm
%===================== ����1���� ===============%
direction11x_interfere=0*pi/180;%����1��λ�� 
direction12x_interfere=10*pi/180;%����1������
%==================== ����2���� ================%
direction21x_interfere=20*pi/180;%����2��λ��
direction22x_interfere=10*pi/180;%����2������
%====================Ŀ���źŷ���===============%Fs
direction1x_signal=10*pi/180; %�����źŷ�λ��
direction2x_signal=10*pi/180; %�����źŸ�����
dx=wavelength/2; %����Ԫ���
dy=wavelength/2; %����Ԫ���

ganraofw=[direction11x_interfere direction21x_interfere];
ganraofy=[direction12x_interfere direction22x_interfere];
figure(2)
plot(ganraofy*180/pi,ganraofw*180/pi,'ro');
hold on;grid on;
plot(direction2x_signal*180/pi,direction1x_signal*180/pi,'b*');
axis([0 60 -60 60 ]);
title('Ŀ��Ƕ�����ŽǶȵķֲ�ʾ��');
xlabel('������/��');
ylabel('��λ��/��');

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
title('������');xlabel('���루m��');ylabel('���루m��');hold on;
plot(X1,Y1,'r*');
legend('������','��������');

steer_interfer11=exp(j*(2*pi/wavelength*X*theta11));%�������еĸ���1����ʸ��:��������ʸ����������� 
steer_interfer13=exp(j*(2*pi/wavelength*Y*theta12));
steer_interfer21=exp(j*(2*pi/wavelength*X*theta21));%�������еĸ���2����ʸ��
steer_interfer23=exp(j*(2*pi/wavelength*Y*theta22));
steer_signal1=exp(j*(2*pi/wavelength*X*theta1)); %�������źŵĵ���ʸ��  
steer_signal3=exp(j*(2*pi/wavelength*Y*theta2));
steer_interfer12=kron(steer_interfer11,steer_interfer13);%��ȡ�к�ȡ��
steer_interfer22=kron(steer_interfer21,steer_interfer23);
steer_signal2=kron(steer_signal1,steer_signal3);

auxsteer_interfer11=exp(j*(2*pi/wavelength*X1*theta11));%�����еĸ���1����ʸ�� 
auxsteer_interfer13=exp(j*(2*pi/wavelength*Y1*theta12));
auxsteer_interfer21=exp(j*(2*pi/wavelength*X1*theta21));%���������еĸ���2����ʸ��
auxsteer_interfer23=exp(j*(2*pi/wavelength*Y1*theta22));
auxsteer_signal1=exp(j*(2*pi/wavelength*X1*theta1)); %���������е��źŵ���ʸ��
auxsteer_signal3=exp(j*(2*pi/wavelength*Y1*theta2));

aux_Gain=1;
auxsteer_interfer12=aux_Gain*auxsteer_interfer11.*auxsteer_interfer13;
auxsteer_interfer22=aux_Gain*auxsteer_interfer21.*auxsteer_interfer23;
auxsteer_signal2=aux_Gain*auxsteer_signal1.*auxsteer_signal3;
% auxsteer_interfer12=kron(auxsteer_interfer11,auxsteer_interfer13);
% auxsteer_interfer22=kron(auxsteer_interfer21,auxsteer_interfer23);
% auxsteer_signal2=kron(auxsteer_signal1,auxsteer_signal3);
%
interfer1=10^(INR1/20)*Interfer1; %������20*log10(x)
interfer2=10^(INR2/20)*Interfer2; %������10*log10(x)
signal=10^(SNR/20)*signal;
% signal_fft=fft(signal+interfer1+interfer2);
% figure(2);
% plot(Fs/number*(0:number-1),abs(signal_fft));grid on;
% xlabel('Ƶ��(��λ��Hz��'),ylabel('����'),title('���Ե�Ƶ�ź�Ƶ��');
%==========================================================================%
%              interference  echo   
%==========================================================================%
input_main=steer_interfer12.'*interfer1 + steer_interfer22.'*interfer2 +...
           DBFSystemNoise; % �����߽��յ����ź�(����1+����2+������
input_main1=steer_interfer12.'*interfer1 + steer_interfer22.'*interfer2 +...
           steer_signal2.'*signal+DBFSystemNoise; % �����߽��յ����ź�(����1+����2+�ź�+������
input_aux=auxsteer_interfer12.'*interfer1 + auxsteer_interfer22.'*interfer2 +...
           SLCSystemNoise; % �������߽��յ����ź�(����1+����2+������   
input_aux1=auxsteer_interfer12.'*interfer1 + auxsteer_interfer22.'*interfer2 +...
           auxsteer_signal2.'*signal+SLCSystemNoise;% �������߽��յ����ź�(����1+����2+�ź�+������
steer_signal1=taywin(N1,-30,8)'.*exp(j*(2*pi/wavelength*X*theta1)); %�������źŵĵ���ʸ��  
steer_signal3=taywin(N2,-30,8)'.*exp(j*(2*pi/wavelength*Y*theta2));
normal_W=kron(steer_signal1,steer_signal3);
% normal_W=steer_signal2.';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%========================������ͨ�����Ÿ����========================================
% S1_main=norm(signal).^2/number;
% SNR1_main=10*log10(S1_main/N1_main);
S_main=norm((normal_W.')'*(steer_signal2.'*signal)).^2/number;
I_main=norm((normal_W.')'*(steer_interfer12.'*interfer1 + steer_interfer22.'*interfer2)).^2/number;
N_main=norm((normal_W.')'*DBFSystemNoise).^2/number;
SNR_main=10*log10(S_main/N_main);
INR_main=10*log10(I_main/N_main)
%========================������ͨ�����Ÿ����=======================================

% %==============================================================================%
% %                       ASLC                                                  %
% %==============================================================================%
main= (normal_W.')'*input_main; %�����е����
main1=(normal_W.')'*input_main1; %���������

%%%%%%%%%%%%%%%����ͨ�����ݵ�����%%%%%%%%%%%%%%%
Maximum1=[max(abs(real(main))),max(abs(imag(main))),max(abs(real(main1))),max(abs(imag(main1))),max(abs(real(input_aux(1,:)))),max(abs(real(input_aux(2,:)))),max(abs(imag(input_aux(1,:)))),max(abs(imag(input_aux(2,:)))),max(abs(real(input_aux1(1,:)))),max(abs(real(input_aux1(2,:)))),max(abs(imag(input_aux1(1,:)))),max(abs(imag(input_aux1(2,:))))];
Maximum=max(Maximum1);
main=round(main./Maximum*(2^15-1));
main1=round(main1./Maximum*(2^15-1));
input_aux=round(input_aux./Maximum*(2^15-1));
input_aux1=round(input_aux1./Maximum*(2^15-1));
%%%%%%%%%%%%%%%����ͨ�����ݵ�����%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%����ͨ�����ݵĵ�ֱ��ȡ��%%%%%%%%%%%%%
% main=round(main);
% main1=round(main1);
% input_aux=round(input_aux);
% input_aux1=round(input_aux1);
% %%%%%%%%%%%%%%����ͨ�����ݵĵ�ֱ��ȡ��%%%%%%%%%%%%%

M=50;%:20:1000;%������
% % NN=number/M;%ÿ�����ĵĵ���
% % ii=1:length(M)
R=zeros(auN1*auN2,auN1*auN2);
r_xd=zeros(auN1*auN2,1);

Xi=main(1,1:M);
Yi=input_aux(:,1:M);
R=Yi*Yi'/M;
r_xd=Yi*Xi'/M;
% SLCW=[50+j*31;30-j*50];
SLCW=inv(R)*r_xd;%��������
out_aux=SLCW'*input_aux1; %�����������
out_main=main1-out_aux;

out_aux1=SLCW'*input_aux;
out_main1=main-out_aux1;

number=length(out_main);
figure(4)
plot(Fs/number*(0:number-1),abs(fft(out_main)));grid on;hold on;
plot(Fs/number*(0:number-1),abs(fft(main1)),'r');
xlabel('Ƶ��(��λ��Hz��'),ylabel('����')
title('����ǰ���ź�Ƶ��');
legend('������','����ǰ');
figure(105)
subplot(2,1,1),plot(real((main1)));grid on;hold on;xlabel('����');ylabel('����');title('����ǰ����źŵ�ʵ��');
subplot(2,1,2),plot(real((out_main)));xlabel('����'),ylabel('����');grid on;title('����������źŵ�ʵ��');


%===============================================================================%
%                    ����ɨ��  Forming beam pattern                             %
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
title('����ǰ����ͼ');
xlabel('������/��');ylabel('��λ��/��');zlabel('����/dB');
% axis([0 60 -60 60 -40 40]);
figure(100)
mesh(p1,a1,y_fuzhu);
title('�������߷���ͼ');
xlabel('������/��');ylabel('��λ��/��');zlabel('����/dB');
figure(6)
plot(p1,y_qian(1+direction1x_signal*180/pi+60,:));
xlabel('������/��');ylabel('����/dB');
grid on;hold on;
figure(7)
plot(a1,y_qian(:,1+direction2x_signal*180/pi/0.5));
xlabel('��λ��/��');ylabel('����/dB');
grid on;hold on;


%% ���Ļ�ͼʹ��
hfig=figure(8);
mesh(p1,a1,y_cancel);
view([90,0]);
title('\fontname{����}����Ӧ�԰���������߷���ͼ', 'FontSize',10.5);
xlabel('\fontname{����}������\fontname{Times new roman}/��');
ylabel('\fontname{����}��λ��\fontname{Times new roman}/��');
zlabel('\fontname{����}����\fontname{Times new roman}/dB');
set(gca,'FontName','Times New Roman','FontSize',10.5);
grid on;
axis tight;
% ����ͼƬ��С��Ӧ��ͨ������õ���Ⱥ͸߶ȡ�
% ����A4ֽ�Ŀ��Ϊ21cm���������ҳ�߾��Ϊ2.5cm�����ĵĿ�Ⱦ���16cm��
% ���ҳ���ǵ�������Ϊ���14cm���߶�8.6cm�ȽϺ��ʣ��ƽ����0.618����
% ˫���Ļ�����7cm����4.3cm���ɡ�
% figWidth = 12;  % ����ͼƬ���
% figHeight = 7.416;  % ����ͼƬ�߶�
% set(hfig,'PaperUnits','centimeters'); % ͼƬ�ߴ����õ�λ��inches��Ӣ�磬centimeters������
% set(hfig,'PaperPosition',[0 0 figWidth figHeight]);
% fileout = 'SLC��������߷���ͼ.'; % ���ͼƬ���ļ���
% print(hfig,[fileout,'tif'],'-r600','-dtiff'); % ����ͼƬ��ʽ���ֱ���
%%


% axis([0 60 -60 60 -40 40]);
figure(9)
plot(p1,y_cancel(1+direction1x_signal*180/pi+60,:),'r');
xlabel('������/��');ylabel('����/dB');
title('��������ͨ������ά����')
grid on
figure(10)
plot(a1,y_cancel(:,1+direction2x_signal*180/pi/0.5),'r');
xlabel('��λ��/��');ylabel('����/dB');
title('��������ͨ����λά����')
grid on
% =====================================================================%
figure(11);
subplot(2,1,1);plot(real(main));grid on;zoom on; %����ǰ�����
ylabel('����');title('����ǰ����');
subplot(2,1,2);plot(real(out_main));grid on;zoom on; %����������
 xlabel('��������'); ylabel('����');title('���������');
%  subplot(3,1,3);plot(abs(SMI_OUT));grid on;zoom on; %SMI���
% ylabel('����');title('SMI���');





