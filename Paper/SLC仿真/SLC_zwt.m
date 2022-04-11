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
N=3;   %��������Ԫ��  
auN=2; %�������߸���
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

figure(1);
subplot(2,1,1);
plot(real(signal));grid on;zoom on;
xlabel('��������'),ylabel('����'),title('ʵ��');
subplot(2,1,2);
plot(imag(signal));grid on;zoom on;
xlabel('��������'),ylabel('����'),title('�鲿');
% signal_fft=fft(signal+Interfer1+Interfer2);
% figure(2);
% plot(Fs/number*(0:number-1),abs(signal_fft));grid on;
% xlabel('Ƶ��(��λ��Hz��'),ylabel('����'),title('���Ե�Ƶ�ź�Ƶ��');
%====================================================================================%
%                                   ����ϵͳ�����ź�                                 %
%====================================================================================%
Systemnoiseamp=1;
SLCSystemNoise=Systemnoiseamp*randn(auN,number); %  mean value is 0, ������������
DBFSystemNoise=Systemnoiseamp*randn(N,number);  %����������
%====================================================================================%
%                                 �γ�ͨ������                                       %
%=====================================================================================%
INR1=20;   %dB ����� ����
INR2=20;   %dB
SNR=10;    %dB�����
% f=200e6; %�ز�����Ƶ��
wavelength=0.1;%�ز��źŲ���:10cm
%===================== ����1���� ===============%
direction1x_interfere=41*pi/180; 
%==================== ����2���� ================%
direction2x_interfere=-30*pi/180;
%====================Ŀ���źŷ���===============%
directionx_signal=10*pi/180; %��λ��
dx=wavelength/2; %��Ԫ���
dass=[-1*dx 0*dx]; %�������ߵ�����    ��ʾ
% d=[1*dx 2*dx 3*dx 4*dx 5*dx 6*dx 7*dx 8*dx 9*dx 10*dx 11*dx 12*dx];% 13*dx 14*dx...
  %15*dx 16*dx 17*dx 18*dx 19*dx 20*dx 21*dx];% �����ߵ������ʾ
d=dx:dx:N*dx;
k1=sin(direction1x_interfere);%����1
steer_interfer1=exp(j*(2*pi/wavelength*k1*d));%�������еĸ���1����ʸ��:��������ʸ����������� 
auxsteer_interfer1=exp(j*(2*pi/wavelength*k1*dass));%�����еĸ���1����ʸ��  
k2=sin(direction2x_interfere);%����2�Ĵ�������ʸ��
steer_interfer2=exp(j*(2*pi/wavelength*k2*d));%�������еĸ���2����ʸ��  
auxsteer_interfer2=exp(j*(2*pi/wavelength*k2*dass));%���������еĸ���2����ʸ��
k3=sin(directionx_signal);%�źŵĴ�������ʸ��
steer_signal=exp(j*(2*pi/wavelength*k3*d)); %�������źŵĵ���ʸ��
auxsteer_signal=exp(j*(2*pi/wavelength*k3*dass)); %���������е��źŵ���ʸ��
%
interfer1=10^(INR1/10)*Interfer1; %������20*log10(x)
interfer2=10^(INR2/10)*Interfer2; %������10*log10(x)
signal=10^(SNR/10)*signal;
% signal_fft=fft(signal+interfer1+interfer2);
% figure(2);
% plot(Fs/number*(0:number-1),abs(signal_fft));grid on;
% xlabel('Ƶ��(��λ��Hz��'),ylabel('����'),title('���Ե�Ƶ�ź�Ƶ��');
%==========================================================================%
%              interference  echo   
%==========================================================================%
input_main=steer_interfer1.'*interfer1 + steer_interfer2.'*interfer2 +...
           DBFSystemNoise; % �����߽��յ����ź�(����1+����2+������
input_aux=auxsteer_interfer1.'*interfer1 + auxsteer_interfer2.'*interfer2 +...
         SLCSystemNoise; % �������߽��յ����ź�     
input_main1=steer_interfer1.'*interfer1 + steer_interfer2.'*interfer2 +...
           steer_signal.'*signal+DBFSystemNoise; % �����߽��յ����ź�(����1+����2+�ź�+������
input_aux1=auxsteer_interfer1.'*interfer1 + auxsteer_interfer2.'*interfer2 +...
           auxsteer_signal.'*signal+SLCSystemNoise;%     
%========================SMI===============================================
Rxx=input_main*input_main'/number;
SMI_W=inv(Rxx)*steer_signal.';%Ȩʸ��������Լ����С����(LCMV)׼�� 
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
% main= (steer_signal.'.*tay)'*input_main; %�����е����
main= (normal_W)'*input_main; %�����е����
[size_x,size_y]=size(input_aux); %���ظ������ߵĸ����Ͳ�������
r_xd=zeros(1,size_y); R=zeros(size_x,size_x);
% R=input_aux*input_aux'/number; %������������ؾ���
R=zeros(auN,auN);
r_xd=zeros(auN,1);
SLCW=zeros(auN,1);
M=50;%:20:1000;%������
% NN=number/M;%ÿ�����ĵĵ���
% ii=1:length(M)
for i=1:M;
    R=input_aux(:,i)*input_aux(:,i)'+R;%%%%%%%%%%�Ǽ�������ͨ����ֵһ�����������غ������Ǹ������߸��Ե�ֵ�����𣿣���
    r_xd=input_aux(:,i)*main(:,i)'+r_xd;
end
    R=R/M;
    r_xd=r_xd/M;
% r_xd=input_aux*main'/number; %�������븨�����ߵĻ��������
SLCW=inv(R)*r_xd; %��������
% r_xd=input_aux*main'/number; %�������븨�����ߵĻ��������
% SLCW=inv(R)*r_xd; %�����԰��������Ȩ(����������ź�d������������) 
out_aux=SLCW'*input_aux1; %�����������
main=(normal_W)'*input_main1; %���������   
% main2=(steer_signal.')'*(steer_signal.'*signal+DBFSystemNoise);
out_main=main-out_aux;  %���������

out_aux1=SLCW'*input_aux;
main1=(normal_W)'*input_main;
out_main1=main1-out_aux1;

p_a=norm(main1).^2/number;%����ǰƽ������
Es0=10*log10(p_a);
p_b=norm(out_main1).^2/number;  %������ƽ������
Es1=10*log10(p_b);
CR1=10*log10(p_a/p_b)%������

% figure(100)
% plot(real(main1));grid on;hold on;
% plot(real(aux),'k');hold on;
% plot(real(out_main),'r');
% legend('����ǰ������','��������','������������');

figure(2);
plot(Fs/number*(0:number-1),abs(fft(main)),'r');grid on;hold on;
plot(Fs/number*(0:number-1),(abs(fft(out_main)))-850*ones(1,length(out_main)));grid on;
xlabel('Ƶ�ʣ�Hz��');
ylabel('����');
title('����ǰ���ź�Ƶ��');
legend('����ǰ','������');


%========================������ͨ�����Ÿ����========================================
S_main=norm((normal_W)'*(steer_signal.'*signal)).^2/number;
I_main=norm((normal_W)'*(steer_interfer1.'*interfer1 + steer_interfer2.'*interfer2)).^2/number;
N_main=norm((normal_W)'*DBFSystemNoise).^2/number;
SNR_main=10*log10(S_main/N_main);
INR_main=10*log10(I_main/N_main);
%========================������ͨ�����Ÿ����=======================================

%===============================================================================%
%                    ����ɨ��  Forming beam pattern                             %
%===============================================================================%


row=0; %��
for loop_v=-90*pi/180:0.2*pi/180:90*pi/180       
     row=row+1
     k= sin(loop_v); %��������ʸ��
     steer_inter=exp(j*(2*pi/wavelength*k*d)).'; %ɨ�赼��ʸ��,searching variable
%      SMI_W=SMI_W.*hamming(length(SMI_W));
     P_SMI(row)=abs((SMI_W)'*steer_inter); % SMI DBF(����ͼ:�������źŴ���3)
     P_normal(row)=((normal_W)'*steer_inter);
%    P_normal(row)=abs((normal_W.*tay)'*steer_inter);  % Usual DBF 
end
u=[-90*pi/180:0.2*pi/180:90*pi/180]*180/pi;
%%==========================SMI�����߷���ͼ====================================%
%%%%%%%%%%%%%%%%%%%%%%%%%%%SMI�����߷���ͼ��ʲô������%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3);
P_SMI_1= P_SMI./(max(P_SMI)); %��һ��
P_SMI_1=20*log10(P_SMI_1); %dB
plot(u,P_SMI_1);grid on;zoom on; %����һ�����SMI��ͼ��ָ���ź�
axis([-100 100 -60 0]);
xlabel('�Ƕ�');ylabel('����(dB)');
title('SMI����ͼ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%SMI�����߷���ͼ��ʲô������%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%========================����ǰ�ķ���ͼ=======================================%
figure(4);
% P_normal_1= P_normal./(max( P_normal));%��һ��
P_normal_1=abs( P_normal);
% P_normal_1=20*log10( P_normal_1./(max( P_normal_1))); %dB
P_normal_1=20*log10( P_normal_1); %dB
plot(u, P_normal_1);grid on;zoom on; %������ǰ�ķ���ͼ
% axis([-100 100 -60 0]); %�޶�����
xlabel('�Ƕ�');ylabel('����(dB)');
title('����ǰ����ͼ');
% legend('3�ȷ���ͼ');
% %============================================================================%
% %                     Figure  SLC                                            %
% %============================================================================%
figure(5);
subplot(2,2,1);plot(real(main));grid on;zoom on; %����ǰ�����
xlabel('��������');ylabel('����');title('����ǰ����źŵ�ʵ��');
subplot(2,2,3);plot(real(out_main));grid on;zoom on; %����������
xlabel('��������');ylabel('����');title('����������źŵ�ʵ��');
subplot(2,2,2);plot(imag(main));grid on;zoom on; %����ǰ�����
xlabel('��������');ylabel('����');title('����ǰ����źŵ��鲿');
subplot(2,2,4);plot(imag(out_main));grid on;zoom on; %����������
xlabel('��������');ylabel('����');title('����������źŵ��鲿');
 
figure(8);
subplot(2,1,1);plot(20*log10(abs(main)));grid on;zoom on; %����ǰ�����
xlabel('��������');ylabel('����(dB)');title('����ǰ����źŵķ���');
subplot(2,1,2);plot(20*log10(abs(out_main)));grid on;zoom on; %����������
 xlabel('��������'); ylabel('����(dB)');title('����������źŵķ���');

%  subplot(3,1,3);plot(abs(SMI_OUT));grid on;zoom on; %SMI���
% ylabel('����');title('SMI���');
%==============================������============================================%
% CR=max(abs(main))./max(abs(out_main));
% CR=20*log10(CR)
p_a=norm(main).^2/number;%����ǰƽ������
Es0=10*log10(p_a)
p_b=norm(out_main).^2/number;  %������ƽ������
Es1=10*log10(p_b)
CR=10*log10(p_a/p_b)%������
%%=============================�������߷���ͼ====================================%
row=0; %��
for loop_v=-90*pi/180:0.2*pi/180:90*pi/180       
     row=row+1;
     k= sin(loop_v); %��������ʸ��
    steer_inter=exp(j*(2*pi/wavelength*k*dass)).'; %ɨ�赼��ʸ��,searching variable
     P_aux(row)=(SLCW'*steer_inter);%�������߷���ͼ
  end
%   P_aux_1=P_aux./max(P_normal);
  P_aux_1=abs(P_aux);
  P_aux_1=20*log10(P_aux_1);
figure(6);
plot(u,P_aux_1,'r');grid on;zoom on;hold on;
plot(u,P_normal_1);
% axis([-100 100 -60 0]);
xlabel('�Ƕ�');ylabel('����(dB)');
legend('����ǰ����ͼ','�������߷���ͼ')
title('�������߷���ͼ');
%%=============================�����ķ���ͼ====================================%
P =abs(P_normal-P_aux);
% P_1=P./(max(P_normal)); % ��һ��
% P_1=P./(max(P)); % ��һ��
% P_2=P;
P_1=20*log10(P); %dB ������20*log10(x)
% P_2=20*log10(P_2);
figure(7);
plot(u,P_1,'r');grid on;zoom on;hold on;
% plot(u,P_aux_1,'g');grid on;zoom on;hold on;
% plot(u,P_normal_1);
xlabel('�Ƕ�');ylabel('����(dB)');
title('��������ͼ');
% text(30,-3,'���ŷ���Ϊ35��-25��');
% legend('������','��������','����ǰ');
%=================================END=============================================%
%  abc=P_1(35)-P_normal_1(35);
h=conj(fliplr(signal)); 
% %----��ƥ���˲���ϵ�����г�ȡ
% d1=fix(number/a);
% for m=1:d1
%     hh(m)=h(a*(m-1)+1);%----���˲������г�ȡ
% end
% h1=hh./sqrt(sum(abs(hh).^2)); %ƥ���˲���ϵ����һ��
h1=h./sqrt(sum(abs(h).^2));%ƥ���˲���ϵ����һ��
Match_out=conv(h1,out_main);
figure(14)
plot(20*log10(abs(Match_out)/max(abs(Match_out))));grid on;
xlabel('����');
ylabel('����');
title('ƥ���˲�����');

t=linspace(-TimeWidth/2,TimeWidth/2,length(h1));
K1=0.08;n=2;
W=K1+(1-K1)*(cos(pi*t/TimeWidth)).^n;%----������
% W=0.42+0.5*cos(2*pi*t/TimeWidth)+0.08*cos(4*pi*t/TimeWidth);%----Blackman��
Wt=W.*h1;
Wt=Wt./sqrt(sum(abs(Wt).^2));
Window_out=conv(Wt,out_main);
figure(15)
plot(20*log10(abs(Window_out)/max(abs(Window_out))));grid on;
xlabel('����');
ylabel('����');
title('��Hamming������');

RMS=rms1(Window_out);
q=20*log10(max(abs(Match_out))/max(abs(Window_out)));