close all;
clear all;
clc; 
%%�źŲ���
A=1;                  %�����źŵ����
Phi0=0;               %�����źŵ��������
T=20e-6;              %�ź�ʱ��
B=10e6;               %�źŴ���
F0=0e6;               %��ƵƵ�ʣ�����ƵƵ��
Fs=40e6;        %����Ƶ��

%%%%% �źŵĲ������� %%%%%
K=B/T;         %��Ƶб��
Ts=1/Fs;       %��������
N=T/Ts;        %��������

%%%%% �������Ե�Ƶ�����ź� %%%%%
t=linspace(-T/2,T/2,N);
freq=linspace(-Fs/2,Fs/2,N);
sn=A*exp(1j*(2*pi*F0*t+pi*K*t.^2+Phi0));

%%%%% �����ܼ����Ŀ�� %%%%%
L=2e-6; % ��Ŀ����
M=20;  % ��Ŀ�����
k=T/L; % �ֶ���Ŀ
duanNum=k+M-1;  % �ع����ܶ���
T_g=duanNum*L;  % �ع���ʱ�䳤��
gn=zeros(1,T_g*Fs);
gn(1,1:N)=sn;
for i=1:M-1
    cur = gn(1,i*L*Fs+1:i*L*Fs+N);
    gn(1,i*L*Fs+1:i*L*Fs+N)=cur+A*exp(1j*(2*pi*F0*(t)+pi*K*(t).^2+Phi0));
end

%%%%% ���ܼ����Ŀ������źŽ�����ѹ %%%%%
S_rec = gn;
S_ref = A*exp(1j*(2*pi*F0*t+pi*K*t.^2+Phi0));
temp = conv(S_rec,conj(fliplr(S_ref)));
S_pc = temp(length(S_ref):end);
                                     
hfig = figure(1);
subplot(3,1,1);
plot((t+T/2)*1e6,real(sn));
title('\fontsize{10.5}\fontname{����}���Ե�Ƶ�ź�ʵ��');
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10.5);
grid on;
axis tight;
xlabel('\fontsize{10.5}\fontname{����}ʱ��\fontname{Times New Roman}/us', 'FontSize', 10.5);
ylabel('\fontsize{10.5}\fontname{����}����\fontname{Times New Roman}/V', 'FontSize', 10.5)

subplot(3,1,2);
plot(linspace(0,T_g,T_g*Fs)*1e6,abs(gn));
title('\fontsize{10.5}\fontname{����}�ܼ����Ŀ������źŰ���');
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10.5);
grid on;
axis tight;
xlabel('\fontsize{10.5}\fontname{����}ʱ��\fontname{Times New Roman}/us', 'FontSize', 10.5);
ylabel('\fontsize{10.5}\fontname{����}����\fontname{Times New Roman}/V', 'FontSize', 10.5);

subplot(3,1,3);
plot(linspace(0,T_g,T_g*Fs)*1e6,abs(S_pc));
title('\fontsize{10.5}\fontname{����}�ܼ����Ŀ������ź���ѹ���');
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10.5);
grid on;
axis tight;
xlabel('\fontsize{10.5}\fontname{����}ʱ��\fontname{Times New Roman}/us', 'FontSize', 10.5);
ylabel('\fontsize{10.5}\fontname{����}����\fontname{Times New Roman}/V', 'FontSize', 10.5);

% ����ͼƬ�����ʽ
figWidth = 14;
figHeight = 12.9;  % ��������14cm,8.6cm���ƽ����Ϊ0.618
set(hfig, 'PaperUnits', 'centimeters');  %centimeters %inches
set(hfig, 'PaperPosition', [0 0 figWidth figHeight]);
print(hfig, ['�ܼ����Ŀ�����.', 'tif'], '-r600', '-dtiff');   
% print(hfig, ['�ܼ����Ŀ�����.', 'jpg'], '-r600', '-djpeg');  