% ���ת��
%��������
clc;
clear;
close all;

%% ����
c = 3e8;          % ���٣�m/s��
k = 1.38*1e-23;     %��������������J/K��
To = 290;               % ��׼���� 17�棨K �����ģ�.

%% �״�ϵͳ����
M=30;
fo = 1.2*1e9;         % ��Ƶ(Hz)
Pt = 15*1e3;          % �������ֵ����(w��
Gt = 33;              % ��������棨dB��
B = 10*1e6;            % ���ջ�����
fs = 2*B;             %����Ƶ��
fr = 2000;            % �����ظ�Ƶ��(Hz)
Tr = 1/fr;            % �����ظ����� ���룩
Tp = 2e-5;              % ���� ���룩
lambda = c/fo;           % ���� ��m��
d = lambda/2;            % ��Ԫ��� ��m��
Gele = 5 ;              %��Ԫ�������棨dB��
F = 12;                 %���ջ�����ϵ����dB��
Ts = To*(10^(F/10)-1);      %���ջ��¶� (K)
Nn = k*Ts;                 % ���ջ������������ܶ� ��W/Hz��
sigmaN2 = 1;               % �������ʣ�W��
dR  = c/2/B;               % �������� ��m��    

%**************************ԭʼ�ź���ʽ
K = B/Tp;  %���Ե�Ƶб��
t0 = 0:1/fs:Tp;
s0 = exp(1i*pi*K*t0.^2);
t = [t0,Tp+1/fs:1/fs:Tr-1/fs];
s = exp(1i*pi*K*t.^2).*rectpuls(t-Tp/2,Tp);%�������Ե�Ƶ�����ź�
Ns = length(t);            %һ���ڲ�������

figure(1)
subplot(2,1,1)
plot(t,real(s));
subplot(2,1,2)
plot(t(1:Tp/Tr*Ns),real(s(1:Tp/Tr*Ns)))
%****************************Ŀ��ز�
%% Ŀ���ź�
R_t = 35e3; 
vt = 35;
fdt = 2*vt/lambda;
tao = 2*R_t/c;
K = B/Tp;  %���Ե�Ƶб��
sr = exp(1i*pi*K*(t-tao).^2).*rectpuls(t-Tp/2-tao,Tp);

figure(2)
subplot(2,1,1)
plot(t,real(sr));
subplot(2,1,2)
plot(t(tao/Tr*Ns:(tao+Tp)/Tr*Ns),real(sr(tao/Tr*Ns:(tao+Tp)/Tr*Ns)));

%ƥ���˲���
s0 = exp(1i*pi*K*t0.^2);
h = conj(fliplr(s0));
H = fft(h,Ns);

S0 = repmat(sr,1,M);
SigAll = S0.*exp(1i*2*pi*fdt*(0:Ns*M-1)/fs);
SigAll= reshape(SigAll,Ns,M);

%% ��Ъ����ת������

TTs = 4e-6;       %�����ź�pri
tau = TTs/2;      %�����ź�Tp
tau1 = TTs;
pt0 = square(2*pi*(t+tau/2)/tau1,50); %�������ڷ����ź�
pt= (pt0+1)/2;       %���ɲ����ź�

figure(11);
plot(t,real(pt));

s11 = sr.*pt;
s11=[zeros(1,ceil(k*tau*fs)),s11(1:Ns-ceil(k*tau*fs))];

% s11 = s11*1;
S1 = repmat(s11,1,M);
JamAll = S1.*exp(1i*2*pi*fdt*(0:Ns*M-1)/fs);
JamAll= reshape(JamAll,Ns,M);


figure(3)
subplot(2,1,1)
plot(t,real(sr),'k');
xlabel('ʱ�䣨s��','FontSize',19,'FontName','����','Color','k');
ylabel('���ȣ�V��','FontSize',19,'FontName','����','Color','k');
subplot(2,1,2)
plot(t,real(JamAll(:,1)),'k');
xlabel('ʱ�䣨s��','FontSize',19,'FontName','����','Color','k');
ylabel('���ȣ�V��','FontSize',19,'FontName','����','Color','k');
% hold on
% plot(t,real(SigAll(:,1)),'r:')
% ylim([-5 5]);
% figure(2)
% Yt=fftshift(fft(SigAll(:,1)));
% f=linspace(-fs/2,fs/2,length(t));
% plot(f,abs(Yt))


Jam_pp = zeros(Ns,M);
Sig_pp = zeros(Ns,M);
for m=1:M
    Jam_fft = fft(JamAll(:,m).');
    Sig_fft = fft(SigAll(:,m).');
    Jam_pp(:,m) = ifft(H.*Jam_fft);
    Sig_pp(:,m) = ifft(H.*Sig_fft);
end

% save('jxcy_zr','Jam_pp','Sig_pp');

mtd1 = fftshift(fft(Jam_pp,16,2),2);
mtd1 = mtd1/max(max(abs(mtd1)));
mtd2 = fftshift(fft(Sig_pp,16,2),2);
mtd2 = mtd2/max(max(abs(mtd2)));
figure(4);
mesh(abs(mtd1));
hold on
mesh(abs(mtd2));
shading interp
colormap jet
axis tight
xlabel('������ͨ��','FontSize',19,'FontName','����','Color','k');
ylabel('������','FontSize',19,'FontName','����','Color','k');
zlabel('��һ������','FontSize',19,'FontName','����','Color','k');

All = [Sig_pp(:,1);Jam_pp(:,1)];
figure(5)
plot(t,abs(Sig_pp(:,1))/max(abs(All)),'r--');
hold on
plot(t,abs(Jam_pp(:,1)/max(abs(All))),'k');
% axis([1.8e-4 2.5e-4 0 1200])
% legend('Ŀ��ز�','��Ъ��������','FontSize',19)
xlabel('ʱ�䣨s��','FontSize',19,'FontName','����','Color','k');
ylabel('���ȣ�V��','FontSize',19,'FontName','����','Color','k')










