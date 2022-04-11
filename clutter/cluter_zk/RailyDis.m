clear all;close all;clc
azi_num=1920;   %ȡ2000����
fr=12e3;

% sigmav=40*0.03;     %�Ӳ�����
sigmav = 40;
sigmaf=2*sigmav*6e9/3e8;

% fr=1000;
% lamda0=0.05;
% sigmav=0.1;
% sigmaf=2*sigmav/lamda0;

xi=randn(1,azi_num);  %�����Ҷ����ĸ�˹����N(0��1)
xq=randn(1,azi_num);
%�γ��˲���Ƶ����Ӧ

FreXi = fftshift(fft(xi,2048));
f = (-length(FreXi)/2:1:length(FreXi)/2-1)/length(FreXi)*fr;
hf=exp(-f.^2/(2*sigmaf^2));% ��˹��
% hf=1./(1+abs(f/sigmaf).^4);% ȫ����
% hf=exp((-abs(f)/sigmaf));%ָ����3

FreXi = FreXi.*hf;
xxi = ifft(fftshift(FreXi));
xxi = xxi(1:azi_num);
FreXq = fftshift(fft(xq,2048));
FreXq = FreXq.*hf;
xxq = ifft(fftshift(FreXq));
xxq = xxq(1:azi_num);

% coe_num=12;
% for n=0:coe_num
%     coeff(n+1)=2*sigmaf*sqrt(pi)*exp(-4*sigmaf^2*pi^2*n^2/fr^2)/fr;
% end
% for n=1:2*coe_num+1
%     if n<=coe_num+1
%         b(n)=1/2*coeff(coe_num+2-n);
%     else
%         b(n)=1/2*coeff(n-coe_num);
%     end
% end
% %Gaussion clutter generation
% xxi=conv(b,xi);
% xxq=conv(b,xq);

xisigmac=std(xxi);     
ximuc=mean(xxi);       
yyi=(xxi-ximuc)/xisigmac;    
xqsigmac=std(xxq);     
xqmuc=mean(xxq);       
yyq=(xxq-xqmuc)/xqsigmac;    %��һ��
sigmac=1.2 ;           %�Ӳ��ı�׼��
yyi=sigmac*yyi;        %ʹ�����ֲ��Ӳ�����ָ���ı�׼��
yyq=sigmac*yyq;        %ʹ�����ֲ��鲿�Ӳ�
ydata=yyi+j*yyq;       %�����ֲ��Ӳ��γ�
figure;
subplot(211);plot(real(ydata)); %�����ֲ��Ӳ�ʵ��
title('�����Ӳ�ʱ���Σ�ʵ��');  
subplot(212);plot(imag(ydata)); %�����ֲ��Ӳ��鲿
title('�����Ӳ�ʱ���Σ��鲿')

num=100;                   %������ܶȺ����Ĳ���
maxdat=max(abs(ydata));
mindat=min(abs(ydata));
NN=hist(abs(ydata),num);   
xpdf1=num*NN/((sum(NN))*(maxdat-mindat));        %��ֱ��ͼ���Ƶĸ����ܶȺ���
xaxisl=mindat:(maxdat-mindat)/num:maxdat-(maxdat-mindat)/num;  
th_val=(xaxisl./sigmac.^2).*exp(-xaxisl.^2./(2*sigmac.^2));   %�����ܶȺ�������ֵ
figure;
plot(xaxisl,xpdf1);               %�����������ĸ����ܶȺ�������
hold on;plot(xaxisl,th_val,'r:'); %�������۸����ܶȺ�������
title('�Ӳ����ȷֲ�');
xlabel('����');ylabel('�����ܶ�');

signal=ydata;
signal=signal-mean(signal);      %�������ܶȣ���ȥ��ֱ������
figure;M=256;                     %��burg�����ƹ������ܶ�
% psd_dat=pburg(abs(signal),32,M,fr);
% psd_dat=psd_dat/(max(psd_dat));   %��һ������
% freqx=0:0.5*M;
% freqx=freqx*fr/M;
psd_dat = fftshift(fft(signal,2048));
freqx=f;

plot(freqx,abs(psd_dat));title('�Ӳ�Ƶ��');
xlabel('Ƶ��/HZ');ylabel('�������ܶ�');

%���������˹������
% powerf=exp(-freqx.^2/(2*sigmaf^2))*200;%��˹
% powerf=1./(1+abs(f/sigmaf).^3)*200;%ȫ����
powerf=exp((-abs(f)/sigmaf))*200;%ָ����
v= freqx*3e8/2/6e9;
hold on;plot(freqx,powerf,'r:');

Acc=fftshift(fft(ydata));
v = (-length(Acc)/2:1:length(Acc)/2-1)/length(Acc)*fr*3e8/2/6e9;