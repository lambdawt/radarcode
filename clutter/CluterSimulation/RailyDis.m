clear all��close all;
azi_num=2000;   %ȡ2000����
fr=1000;        %�״��ظ�Ƶ��

lamda0=0.05;   %�Ӳ�����
sigmav=1.0;     %�Ӳ�����
sigmaf=2*sigmav/lamda0;  


rand('state',sum(100*clock)); %��������U(0,1)�ֲ����������
d1=rand(1,azi_num);            
rand('state',7*sum(100*clock)+3);
d2=rand(1,azi_num);
xi=2*sqrt(-2*log(d1)).*cos(2*pi*d2);  %�����Ҷ����ĸ�˹����N(0��1)
xq=2*sqrt(-2*log(d1)).*sin(2*pi*d2);
%�γ��˲���Ƶ����Ӧ
coe_num=12;           %���˲���ϵ�����ø���Ҷ����չ����
for n=0:coe_num
    coeff(n+1)=2*sigmaf*sqrt(pi)*exp(-4*sigmaf^2*pi^2*n^2/fr^2)/fr;  
end
for n=1:2*coe_num+1
    if n<=coe_num+1
        b(n)=1/2*coeff(coe_num+2-n);
    else
        b(n)=1/2*coeff(n-coe_num);
    end
end
%���ɸ�˹���Ӳ�
xxi=conv(b,xi);   
xxq=conv(b,xq);   
xxi=xxi(coe_num*2+1:azi_num+coe_num*2);%Ŀ����ȥ����̬��Ӧ
xxq=xxq(coe_num*2+1:azi_num+coe_num*2);
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
% 
% clear;
% azi_num=2000;
% sigmac=1.2 ; 
% rand('state',sum(100*clock)); %��������U(0,1)�ֲ����������
% d1=randn(1,azi_num);            
% rand('state',7*sum(100*clock)+3);
% d2=randn(1,azi_num);
% d1f = fft(d1);d1f(100:1900) = 0;d1 = ifft(d1f);
% d2f = fft(d2);d1f(100:1900) = 0;d2 = ifft(d2f);
% xisigmac=std(d1);     
% ximuc=mean(d1);       
% yyi=(d1-ximuc)/xisigmac;    
% xqsigmac=std(d2);     
% xqmuc=mean(d2);       
% yyq=(d2-xqmuc)/xqsigmac;    %��һ��
% 
% 
% ydata = yyi + 1j* yyq;
% 
% % ydata = randn(1,2000) + 1j* randn(1,2000);
% 
% num=100;                   %������ܶȺ����Ĳ���
% maxdat=max(abs(ydata));
% mindat=min(abs(ydata));
% NN=hist(abs(ydata),num);   
% xpdf1=num*NN/((sum(NN))*(maxdat-mindat));        %��ֱ��ͼ���Ƶĸ����ܶȺ���
% xaxisl=mindat:(maxdat-mindat)/num:maxdat-(maxdat-mindat)/num;  
% th_val=(xaxisl./sigmac.^2).*exp(-xaxisl.^2./(2*sigmac.^2));   %�����ܶȺ�������ֵ
% figure;
% plot(xaxisl,xpdf1);               %�����������ĸ����ܶȺ�������
% hold on;plot(xaxisl,th_val,'r:'); %�������۸����ܶȺ�������
% title('�Ӳ����ȷֲ�');
% xlabel('����');ylabel('�����ܶ�');

signal=ydata;
signal=signal-mean(signal);      %�������ܶȣ���ȥ��ֱ������
figure;M=256;                     %��burg�����ƹ������ܶ�
psd_dat=pburg(real(signal),32,M,fr);
psd_dat=psd_dat/(max(psd_dat));   %��һ������
freqx=0:0.5*M;
freqx=freqx*fr/M;
plot(freqx,psd_dat);title('�Ӳ�Ƶ��');
xlabel('Ƶ��/HZ');ylabel('�������ܶ�');

%���������˹������
powerf=exp(-freqx.^2/(2*sigmaf.^2));
hold on;plot(freqx,powerf,'r:');