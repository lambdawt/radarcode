clear all;close all;
azi_num=2000;
fr=1000;

fr=12e3;
sigmav=40;     %�Ӳ�����
sigmaf=2*sigmav*6e9/3e8;

xi=randn(1,azi_num);  %�����Ҷ����ĸ�˹����N(0��1)
xq=randn(1,azi_num);
%�γ��˲���Ƶ����Ӧ

FreXi = fftshift(fft(xi,2048));
f = (-length(FreXi)/2:1:length(FreXi)/2-1)/length(FreXi)*fr;
hf=exp(-f.^2/(2*sigmaf^2));
FreXi = FreXi.*hf;
xxi = ifft(fftshift(FreXi));
xxi = xxi(1:azi_num);
FreXq = fftshift(fft(xq,2048));
FreXq = FreXq.*hf;
xxq = ifft(fftshift(FreXq));
xxq = xxq(1:azi_num);

% lamda0=0.05;
% sigmav=0.1;
% sigmaf=2*sigmav/lamda0;
% 
% rand('state',sum(100*clock));
% d1=rand(1,azi_num);
% rand('state',7*sum(100*clock)+3);
% d2=rand(1,azi_num);
% xi=2*(sqrt(-2*log(d1)).*cos(2*pi*d2));
% xq=2*sqrt(-2*log(d1)).*sin(2*pi*d2);
% 
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
% xxi=xxi(coe_num*2+1:azi_num+coe_num*2);
% xxq=xxq(coe_num*2+1:azi_num+coe_num*2);

xisigmac=std(xxi);
ximuc=mean(xxi);
yyi=(xxi-ximuc)/xisigmac;
xqsigmac=std(xxq);
xqmuc=mean(xxq);
yyq=(xxq-xqmuc)/xqsigmac;

p=1.5;%��״����
q=2.2;%�߶Ȳ���
sigmac=sqrt((q.^p)/2);
yyi=sigmac*yyi;
yyq=sigmac*yyq;

% xdata=(yyi.*yyi+yyq.*yyq).^(1/p);
xdata=((yyi + 1j*yyq).^2).^(1/p);

figure,plot(abs(xdata));xlabel('ʱ��');ylabel('����');title('Τ�����ֲ��Ӳ�ʱ����');

num=100;
maxdat=max(abs(xdata));
mindat=min(abs(xdata));
NN=hist(abs(xdata),num);
xpdf1=num*NN/((sum(NN)*(maxdat-mindat)));
xaxis1=mindat:(maxdat-mindat)/num:maxdat-(maxdat-mindat)/num;
th_val=p*(xaxis1.^(p-1)).*exp(-(xaxis1/q).^p)./(q.^p);

figure;plot(xaxis1,xpdf1);
hold,plot(xaxis1,th_val,':r');
title('�Ӳ����ȷֲ�');xlabel('����');ylabel('�����ܶ�');
signal=xdata;
signal=signal-mean(signal);

figure,M=256;
psd_dat=pburg(real(signal),16,M,fr);
psd_dat=psd_dat/(max(psd_dat));
freqx=0:0.5*M;
freqx=freqx*fr/M;
psd_dat = fftshift(fft(signal,2048));
freqx=f;

plot(freqx,abs(psd_dat));title('�Ӳ�Ƶ��');xlabel('Ƶ��(Hz)');ylabel('�������ܶ�');
powerf=exp(-freqx.^2/(2*sigmaf.^2))*200;
hold;plot(freqx,powerf,':r');