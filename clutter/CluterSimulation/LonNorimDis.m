clear all;close all;
azi_num=2000;
fr=1000;

lamda0=0.05;
sigmav=1.0;
sigmaf=2*sigmav/lamda0;

rand('state',sum(100*clock));
d1=rand(1,azi_num);
rand('state',7*sum(100*clock)+3);
d2=rand(1,azi_num);
xi=2*(sqrt(-2*log(d1)).*cos(2*pi*d2));
xq=2*sqrt(-2*log(d1)).*sin(2*pi*d2);

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
%Gaussion clutter generation
% xxi=conv(b,xi);
% xxi=xxi(coe_num*2+1:azi_num+coe_num*2);
% 
% 
% xisigmac=std(xxi);
% ximuc=mean(xxi);
% yyi=(xxi-ximuc)/xisigmac;
% 
% xxi=conv(b,xq);
% xxi=xxi(coe_num*2+1:azi_num+coe_num*2);
% xisigmac=std(xxi);
% ximuc=mean(xxi);
% yyq=(xxi-ximuc)/xisigmac;
yyi = xi;yyq= xq;

yyi = yyi + 1j*yyq;

% yyi = randn(1,2048) +1j * randn(1,2048);

muc=1.5;
sigmac=0.6;
yyi=sigmac*yyi+log(muc);
xdata=exp(yyi);
figure,plot(abs(xdata));xlabel('ʱ��');ylabel('����');title('������̬�ֲ��Ӳ�ʱ����');
num=100;
maxdat=max(abs(xdata));
mindat=min(abs(xdata));
NN=hist(abs(xdata),num);
xpdf1=num*NN/((sum(NN)*(maxdat-mindat)));
xaxis1=mindat:(maxdat-mindat)/num:maxdat-(maxdat-mindat)/num;
th_val=lognpdf(xaxis1,log(muc),sigmac);

figure;plot(xaxis1,xpdf1);
hold,plot(xaxis1,th_val,':r');
title('�Ӳ����ȷֲ�');xlabel('����');ylabel('�����ܶ�');
signal=xdata;
signal=signal-mean(signal);

figure,M=128;
psd_dat=pburg(real(signal),16,M,fr);
psd_dat=psd_dat/(max(psd_dat));
freqx=0:0.5*M;
freqx=freqx*fr/M;
plot(freqx,psd_dat);title('�Ӳ�Ƶ��');xlabel('Ƶ��(Hz)');ylabel('�������ܶ�');
powerf=exp(-freqx.^2/(2*sigmaf.^2));
hold;plot(freqx,powerf,':r');