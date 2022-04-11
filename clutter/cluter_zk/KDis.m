clear all;
close all;
clc
azi_num=2048;


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
% sigmav=1.0;
% sigmaf=2*sigmav/lamda0;
% 
% % rand('state',sum(100*clock));
% d1=rand(1,azi_num);
% % rand('state',7*sum(100*clock)+3);
% d2=rand(1,azi_num);
% xi=2*(sqrt(-2*log(d1)).*cos(2*pi*d2));
% xq=2*(sqrt(-2*log(d1)).*sin(2*pi*d2));
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
% xxi=conv(b,xi);
% xxi=xxi(coe_num*2+1:azi_num+coe_num*2);
% xxq=conv(b,xq);
% xxq=xxq(coe_num*2+1:azi_num+coe_num*2);
xisigmac=std(xxi);
ximuc=mean(xxi);
xxi=(xxi-ximuc)/xisigmac;
xqsigmac=std(xxq);
xqmuc=mean(xxq);
xxq=(xxq-xqmuc)/xqsigmac;


xdata=xxi+1j*xxq;
tmpdat=randn(1,azi_num)+1j*randn(1,azi_num);

% [b,a]=butter(5,0.01);
% sk_dat=filter(b,a,tmpdat);
tmpdatf = fft(tmpdat);
tmpdatf(100:1900) = 0;
sk_dat = ifft(tmpdatf);

sk_dat=sk_dat/std(sk_dat);
% figure;plot(abs(sk_dat));
%%%%%%%%%%%%%%����ĳ��������Է���%%%%%%%%%%%%%%%
max_z=6;
step=0.005;
table_z=0:step:max_z;
vmuc=2;
table_s=nonline_eq_sirp(table_z,vmuc);
% figure;plot(abs(table_s));
for n=1:azi_num
   index(n)=floor(abs(sk_dat(n))/max_z*length(table_z)+1);%length���鳤��
   sk_dat(n)=table_s(index(n));
end
ydata=xdata.*sk_dat;
% figure;plot(abs(index));
figure,subplot(2,1,1),plot(real(ydata)), title('K�ֲ��Ӳ�ʱ���Σ�ʵ��');
subplot(2,1,2),plot(imag(ydata)),title('K�ֲ��Ӳ�ʱ���Σ��鲿');
%%%%%%%%������ܶȺ����Ĳ���%%%%%%%%%%%%%%%%%%
num=100;
maxdat=max(abs(ydata));
mindat=min(abs(ydata));
NN=hist(abs(ydata),num);
xpdf1=num*NN/((sum(NN))*(maxdat-mindat));
xaxis1=mindat:(maxdat-mindat)/num:maxdat-(maxdat-mindat)/num;
alpha=sqrt(std(xdata).^2./(2*vmuc));%std()��xdata��׼��
th_val=lognpdf(xaxis1,xpdf1);

% xpdf1=getnpdf(abs(xdata),num,maxdat,mindat);
xaxis1=mindat:(maxdat-mindat)/num:maxdat-(maxdat-mindat)/num;
th_val=2*((xaxis1/(2*alpha)).^vmuc).*besselk((vmuc-1),xaxis1/alpha)./(alpha*gamma(vmuc));
figure,plot(xaxis1,xpdf1);
hold,plot(xaxis1,th_val,':r'),title('�Ӳ����ȷֲ�'),xlabel('����'),
ylabel('�����ܶ�');
signal=ydata;
signal=signal-mean(ydata);
M=256;
psd_dat=pburg(real(signal),16,M,fr);
psd_dat=psd_dat/(max(psd_dat));
freqx=0:0.5*M;
freqx=freqx*fr/M;
psd_dat = fftshift(fft(signal,2048));
freqx=f;

figure,plot(freqx,abs(psd_dat));
title('�Ӳ�Ƶ��'),xlabel('Ƶ��/Hz'),
ylabel('�������ܶ�');
%%%%%%%% �����˹����%%%%%%%%%%%%
powerf=exp(-freqx.^2/(2*sigmaf.^2))*200;
hold;
plot(freqx,powerf,':r');