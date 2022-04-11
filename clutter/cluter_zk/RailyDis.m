clear all;close all;clc
azi_num=1920;   %取2000个点
fr=12e3;

% sigmav=40*0.03;     %杂波方差
sigmav = 40;
sigmaf=2*sigmav*6e9/3e8;

% fr=1000;
% lamda0=0.05;
% sigmav=0.1;
% sigmaf=2*sigmav/lamda0;

xi=randn(1,azi_num);  %正交且独立的高斯序列N(0，1)
xq=randn(1,azi_num);
%形成滤波器频率响应

FreXi = fftshift(fft(xi,2048));
f = (-length(FreXi)/2:1:length(FreXi)/2-1)/length(FreXi)*fr;
hf=exp(-f.^2/(2*sigmaf^2));% 高斯谱
% hf=1./(1+abs(f/sigmaf).^4);% 全极谱
% hf=exp((-abs(f)/sigmaf));%指数谱3

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
yyq=(xxq-xqmuc)/xqsigmac;    %归一化
sigmac=1.2 ;           %杂波的标准差
yyi=sigmac*yyi;        %使瑞利分布杂波具有指定的标准差
yyq=sigmac*yyq;        %使瑞利分布虚部杂波
ydata=yyi+j*yyq;       %瑞利分布杂波形成
figure;
subplot(211);plot(real(ydata)); %瑞利分布杂波实部
title('瑞利杂波时域波形，实部');  
subplot(212);plot(imag(ydata)); %瑞利分布杂波虚部
title('瑞利杂波时域波形，虚部')

num=100;                   %求概率密度函数的参数
maxdat=max(abs(ydata));
mindat=min(abs(ydata));
NN=hist(abs(ydata),num);   
xpdf1=num*NN/((sum(NN))*(maxdat-mindat));        %用直方图估计的概率密度函数
xaxisl=mindat:(maxdat-mindat)/num:maxdat-(maxdat-mindat)/num;  
th_val=(xaxisl./sigmac.^2).*exp(-xaxisl.^2./(2*sigmac.^2));   %概率密度函数理论值
figure;
plot(xaxisl,xpdf1);               %做出仿真结果的概率密度函数曲线
hold on;plot(xaxisl,th_val,'r:'); %做出理论概率密度函数曲线
title('杂波幅度分布');
xlabel('幅度');ylabel('概率密度');

signal=ydata;
signal=signal-mean(signal);      %求功率谱密度，先去掉直流分量
figure;M=256;                     %用burg法估计功率谱密度
% psd_dat=pburg(abs(signal),32,M,fr);
% psd_dat=psd_dat/(max(psd_dat));   %归一化处理
% freqx=0:0.5*M;
% freqx=freqx*fr/M;
psd_dat = fftshift(fft(signal,2048));
freqx=f;

plot(freqx,abs(psd_dat));title('杂波频谱');
xlabel('频率/HZ');ylabel('功率谱密度');

%做出理想高斯谱曲线
% powerf=exp(-freqx.^2/(2*sigmaf^2))*200;%高斯
% powerf=1./(1+abs(f/sigmaf).^3)*200;%全极谱
powerf=exp((-abs(f)/sigmaf))*200;%指数谱
v= freqx*3e8/2/6e9;
hold on;plot(freqx,powerf,'r:');

Acc=fftshift(fft(ydata));
v = (-length(Acc)/2:1:length(Acc)/2-1)/length(Acc)*fr*3e8/2/6e9;