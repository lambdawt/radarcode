clear all；close all;
azi_num=2000;   %取2000个点
fr=1000;        %雷达重复频率

lamda0=0.05;   %杂波波长
sigmav=1.0;     %杂波方差
sigmaf=2*sigmav/lamda0;  


rand('state',sum(100*clock)); %产生服从U(0,1)分布的随机序列
d1=rand(1,azi_num);            
rand('state',7*sum(100*clock)+3);
d2=rand(1,azi_num);
xi=2*sqrt(-2*log(d1)).*cos(2*pi*d2);  %正交且独立的高斯序列N(0，1)
xq=2*sqrt(-2*log(d1)).*sin(2*pi*d2);
%形成滤波器频率响应
coe_num=12;           %求滤波器系数，用傅里叶级数展开法
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
%生成高斯谱杂波
xxi=conv(b,xi);   
xxq=conv(b,xq);   
xxi=xxi(coe_num*2+1:azi_num+coe_num*2);%目的是去掉暂态响应
xxq=xxq(coe_num*2+1:azi_num+coe_num*2);
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
% 
% clear;
% azi_num=2000;
% sigmac=1.2 ; 
% rand('state',sum(100*clock)); %产生服从U(0,1)分布的随机序列
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
% yyq=(d2-xqmuc)/xqsigmac;    %归一化
% 
% 
% ydata = yyi + 1j* yyq;
% 
% % ydata = randn(1,2000) + 1j* randn(1,2000);
% 
% num=100;                   %求概率密度函数的参数
% maxdat=max(abs(ydata));
% mindat=min(abs(ydata));
% NN=hist(abs(ydata),num);   
% xpdf1=num*NN/((sum(NN))*(maxdat-mindat));        %用直方图估计的概率密度函数
% xaxisl=mindat:(maxdat-mindat)/num:maxdat-(maxdat-mindat)/num;  
% th_val=(xaxisl./sigmac.^2).*exp(-xaxisl.^2./(2*sigmac.^2));   %概率密度函数理论值
% figure;
% plot(xaxisl,xpdf1);               %做出仿真结果的概率密度函数曲线
% hold on;plot(xaxisl,th_val,'r:'); %做出理论概率密度函数曲线
% title('杂波幅度分布');
% xlabel('幅度');ylabel('概率密度');

signal=ydata;
signal=signal-mean(signal);      %求功率谱密度，先去掉直流分量
figure;M=256;                     %用burg法估计功率谱密度
psd_dat=pburg(real(signal),32,M,fr);
psd_dat=psd_dat/(max(psd_dat));   %归一化处理
freqx=0:0.5*M;
freqx=freqx*fr/M;
plot(freqx,psd_dat);title('杂波频谱');
xlabel('频率/HZ');ylabel('功率谱密度');

%做出理想高斯谱曲线
powerf=exp(-freqx.^2/(2*sigmaf.^2));
hold on;plot(freqx,powerf,'r:');