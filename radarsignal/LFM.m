%%线性调频信号
T=10e-6;                                  %p脉冲持续时间10us
B=30e6;                                   %线性调频信号的频带宽度30MHz
K=B/T;                                      %调频斜率
Fs=2*B;Ts=1/Fs;                      %采样频率和采样间隔
N=T/Ts;
t=linspace(-T/2,T/2,N);
St=exp(1i*pi*K*t.^2);                    %线性调频信号

hfig = figure;
subplot(211)
plot(t*1e6,real(St));
xlabel('\fontname{宋体}时间\fontname{Times new roman}/us', 'FontSize',10.5);
ylabel('\fontname{宋体}幅度\fontname{Times new roman}/V', 'FontSize',10.5);
title('\fontname{宋体}线性调频信号的实部', 'FontSize',10.5);
set(gca,'FontName','Times New Roman','FontSize',10.5);
grid on;
axis tight;

subplot(212)
freq=linspace(-Fs/2,Fs/2,N);
plot(freq*1e-6,fftshift(abs(fft(St))));
xlabel('\fontname{宋体}频率\fontname{Times new roman}/MHz', 'FontSize',10.5);
ylabel('\fontname{宋体}幅度\fontname{Times new roman}/V', 'FontSize',10.5);
title('\fontname{宋体}线性调频信号的幅频特性', 'FontSize',10.5);
set(gca,'FontName','Times New Roman','FontSize',10.5);
grid on;
axis tight;


% 关于图片大小，应该通过计算得到宽度和高度。
% 例如A4纸的宽度为21cm，如果左右页边距各为2.5cm，正文的宽度就是16cm。
% 如果页面是单栏，设为宽度14cm，高度8.6cm比较合适（黄金比例0.618）。
% 双栏的话，宽7cm，高4.3cm即可。
figWidth = 12;  % 设置图片宽度
figHeight = 7.416;  % 设置图片高度
set(hfig,'PaperUnits','centimeters'); % 图片尺寸所用单位，inches—英寸，centimeters—厘米
set(hfig,'PaperPosition',[0 0 figWidth figHeight]);
fileout = 'test1.'; % 输出图片的文件名
print(hfig,[fileout,'tif'],'-r600','-dtiff'); % 设置图片格式、分辨率
