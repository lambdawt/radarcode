%%���Ե�Ƶ�ź�
T=10e-6;                                  %p�������ʱ��10us
B=30e6;                                   %���Ե�Ƶ�źŵ�Ƶ�����30MHz
K=B/T;                                      %��Ƶб��
Fs=2*B;Ts=1/Fs;                      %����Ƶ�ʺͲ������
N=T/Ts;
t=linspace(-T/2,T/2,N);
St=exp(1i*pi*K*t.^2);                    %���Ե�Ƶ�ź�

hfig = figure;
subplot(211)
plot(t*1e6,real(St));
xlabel('\fontname{����}ʱ��\fontname{Times new roman}/us', 'FontSize',10.5);
ylabel('\fontname{����}����\fontname{Times new roman}/V', 'FontSize',10.5);
title('\fontname{����}���Ե�Ƶ�źŵ�ʵ��', 'FontSize',10.5);
set(gca,'FontName','Times New Roman','FontSize',10.5);
grid on;
axis tight;

subplot(212)
freq=linspace(-Fs/2,Fs/2,N);
plot(freq*1e-6,fftshift(abs(fft(St))));
xlabel('\fontname{����}Ƶ��\fontname{Times new roman}/MHz', 'FontSize',10.5);
ylabel('\fontname{����}����\fontname{Times new roman}/V', 'FontSize',10.5);
title('\fontname{����}���Ե�Ƶ�źŵķ�Ƶ����', 'FontSize',10.5);
set(gca,'FontName','Times New Roman','FontSize',10.5);
grid on;
axis tight;


% ����ͼƬ��С��Ӧ��ͨ������õ���Ⱥ͸߶ȡ�
% ����A4ֽ�Ŀ��Ϊ21cm���������ҳ�߾��Ϊ2.5cm�����ĵĿ�Ⱦ���16cm��
% ���ҳ���ǵ�������Ϊ���14cm���߶�8.6cm�ȽϺ��ʣ��ƽ����0.618����
% ˫���Ļ�����7cm����4.3cm���ɡ�
figWidth = 12;  % ����ͼƬ���
figHeight = 7.416;  % ����ͼƬ�߶�
set(hfig,'PaperUnits','centimeters'); % ͼƬ�ߴ����õ�λ��inches��Ӣ�磬centimeters������
set(hfig,'PaperPosition',[0 0 figWidth figHeight]);
fileout = 'test1.'; % ���ͼƬ���ļ���
print(hfig,[fileout,'tif'],'-r600','-dtiff'); % ����ͼƬ��ʽ���ֱ���
