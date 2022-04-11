% 2021-05-14 
% ������������

% 2022-02-18 
% ��Ъ����ת�� ֱ�� ת������ 
% ���������

% 2022-02-19
% ��Ъ����ת�� �ظ� ת������
% ���������

close all;
clear;
clc;
%% 
Tp = 20e-6;                  % ����
B = 10e6;                   % ����
gama = B/Tp;



C = 3e8;                    % ���٣�  
fs = 4*B;                  % �����ʣ� 
R_wave_gate = 10e3;         % �������ģ�
T_wave_gate = Tp*2;         % ����ʱ��
nrn = round(fs*T_wave_gate);% ���������������
Tstart = 2*R_wave_gate/C - nrn/2/fs;    % �����������ʼʱ�̣�
Tend = 2*R_wave_gate/C + nrn/2/fs - 1/fs;   % �����������ֹʱ�̣�
tnrn = Tstart:1/fs:Tend;    % ���������ʱ�̣�


%%-------Ŀ�����
tar_num = 1;                % Ŀ�������
% R_tar = [10e3 10e3+500 10e3-100];
R_tar = [8e3 10e3+500 10e3-100];
R_tar = R_tar(1:tar_num);   % Ŀ�꾶����루��ʼʱ�̣���
V_tar = [50 100 150 200];
V_tar = V_tar(1:tar_num);   % Ŀ�꾶���ٶȣ���ʼʱ�̣���Ŀ���˶�״̬����������ֱ���˶���

tao_tar = 2*R_tar/C;        % ��һ�����壬Ŀ���ʱ�ӣ�

%% ----Ŀ��ز�
S_tar = zeros(1,nrn);
for num = 1:tar_num
    win = (tnrn-tao_tar(1,num)-Tp/2)>-Tp/2 & (tnrn-tao_tar(1,num)-Tp/2)<Tp/2;
    S_tar = S_tar + win.*exp(1j*pi*gama*(tnrn-tao_tar(1,num)-Tp/2).^2);
end
%------�ο��źţ��״﷢���źţ���
S_ref = exp(1j*pi*gama*(-Tp/2:1/fs:(Tp/2-1/fs)).^2);  


%% ------��Ъ�����������ɣ�ֱ��ת��or�ظ�ת����ͨ��repeat_num������������ƣ�
% �ظ�ת�����ţ�����1��
jam_amp = 1;
sampl_time = 2e-6;     % ����ʱ����
sampl_period = 6e-6;   % �������ڣ�
ISRJ_signal_1 = ISRJ(S_tar,jam_amp,sampl_time,sampl_period,fs);

% figure,subplot(2,1,1),plot(real(S_tar));title('�״﷢���ź�');
% subplot(2,1,2),plot(real(ISRJ_signal_1));title('��Ъ�����ظ�ת�������ź�');

% ���Ļ�ͼ
[c, r] = size(S_tar);
x = linspace(0, r/fs, r) * 1e3;% ms
F = linspace(-fs/2, fs/2, r);% Hz
hfig = figure(1);
subplot(3, 1, 1);
plot(x*1e3, real(S_tar));
title('\fontsize{10.5}\fontname{����}�״﷢���ź�ʱ����');
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10.5);
grid on;
axis tight;
xlabel('\fontname{����}ʱ��\fontname{Times New Roman}/us', 'FontSize', 10.5);
ylabel('\fontname{����}����\fontname{Times New Roman}/V', 'FontSize', 10.5);

subplot(3, 1, 2);
plot(x*1e3, real(ISRJ_signal_1));
title('\fontsize{10.5}\fontname{����}��Ъ����ת�������ź�ʱ����');
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10.5);
grid on;
axis tight;
xlabel('\fontname{����}ʱ��\fontname{Times New Roman}/us', 'FontSize', 10.5);
ylabel('\fontname{����}����\fontname{Times New Roman}/V', 'FontSize', 10.5);

subplot(3, 1, 3);
plot(F/1e6, abs(fftshift(fft(ISRJ_signal_1))));
% plot(F/1e6, db(abs(fftshift(fft(ISRJ_signal_1))) / max(abs(fftshift(fft(ISRJ_signal_1))))));
title('\fontsize{10.5}\fontname{����}��Ъ����ת�������ź�Ƶ��');  
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10.5);
grid on;
% axis tight;
xlabel('\fontname{����}Ƶ��\fontname{Times New Roman}/MHz', 'FontSize', 10.5);
ylabel('\fontname{����}����\fontname{Times New Roman}/V', 'FontSize', 10.5);


% ����ͼƬ�����ʽ
% figWidth = 14;
% figHeight = 12.9;  % ��������14cm,8.6cm���ƽ����Ϊ0.618
% set(hfig, 'PaperUnits', 'centimeters');  %centimeters %inches
% set(hfig, 'PaperPosition', [0 0 figWidth figHeight]);
% print(hfig, ['��Ъ����ת������.', 'tif'], '-r600', '-dtiff');   