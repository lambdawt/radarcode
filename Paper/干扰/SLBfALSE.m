close all;
clear all;
clc; 
Tp = 20e-6;                  % ����
B = 10e6;                   % ����
gama = B/Tp;                % ��Ƶб��
C = 3e8;                    % ���٣�  
fs = 4*B;                  % �����ʣ� 
R_wave_gate = 10e3;         % �������ģ�
T_wave_gate = Tp*5;         % ����ʱ��
nrn = round(fs*T_wave_gate);% ���������������
Tstart = 2*R_wave_gate/C - nrn/2/fs;    % �����������ʼʱ�̣�
Tend = 2*R_wave_gate/C + nrn/2/fs - 1/fs;   % �����������ֹʱ�̣�
tnrn = Tstart:1/fs:Tend;    % ���������ʱ�̣�


%%-------Ŀ�����
tar_num = 2;                % Ŀ�������
R_tar = [10e3 12e3];
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

figure;
plot(tnrn*1e6, abs(S_tar));
xlabel('\fontsize{10.5}\fontname{����}ʱ��\fontname{Times New Roman}/us', 'FontSize', 10.5);
ylabel('\fontsize{10.5}\fontname{����}����\fontname{Times New Roman}/V', 'FontSize', 10.5);


%%%%% �����ܼ����Ŀ�� %%%%%
% L=2e-6; % ��Ŀ����
% M=20;  % ��Ŀ�����
% k=T/L; % �ֶ���Ŀ
% duanNum=k+M-1;  % �ع����ܶ���
% T_g=duanNum*L;  % �ع���ʱ�䳤��
% gn=zeros(1,T_g*Fs);
% gn(1,1:N)=sn;
% for i=1:M-1
%     cur = gn(1,i*L*Fs+1:i*L*Fs+N);
%     gn(1,i*L*Fs+1:i*L*Fs+N)=cur+A*exp(1j*(2*pi*F0*(t)+pi*K*(t).^2+Phi0));
% end
% 
% %%%%% ���ܼ����Ŀ������źŽ�����ѹ %%%%%
% S_rec = gn;
% S_ref = A*exp(1j*(2*pi*F0*t+pi*K*t.^2+Phi0));
% temp = conv(S_rec,conj(fliplr(S_ref)));
% S_pc = temp(length(S_ref):end);
%                                      
% hfig = figure(1);
% subplot(3,1,1);
% plot((t+T/2)*1e6,real(sn));
% title('\fontsize{10.5}\fontname{����}���Ե�Ƶ�ź�ʵ��');
% set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10.5);
% grid on;
% axis tight;
% xlabel('\fontsize{10.5}\fontname{����}ʱ��\fontname{Times New Roman}/us', 'FontSize', 10.5);
% ylabel('\fontsize{10.5}\fontname{����}����\fontname{Times New Roman}/V', 'FontSize', 10.5)
% 
% subplot(3,1,2);
% plot(linspace(0,T_g,T_g*Fs)*1e6,abs(gn));
% title('\fontsize{10.5}\fontname{����}�ܼ����Ŀ������źŰ���');
% set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10.5);
% grid on;
% axis tight;
% xlabel('\fontsize{10.5}\fontname{����}ʱ��\fontname{Times New Roman}/us', 'FontSize', 10.5);
% ylabel('\fontsize{10.5}\fontname{����}����\fontname{Times New Roman}/V', 'FontSize', 10.5);
% 
% subplot(3,1,3);
% plot(linspace(0,T_g,T_g*Fs)*1e6,abs(S_pc));
% title('\fontsize{10.5}\fontname{����}�ܼ����Ŀ������ź���ѹ���');
% set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10.5);
% grid on;
% axis tight;
% xlabel('\fontsize{10.5}\fontname{����}ʱ��\fontname{Times New Roman}/us', 'FontSize', 10.5);
% ylabel('\fontsize{10.5}\fontname{����}����\fontname{Times New Roman}/V', 'FontSize', 10.5);