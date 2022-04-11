%������
function []=BKmainJO_tx(fz,B,Gt,Gr,F,Pt,L,Pfa,R,v,rcsk,sigma0,Te,tau,tr,f0,f1,fs,frame,num_jilei,num_tongdao,num_cankao,num_baohu,model,number1,code) 
% clear all;
%  close all;
% clc;
global c;
c=3e8;
%%%%%%%%%1�״����%%%%%%%%%%%%

%%%%%%%1.1���������%%%%%%%%
% Gt=30;%�״﷢����������
% fz=10e9;%�״�����Ƶ��
% Rmax=18e3; %����������

lamta=c/fz;%0.03;%c/fz;%0.03;%����

% load data/data_DeceptionJammingParameter;
load data/data_BFParameter;
load data/data_target0Parameter;
 load data/data_tx
%  Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
%   Aj=sqrt(Prj);
%   f_doppler1=2*v1/lamta;
sigma=rcs(rcsk,sigma0);%%%%%%%%%%%%%%%%%%%%%%%%%%Ŀ���ɢ�������ע����������Ŀ��ɢ������ĺ���
f_doppler=2*v/lamta;%��Ŀ�������Ƶ��

%%%%%%%%%%1.2�������%%%%%%%%%
%%%%%%%%%%1.2.1��������%%%%%
% tau=1e-6;%������
% tr=40e-6; %�����ظ�����PRI

fr=1/tr;%�����ظ�Ƶ��
%%%%%%%%%1.2.2�ź�ѡ����λ���룩����%%%%%%%
%model=2;
% % % model=1;
% % % number1=13;
% % % number2=2;
% % % %code=[0,1,2,3,0,1,2,3,0,1,2,3,0];
% % % code=[1,1,1,1,1,-1,-1,1,1,-1,1,-1,1];
% % % %%%%%%%%%%%2Ŀ�����%%%%%%%%%%%
% % % %%%%%%%%%%%2.1����Ŀ�����%%%%
% % % R=2e3;    %��Ŀ�����
% % % v=300;     %��Ŀ���ٶ�
% % % sigma0=2;%Ŀ���ƽ�������
% % % rcsk=0;%���ò������swerling4��ģ�ͣ�ȡֵ��Χ��0~2��

sigma=rcs(rcsk,sigma0);%%%%%%%%%%%%%%%%%%%%%%%%%%Ŀ���ɢ�������ע����������Ŀ��ɢ������ĺ���
f_doppler=2*v/lamta;%��Ŀ�������Ƶ��



%%%%%%%%%%%1.3���ջ�����%%%%%%%%%%%
%%%%%%%%%%%1.3.1������Ĳ���%%%%%%%
% % F=4;%����ϵ��
% % L=3;%�ۺ����
% % Te=1143;%�¶�
% % Gr=30;%�״��Ч�����������棬��λ��dB��
% % 
% % 
% % 
% % %%%%%%%%%%1.3.2���ջ����ܲ���%%%%%%
% % B=20e6;%Ƶ�����
% % f1=10e6;  %����Ƶ��


Srmin=-114+F+10*log10(B/1e6)+13;                                                                   %�״���ջ�������
% Pt=((Rmax)^4)*(((4*pi)^3)*(10^(Srmin/10)))/(sigma*(lamta^2)*(10^((Gt/10)*2))*1e3);%�״��ֵ����
An=10*log10((1.382e-23)*Te*B*10^(F/10));                                                           %����ǿ��

%%%%%%%%%%%2.2Ŀ��ز�����%%%%%%%%%%
Prs=((Pt*(10^((Gt/10)))*(10^((Gr/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L/10))); %Ŀ��ز��źŹ���
A=sqrt(Prs);                                                                                       %�ز��źŷ���



%%%%%%1.4����AD����%%%%%%
% f0=20e6; %��������Ƶ��
% fs=160e6; %����Ƶ��

ts=1/fs; %�������
tm=0:1/fs:tr-1/fs;%һ�������ظ����ڲ�������
N=length(tm);%һ�������ظ����ڲ�����������

%%%%%%%1.5�źż�����%%%%%%%%
% % frame=64;%�����������
% % num_jilei=32;%������۸���
% % num_tongdao=32;%��ɻ���ʱ fft����
% % num_cankao=32;%���龯�ο���Ԫ
% % num_baohu=5;%���龯������Ԫ
% % Pfa=1e-6;%���龯��

%%%%%%%%%%%�źż�����%%%%%%%%%%
%%%%%%%%%%%%%1�ź�����%%%%%%%%%%%%%

%%%%%%%%%%1.1�������ɣ�%%%%%%%%

%%%%%%%%%1.1.1��λ�����ź�%%%%%%%%
[y,y1,D]=shengchengBKxinhao(tau,fs,f0,model,number1,code,Pt,tr,ts);
figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(��λ����)'), ylabel('y(��λ����)'),title('��λ���뷢���ź�');
figure,plot((0:fs/length(abs(fft(abs(fftshift(y))))):fs-fs/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����/Hz)'),title('��λ���뷢���źŵ�Ƶ��');
% [y,A,Aj,An,ts,f_doppler,f_doppler1,O,N,fr]=shengchengLFMxinhao(B,Gt,Gr,lamta,F,Srmin,Pt,L,Pfa,R,v,R1,v1,sigma,Rmax,Rmin,Rf0,Gj,Gjr,rj,Kj,Te,Pj,tau,tr,f0,f1,fs,time_jam,frame,num_jilei,num_tongdao,num_cankao,num_baohu,number,deltf,deltt,Rj);
% y_fft_result=(fft(y));
%figure,plot(0:ts:(tau/ts-1)*ts,real(y)),xlabel('t(��λ����)'), ylabel('y(��λ����)'),title('�����ź�');
%figure,plot((0:fs/length(abs(fft(abs(fftshift(y_fft_result))))):fs-fs/length(abs(fft(abs(fftshift(y_fft_result)))))),abs((y_fft_result))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�����źŵ�Ƶ��');

%%%%%%%%%%%1.1.2������������ѹ��ϵ��
[M,match_filter_fft]=maiyaxishu(f0,fs,y/sqrt(Pt),tr,ts,N);

%%%%%%%%%%%%%%1.2�����ź�%%%%%%%%%%%%%
% [s_fashe]=fashexinhao(y,tr,ts,Pt,frame);
% figure,plot(0:ts:N*ts,real(s_fashe)),xlabel('t(��λ����)'), ylabel('y(��λ����)'),title('�����ź�');
% %xlim([0 4e-5])
% figure,plot((0:fs/length(abs(fft(abs(fftshift(s_fashe))))):fs-fs/length(abs(fft(abs(fftshift(s_fashe)))))),abs(fft(s_fashe))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�����źŵ�Ƶ��');

%%%%%%%%%%%%%1.3���ɻز��ź�%%%%%%%%%%%%
[s_echo_2,echo]=BKhuiboxinhao(y1,R,tr,ts,A,N,frame,fs,f_doppler,tau); 
%��������
[s_noise]=zaosheng(frame,N,An,B,fs);
%Ŀ��ز��źš���Ŀ���źš�����������һ��������ջ�
%%%%�����ź�  (D,y,R1,tr,ts,c,Aj,N,frame,fs,f_doppler1,tm,f0,tau,congmubiao,y1)
[noise_tx] =zaoshengtiaoxiang(fs_tx,Bj_tx,fj_tx,Prj_tx,Tr_tx,frame_tx);     
%Ŀ��ز��źš���Ŀ���źš�����������һ��������ջ�
s_echo_1=s_echo_2+s_noise+noise_tx;
t=0:1/fs:frame*tr-1/fs; 
s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau/2,tau);

%���ɻز��ź�
figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('�״�����ź�');

%�߷�
[s_echo_1]=gaofang(f0,B,fs,s_echo_1);
%��Ƶ
[s_echo_1,f0]=hunpin(s_echo_1,N,frame,f1,fs,f0);
%�첨�����ֻ�
[s_echo_mf]=jianbo(s_echo_1,N,frame,f0,fs);
%����ѹ����������
[pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame,match_filter_fft,tau,D,ts);
figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(��λ��s)'), ylabel('y(��λ��dB)'),title('��λ����ƥ���˲�');
% s_pc_result=reshape(pc_result,1,M1*frame);
%��Ŀ����
[s_mtd]=mtd(pc_result1.',M1,num_jilei,num_tongdao);
figure,mesh(1:fr/num_tongdao:fr,0:ts*(tau/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau/D/ts)-ts*(tau/D/ts))*c/2,abs(s_mtd)),xlabel('������Ƶ�ʣ���λ��Hz'),ylabel('���룬��λ����'),zlabel('y(��λ����)'),title('��λ����MTD���');
%���龯����
hengxujing(M1,Pfa,s_mtd,num_cankao,num_tongdao,num_baohu,ts,c,tau,D);
end
