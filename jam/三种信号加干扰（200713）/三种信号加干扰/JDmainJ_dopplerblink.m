%������������
function []=JDmainJ_dopplerblink(fz,B,Gt,Gr,F,Pt,L,Pfa,R,v,rcsk,sigma0,Te,tau,tr,f0,f1,fs,frame,num_jilei,num_tongdao,num_cankao,num_baohu) 
% clear all;
% close all;
% clc;
global c;
c=3e8;
%1�״����
%1.1���������
% % % Gt=30;%�״﷢����������
% % % fz=10e9;%�״�����Ƶ��
% % % Rmax=18e3; %����������

lamta=c/fz;%����
load data/data_dopplerblink;
load data/data_JDParameter;
load data/data_target0Parameter;
R0=rand(100);
sigma=rcs(rcsk,sigma0);%%%%%%%%%%%%%%%%%%%%%%%%%%Ŀ���ɢ�������ע����������Ŀ��ɢ������ĺ���
f_doppler=2*v/lamta;%��Ŀ�������Ƶ��
%1.2�������
%1.2.1��������
% % tau=1e-6;%������
% % tr=40e-6; %�����ظ�����PRI
% % % B1=20e6;%���Ե�Ƶ�źŴ���

% k=B1/tau;%���Ե�Ƶ�źŵ���ϵ��
fr=1/tr;%�����ظ�Ƶ��
% D=B1*tau;%��ѹ��
%2Ŀ�����
%2.1����Ŀ�����
% % % R=2e3;    %��Ŀ�����
% % % v=300;     %��Ŀ���ٶ�
% % % sigma0=2;%Ŀ���ƽ�������
% % % rcsk=0;%���ò������swerling4��ģ�ͣ�ȡֵ��Χ��0~2��

sigma=rcs(rcsk,sigma0);%%%%%%%%%%%%%%%%%%%%%%%%%%Ŀ���ɢ�������ע����������Ŀ��ɢ������ĺ���
f_doppler=2*v/lamta;%��Ŀ�������Ƶ��

%1.3���ջ�����
%1.3.1����
% % % F=4;%����ϵ��
% % % L=3;%�ۺ����
% % % Te=1143;%�¶�
% % % Gr=30;%�״�������ߵ�Ч����
% % % 
% % % %1.3.2���ջ����ܲ���
% % % B=20e6;%���ջ�����
% % % f1=10e6;  %����Ƶ��

Srmin=-114+F+10*log10(B/1e6)+13;%�״���ջ�������
% % Pt=((Rmax)^4)*(((4*pi)^3)*(10^(Srmin/10)))/(sigma*(lamta^2)*(10^((Gt/10)*2))*1e3);%�״��ֵ����
An=10*log10((1.382e-23)*Te*B*10^(F/10));%����ǿ��

%1.4�������



%2.2Ŀ��ز�����
Prs=((Pt*(10^((Gt/10)))*(10^((Gr/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L/10))); %Ŀ��ز��źŹ���
A=sqrt(Prs);%�ز��źŷ���

%�������
% % f0=20e6; %����Ƶ��
% % fs=160e6; %����Ƶ��
ts=1/fs;%�������
tm=0:1/fs:tr-1/fs;%һ�������ظ����ڲ�������
N=length(tm);%һ�������ظ����ڲ�����������

% 1.5�źż�����
% % % frame=64;%�����������
% % % num_jilei=32;%������۸���
% % % num_tongdao=32;%��ɻ���ʱ fft����
% % % num_cankao=16;%���龯�ο���Ԫ
% % % num_baohu=5;%���龯������Ԫ
% % % Pfa=1e-6;%���龯��

%���ɷ����ź�
[y,D]=shengchengJDxinhao(Pt,tau,f0,tm);
figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(��λ����)'), ylabel('y(��λ����)'),title('�����巢���ź�');
figure,plot((0:fs/length(abs(fft(abs(fftshift(y))))):fs-fs/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����/Hz)'),title('�����巢���źŵ�Ƶ��');
%������������ѹ��ϵ��
[M,match_filter_fft]=maiyaxishu(f0,fs,y/sqrt(Pt),tr,ts,N);

%���ɻز��ź�
[s_echo_2,echo]=JDhuiboxinhao(R,c,A,N,frame,fs,f_doppler,tm,f0,tau); 
%��������
[s_noise]=zaosheng(frame,N,An,B,fs);
%Ŀ��ز��źš���Ŀ���źš�����������һ��������ջ�
%%%%�����ź� D,R1,tr,c,Aj,N,frame,fs,f_doppler1,tm,f0,tau,congmubiao
[ sig_jam,t_jam ] = jam_dopplerblink( fd_dopplerblink,Td_dopplerblink,R0,s_echo_2,fs_dopplerblink,Pj_dopplerblink,flagT_dopplerblink );     
%Ŀ��ز��źš���Ŀ���źš�����������һ��������ջ�
s_echo_1=s_echo_2+s_noise+sig_jam;
t=0:1/fs:frame*tr-1/fs; 
s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau/2,tau);
%Ŀ��ز��źš���Ŀ���źš�����������һ��������ջ�

%�Ӹ��Żز��ź�
figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('�״�����ź�');

%�߷�
[s_echo_1]=gaofang(f0,B,fs,s_echo_1);
%��Ƶ
[s_echo_1,f0]=hunpin(s_echo_1,N,frame,f1,fs,f0);
%�첨�����ֻ�
[s_echo_mf]=jianbo(s_echo_1,N,frame,f0,fs);
%����ѹ����������
[pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame,match_filter_fft,tau,D,ts);
figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(��λ��s)'), ylabel('y(��λ��dB)'),title('������ƥ���˲�');
% s_pc_result=reshape(pc_result,1,M1*frame);
%��Ŀ����
[s_mtd]=mtd(pc_result1.',M1,num_jilei,num_tongdao);
figure,mesh(1:fr/num_tongdao:fr,0:ts*(tau/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau/D/ts)-ts*(tau/D/ts))*c/2,abs(s_mtd)),xlabel('������Ƶ�ʣ���λ��Hz'),ylabel('���룬��λ����'),zlabel('y(��λ����)'),title('������MTD���');
%���龯����
hengxujing(M1,Pfa,s_mtd,num_cankao,num_tongdao,num_baohu,ts,c,tau,D);
end