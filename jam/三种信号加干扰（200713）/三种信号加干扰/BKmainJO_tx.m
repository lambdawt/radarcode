%主函数
function []=BKmainJO_tx(fz,B,Gt,Gr,F,Pt,L,Pfa,R,v,rcsk,sigma0,Te,tau,tr,f0,f1,fs,frame,num_jilei,num_tongdao,num_cankao,num_baohu,model,number1,code) 
% clear all;
%  close all;
% clc;
global c;
c=3e8;
%%%%%%%%%1雷达参数%%%%%%%%%%%%

%%%%%%%1.1发射机参数%%%%%%%%
% Gt=30;%雷达发射天线增益
% fz=10e9;%雷达中心频率
% Rmax=18e3; %干扰最大距离

lamta=c/fz;%0.03;%c/fz;%0.03;%波长

% load data/data_DeceptionJammingParameter;
load data/data_BFParameter;
load data/data_target0Parameter;
 load data/data_tx
%  Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
%   Aj=sqrt(Prj);
%   f_doppler1=2*v1/lamta;
sigma=rcs(rcsk,sigma0);%%%%%%%%%%%%%%%%%%%%%%%%%%目标的散射面积，注意这里有求目标散射面积的函数
f_doppler=2*v/lamta;%真目标多普勒频率

%%%%%%%%%%1.2脉冲参数%%%%%%%%%
%%%%%%%%%%1.2.1基本参数%%%%%
% tau=1e-6;%脉冲宽度
% tr=40e-6; %脉冲重复周期PRI

fr=1/tr;%脉冲重复频率
%%%%%%%%%1.2.2信号选择（相位编码）参数%%%%%%%
%model=2;
% % % model=1;
% % % number1=13;
% % % number2=2;
% % % %code=[0,1,2,3,0,1,2,3,0,1,2,3,0];
% % % code=[1,1,1,1,1,-1,-1,1,1,-1,1,-1,1];
% % % %%%%%%%%%%%2目标参数%%%%%%%%%%%
% % % %%%%%%%%%%%2.1基本目标参数%%%%
% % % R=2e3;    %真目标距离
% % % v=300;     %真目标速度
% % % sigma0=2;%目标的平均截面积
% % % rcsk=0;%采用不起伏或swerling4种模型，取值范围（0~2）

sigma=rcs(rcsk,sigma0);%%%%%%%%%%%%%%%%%%%%%%%%%%目标的散射面积，注意这里有求目标散射面积的函数
f_doppler=2*v/lamta;%真目标多普勒频率



%%%%%%%%%%%1.3接收机参数%%%%%%%%%%%
%%%%%%%%%%%1.3.1噪声损耗参数%%%%%%%
% % F=4;%噪声系数
% % L=3;%综合损耗
% % Te=1143;%温度
% % Gr=30;%雷达等效接收天线增益，单位（dB）
% % 
% % 
% % 
% % %%%%%%%%%%1.3.2接收机性能参数%%%%%%
% % B=20e6;%频带宽度
% % f1=10e6;  %本振频率


Srmin=-114+F+10*log10(B/1e6)+13;                                                                   %雷达接收机灵敏度
% Pt=((Rmax)^4)*(((4*pi)^3)*(10^(Srmin/10)))/(sigma*(lamta^2)*(10^((Gt/10)*2))*1e3);%雷达峰值功率
An=10*log10((1.382e-23)*Te*B*10^(F/10));                                                           %噪声强度

%%%%%%%%%%%2.2目标回波参数%%%%%%%%%%
Prs=((Pt*(10^((Gt/10)))*(10^((Gr/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L/10))); %目标回波信号功率
A=sqrt(Prs);                                                                                       %回波信号幅度



%%%%%%1.4仿真AD参数%%%%%%
% f0=20e6; %仿真中心频率
% fs=160e6; %采样频率

ts=1/fs; %采样间隔
tm=0:1/fs:tr-1/fs;%一个脉冲重复周期采样序列
N=length(tm);%一个脉冲重复周期采样点数长度

%%%%%%%1.5信号检测参数%%%%%%%%
% % frame=64;%仿真脉冲个数
% % num_jilei=32;%脉冲积累个数
% % num_tongdao=32;%相干积累时 fft点数
% % num_cankao=32;%恒虚警参考单元
% % num_baohu=5;%恒虚警保护单元
% % Pfa=1e-6;%恒虚警率

%%%%%%%%%%%信号级仿真%%%%%%%%%%
%%%%%%%%%%%%%1信号生成%%%%%%%%%%%%%

%%%%%%%%%%1.1波形生成：%%%%%%%%

%%%%%%%%%1.1.1相位编码信号%%%%%%%%
[y,y1,D]=shengchengBKxinhao(tau,fs,f0,model,number1,code,Pt,tr,ts);
figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(单位：秒)'), ylabel('y(单位：伏)'),title('相位编码发射信号');
figure,plot((0:fs/length(abs(fft(abs(fftshift(y))))):fs-fs/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏/Hz)'),title('相位编码发射信号的频谱');
% [y,A,Aj,An,ts,f_doppler,f_doppler1,O,N,fr]=shengchengLFMxinhao(B,Gt,Gr,lamta,F,Srmin,Pt,L,Pfa,R,v,R1,v1,sigma,Rmax,Rmin,Rf0,Gj,Gjr,rj,Kj,Te,Pj,tau,tr,f0,f1,fs,time_jam,frame,num_jilei,num_tongdao,num_cankao,num_baohu,number,deltf,deltt,Rj);
% y_fft_result=(fft(y));
%figure,plot(0:ts:(tau/ts-1)*ts,real(y)),xlabel('t(单位：秒)'), ylabel('y(单位：伏)'),title('发射信号');
%figure,plot((0:fs/length(abs(fft(abs(fftshift(y_fft_result))))):fs-fs/length(abs(fft(abs(fftshift(y_fft_result)))))),abs((y_fft_result))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('发射信号的频谱');

%%%%%%%%%%%1.1.2生成理想脉冲压缩系数
[M,match_filter_fft]=maiyaxishu(f0,fs,y/sqrt(Pt),tr,ts,N);

%%%%%%%%%%%%%%1.2发射信号%%%%%%%%%%%%%
% [s_fashe]=fashexinhao(y,tr,ts,Pt,frame);
% figure,plot(0:ts:N*ts,real(s_fashe)),xlabel('t(单位：秒)'), ylabel('y(单位：伏)'),title('发射信号');
% %xlim([0 4e-5])
% figure,plot((0:fs/length(abs(fft(abs(fftshift(s_fashe))))):fs-fs/length(abs(fft(abs(fftshift(s_fashe)))))),abs(fft(s_fashe))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('发射信号的频谱');

%%%%%%%%%%%%%1.3生成回波信号%%%%%%%%%%%%
[s_echo_2,echo]=BKhuiboxinhao(y1,R,tr,ts,A,N,frame,fs,f_doppler,tau); 
%生成噪声
[s_noise]=zaosheng(frame,N,An,B,fs);
%目标回波信号、假目标信号、噪声叠加在一起送入接收机
%%%%干扰信号  (D,y,R1,tr,ts,c,Aj,N,frame,fs,f_doppler1,tm,f0,tau,congmubiao,y1)
[noise_tx] =zaoshengtiaoxiang(fs_tx,Bj_tx,fj_tx,Prj_tx,Tr_tx,frame_tx);     
%目标回波信号、假目标信号、噪声叠加在一起送入接收机
s_echo_1=s_echo_2+s_noise+noise_tx;
t=0:1/fs:frame*tr-1/fs; 
s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau/2,tau);

%生成回波信号
figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('雷达接收信号');

%高放
[s_echo_1]=gaofang(f0,B,fs,s_echo_1);
%混频
[s_echo_1,f0]=hunpin(s_echo_1,N,frame,f1,fs,f0);
%检波及数字化
[s_echo_mf]=jianbo(s_echo_1,N,frame,f0,fs);
%脉冲压缩及降采样
[pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame,match_filter_fft,tau,D,ts);
figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(单位：s)'), ylabel('y(单位：dB)'),title('相位编码匹配滤波');
% s_pc_result=reshape(pc_result,1,M1*frame);
%动目标检测
[s_mtd]=mtd(pc_result1.',M1,num_jilei,num_tongdao);
figure,mesh(1:fr/num_tongdao:fr,0:ts*(tau/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau/D/ts)-ts*(tau/D/ts))*c/2,abs(s_mtd)),xlabel('多普勒频率，单位：Hz'),ylabel('距离，单位：米'),zlabel('y(单位：伏)'),title('相位编码MTD结果');
%恒虚警处理
hengxujing(M1,Pfa,s_mtd,num_cankao,num_tongdao,num_baohu,ts,c,tau,D);
end
