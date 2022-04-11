%简单脉冲主函数
function []=JDmainJ_dopplerblink(fz,B,Gt,Gr,F,Pt,L,Pfa,R,v,rcsk,sigma0,Te,tau,tr,f0,f1,fs,frame,num_jilei,num_tongdao,num_cankao,num_baohu) 
% clear all;
% close all;
% clc;
global c;
c=3e8;
%1雷达参数
%1.1发射机参数
% % % Gt=30;%雷达发射天线增益
% % % fz=10e9;%雷达中心频率
% % % Rmax=18e3; %干扰最大距离

lamta=c/fz;%波长
load data/data_dopplerblink;
load data/data_JDParameter;
load data/data_target0Parameter;
R0=rand(100);
sigma=rcs(rcsk,sigma0);%%%%%%%%%%%%%%%%%%%%%%%%%%目标的散射面积，注意这里有求目标散射面积的函数
f_doppler=2*v/lamta;%真目标多普勒频率
%1.2脉冲参数
%1.2.1基本参数
% % tau=1e-6;%脉冲宽度
% % tr=40e-6; %脉冲重复周期PRI
% % % B1=20e6;%线性调频信号带宽

% k=B1/tau;%线性调频信号调制系数
fr=1/tr;%脉冲重复频率
% D=B1*tau;%脉压比
%2目标参数
%2.1基本目标参数
% % % R=2e3;    %真目标距离
% % % v=300;     %真目标速度
% % % sigma0=2;%目标的平均截面积
% % % rcsk=0;%采用不起伏或swerling4种模型，取值范围（0~2）

sigma=rcs(rcsk,sigma0);%%%%%%%%%%%%%%%%%%%%%%%%%%目标的散射面积，注意这里有求目标散射面积的函数
f_doppler=2*v/lamta;%真目标多普勒频率

%1.3接收机参数
%1.3.1噪声
% % % F=4;%噪声系数
% % % L=3;%综合损耗
% % % Te=1143;%温度
% % % Gr=30;%雷达接收天线等效增益
% % % 
% % % %1.3.2接收机性能参数
% % % B=20e6;%接收机带宽
% % % f1=10e6;  %本振频率

Srmin=-114+F+10*log10(B/1e6)+13;%雷达接收机灵敏度
% % Pt=((Rmax)^4)*(((4*pi)^3)*(10^(Srmin/10)))/(sigma*(lamta^2)*(10^((Gt/10)*2))*1e3);%雷达峰值功率
An=10*log10((1.382e-23)*Te*B*10^(F/10));%噪声强度

%1.4仿真参数



%2.2目标回波参数
Prs=((Pt*(10^((Gt/10)))*(10^((Gr/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L/10))); %目标回波信号功率
A=sqrt(Prs);%回波信号幅度

%仿真参数
% % f0=20e6; %仿真频率
% % fs=160e6; %采样频率
ts=1/fs;%采样间隔
tm=0:1/fs:tr-1/fs;%一个脉冲重复周期采样序列
N=length(tm);%一个脉冲重复周期采样点数长度

% 1.5信号检测参数
% % % frame=64;%仿真脉冲个数
% % % num_jilei=32;%脉冲积累个数
% % % num_tongdao=32;%相干积累时 fft点数
% % % num_cankao=16;%恒虚警参考单元
% % % num_baohu=5;%恒虚警保护单元
% % % Pfa=1e-6;%恒虚警率

%生成发射信号
[y,D]=shengchengJDxinhao(Pt,tau,f0,tm);
figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(单位：秒)'), ylabel('y(单位：伏)'),title('简单脉冲发射信号');
figure,plot((0:fs/length(abs(fft(abs(fftshift(y))))):fs-fs/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏/Hz)'),title('简单脉冲发射信号的频谱');
%生成理想脉冲压缩系数
[M,match_filter_fft]=maiyaxishu(f0,fs,y/sqrt(Pt),tr,ts,N);

%生成回波信号
[s_echo_2,echo]=JDhuiboxinhao(R,c,A,N,frame,fs,f_doppler,tm,f0,tau); 
%生成噪声
[s_noise]=zaosheng(frame,N,An,B,fs);
%目标回波信号、假目标信号、噪声叠加在一起送入接收机
%%%%干扰信号 D,R1,tr,c,Aj,N,frame,fs,f_doppler1,tm,f0,tau,congmubiao
[ sig_jam,t_jam ] = jam_dopplerblink( fd_dopplerblink,Td_dopplerblink,R0,s_echo_2,fs_dopplerblink,Pj_dopplerblink,flagT_dopplerblink );     
%目标回波信号、假目标信号、噪声叠加在一起送入接收机
s_echo_1=s_echo_2+s_noise+sig_jam;
t=0:1/fs:frame*tr-1/fs; 
s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau/2,tau);
%目标回波信号、假目标信号、噪声叠加在一起送入接收机

%加干扰回波信号
figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('雷达接收信号');

%高放
[s_echo_1]=gaofang(f0,B,fs,s_echo_1);
%混频
[s_echo_1,f0]=hunpin(s_echo_1,N,frame,f1,fs,f0);
%检波及数字化
[s_echo_mf]=jianbo(s_echo_1,N,frame,f0,fs);
%脉冲压缩及降采样
[pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame,match_filter_fft,tau,D,ts);
figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(单位：s)'), ylabel('y(单位：dB)'),title('简单脉冲匹配滤波');
% s_pc_result=reshape(pc_result,1,M1*frame);
%动目标检测
[s_mtd]=mtd(pc_result1.',M1,num_jilei,num_tongdao);
figure,mesh(1:fr/num_tongdao:fr,0:ts*(tau/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau/D/ts)-ts*(tau/D/ts))*c/2,abs(s_mtd)),xlabel('多普勒频率，单位：Hz'),ylabel('距离，单位：米'),zlabel('y(单位：伏)'),title('简单脉冲MTD结果');
%恒虚警处理
hengxujing(M1,Pfa,s_mtd,num_cankao,num_tongdao,num_baohu,ts,c,tau,D);
end