%主函数
function []=LFMtuoyin(fz,B,B1,Gt,Gr,lamta,F,Srmin,Pt,L,Pfa,R,v,rcsk,sigma0,Rmax,Te,tau,tr,f0,f1,fs,frame,num_jilei,num_tongdao,num_cankao,num_baohu,R1,v1,Pj,Gj,Gjr,Rj,congmubiao) 
clear all;
% close all;
clc;
global c;
c=3e8;
%1雷达参数
%1.1发射机参数
Gt=30;%雷达发射天线增益
fz=10e9;%雷达中心频率
Rmax=18e3; %干扰最大距离

lamta=c/fz;%波长

%1.2脉冲参数
%1.2.1基本参数
tau=1e-6;%脉冲宽度
tr=40e-6; %脉冲重复周期PRI
B1=20e6;%线性调频信号带宽

k=B1/tau;%线性调频信号调制系数
fr=1/tr;%脉冲重复频率

%2目标参数
%2.1基本目标参数
R=2e3;    %真目标距离
v=300;     %真目标速度
sigma0=2;%目标的平均截面积
rcsk=0;%采用不起伏或swerling4种模型，取值范围（0~2）
%2.2干扰目标参数
R1=R;
v1=v;
tf=7;%距离拖引
vf=7;%速度拖引
Rj=2e3;
Gj=30;%干扰机发射天线增益
Gjr=10;%雷达接收天线等效增益
Pj=9.2533;%干扰机峰值功率
congmubiao=10;%从目标个数

sigma=rcs(rcsk,sigma0);%%%%%%%%%%%%%%%%%%%%%%%%%%目标的散射面积，注意这里有求目标散射面积的函数
f_doppler=2*v/lamta;%真目标多普勒频率

%1.3接收机参数
%1.3.1噪声
F=4;%噪声系数
L=3;%综合损耗
Te=1143;%温度
Gr=30;%雷达接收天线等效增益

%1.3.2接收机性能参数
B=20e6;%接收机带宽
f1=10e6;  %本振频率

Srmin=-114+F+10*log10(B/1e6)+13;%雷达接收机灵敏度
Pt=((Rmax)^4)*(((4*pi)^3)*(10^(Srmin/10)))/(sigma*(lamta^2)*(10^((Gt/10)*2))*1e3);%雷达峰值功率
An=10*log10((1.382e-23)*Te*B*10^(F/10));%噪声强度

%1.4仿真参数



%2.2目标回波参数
Prs=((Pt*(10^((Gt/10)))*(10^((Gr/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L/10))); %目标回波信号功率
A=sqrt(Prs);%回波信号幅度
%%干扰回波参数

Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L/10)));
Aj=sqrt(Prj);

%速度
f_doppler1=f_doppler;
%密集数
%角度

%仿真参数
f0=20e6; %仿真频率
fs=160e6; %采样频率
ts=1/fs;%采样间隔
tm=0:1/fs:tr-1/fs;%一个脉冲重复周期采样序列
O=tau/ts;%一个脉宽采样点数长度
N=length(tm);%一个脉冲重复周期采样点数长度

% 1.5信号检测参数
frame=129;%仿真脉冲个数
num_jilei=16;%脉冲积累个数
num_tongdao=16;%相干积累时 fft点数
num_cankao=16;%恒虚警参考单元
num_baohu=5;%恒虚警保护单元
Pfa=1e-6;%恒虚警率

%生成发射信号
[y,D]=shengchengLFMxinhao(B1,Pt,tau,f0,tm,k);
figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(单位：秒)'), ylabel('y(单位：伏)'),title('发射信号');
figure,plot((0:fs/length(abs(fft(abs(fftshift(y))))):fs-fs/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('发射信号的频谱');
%生成理想脉冲压缩系数
[M,match_filter_fft]=maiyaxishu(f0,fs,y,tr,ts,N);

%生成回波信号
[s_echo_2,echo]=LFMtuoyinhuiboxinhao(y,R,tr,ts,c,A,N,frame,fs,f_doppler,tm,f0,B1,tau,k); 
%figure,plot(0:ts:(N-1)*ts,real(s_echo_2(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('回波信号');
%figure,plot((0:fs/length(abs(fft(abs(fftshift(echo(1,:)))))):fs-fs/length(abs(fft(abs(fftshift(echo(1,:))))))),abs((fft(echo(1,:))))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('回波信号的频谱');

%生成干扰信号
[s_ft,echo3]=LFMtuoyinganrao(R,frame,tf,Aj,tm,tau,f0,B1,N,f_doppler1,f_doppler,vf,lamta,k,ts,c,fs);

%生成噪声
[s_noise]=zaosheng(frame,N,An,B,fs);
%目标回波信号、假目标信号、噪声叠加在一起送入接收机
s_echo_1=s_echo_2+s_noise+s_ft;
figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('雷达接收信号');
figure,plot((0:fs/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs-fs/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('雷达接收信号的频谱');

%高放
% [s_echo_1]=gaofang(f0,B,fs,s_echo_1);
%混频
[s_echo_1,f0]=hunpin(s_echo_1,N,frame,f1,fs,f0);
%检波及数字化
[s_echo_mf]=jianbo(s_echo_1,N,frame,f0,fs);
%脉冲压缩及降采样
[pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame,match_filter_fft,tau,D,ts);
figure,plot(0:ts*c/2:(M-1)*ts*c/2,20*log10(abs(pc_result1(1,:)))),xlabel('t(单位：s)'), ylabel('y(单位：dB)'),title('匹配滤波');
s_pc_result=reshape(pc_result,1,M1*frame);
%动目标检测
[s_mtd,s_mtd1]=MTD_tuoyin(pc_result,M1,num_jilei,num_tongdao);
% figure,mesh(1:fr/num_tongdao*lamta/2:fr*lamta/2,0:ts*(tau/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau/D/ts)-ts*(tau/D/ts))*c/2,abs(s_mtd)),xlabel('多普勒频率，单位：Hz'),ylabel('距离，单位：米'),zlabel('y(单位：伏)'),title('MTD结果');
subplot(121),mesh(1:fr/num_tongdao*lamta/2:fr*lamta/2,0:ts*(tau/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau/D/ts)-ts*(tau/D/ts))*c/2,abs(s_mtd)),xlabel('速度，单位：米/秒'),ylabel('距离，单位：米'),zlabel('y(单位：伏)'),title('停拖期（1至16个脉冲）的MTD结果');
subplot(122),mesh(1:fr/num_tongdao*lamta/2:fr*lamta/2,0:ts*(tau/D/ts)*c/2:(length(abs(s_mtd1(:,1)))*ts*(tau/D/ts)-ts*(tau/D/ts))*c/2,abs(s_mtd1)),xlabel('速度，单位：米/秒'),ylabel('距离，单位：米'),zlabel('y(单位：伏)'),title('拖引期（113至128个脉冲）的MTD结果');
%恒虚警处理
% hengxujing(M1,Pfa,s_mtd,num_cankao,num_tongdao,num_baohu,ts,c,tau,D);
hengxujing_tuoyin(M1,Pfa,s_mtd,num_cankao,num_tongdao,num_baohu,ts,c,tau,D,s_mtd1);
end