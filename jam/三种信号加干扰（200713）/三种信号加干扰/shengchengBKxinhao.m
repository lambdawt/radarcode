function [y,y1,D]=shengchengBKxinhao(tau,fs,f0,~,number1,code,Pt,tr,ts)
%% 输入 
%   脉宽       tau
%   采样频率   fs
%   中心频率   f0
%   model模式 1巴克码 2自定义
%   number1   输入编码总长度（巴克码模式）
%   code      输入对应数组2相或4相码序列（自定模式）
%   Pt        雷达发射功率
%   tr        脉冲重复周期
%   ts        采样周期
%% 输出
%   y         相位编码信号（PRI内）     
%   y1        相位编码信号（原始产生信号）
%   D         脉压系数
%% 功能 生成相位编码信号
%   可以选择巴克码或自定义
%   可以选择二相或四相码
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rmin=sqrt((Kj*Pt*sigma)/(4*pi*(Pj)*rj))
% f1=-f0+f1;
%Prs=((Pt*(10^((Gt/10)))*(10^((Gr/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L/10))); %目标回波信号功率
%Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L/10)));%干扰信号功率
%A=sqrt(Prs);%回波信号幅度
%Aj=sqrt(Prj);%干扰信号幅度
%An=10*log10((1.382e-23)*Te*B*10^(F/10));%噪声强度
% deltt1=2*(R1-R)/c;%假目标相对真目标延迟时间
%fr=1/tr;%脉冲重复频率
%ts=1/fs; 
%f_doppler=2*v/lamta;%真目标多普勒频率
%f_doppler1=2*v1/lamta;%假目标多普勒频率
%tm=0:1/fs:tr-1/fs;%一个脉冲重复周期采样序列
%O=tau/ts;%一个脉宽采样点数
%ft=linspace(-tau/2,tau/2,O);%一个脉宽采样序列
%N=length(tm);%一个脉冲重复周期采样点数长度
%k=B/tau;%B/fs*2*pi/max(ft);  %调制系数
%y=exp(j*2*pi*((f0+k*ft/2).*ft));%产生和脉宽长度相同的线性调频信号

%number1=length码元个数  ,number2= 2相位 [1,-1][0,pi]   4相位 [0,2*pi]/[0,1,2,3]phase[0,pi/2,pi,3pi/2]
%code=[编码……]
%0巴克码or自定义 
%0.1巴克码  0.1.1选择位数7位 0.1.2选择位数13位 
%0.2自定义  0.2.1选择二相码 0.2.2 选择四相码  0.2.*.1请输入码元个数 0.2.*.2请输入编码       
%model=2;
%number2=4;
% % if model==1
% %     code=[1,1,1,-1,-1,-1,1];
% % else if model==2
% %         code=[1,1,1,1,1,-1,-1,1,1,-1,1,-1,1];%13位巴克码
% %     else
% %         code=code;
% %     end
% % end
        

% code=[1,1,1,1,1,-1,-1,1,1,-1,1,-1,1];%13位巴克码%code=
%code=[1,1,1,-1,-1,-1,1];%7位巴克码
%code=[0,1,2,3,0,1,2,3,0,1,2,3,0];%任意指定输入文本
D=number1;
ncode=length(code);%number1
tao=tau/ncode; %number1=13
t_tao=0:1/fs:tao-1/fs;
Ntao=length(t_tao);
pha=0;
t=0:1/fs:ncode*tao-1/fs;%number1=13
y1=zeros(1,length(t));
for i=1:ncode
   % if code(i)==1
        pha=code(i)*pi/2;
    %else pha=0;
    %end
    y1(1,(i-1)*Ntao+1:i*Ntao)=exp(1j*(2*pi*f0*t_tao+pha));
end
y=sqrt(Pt).*[y1,zeros(1,fix((tr-tau)/ts))];


%figure(1),subplot(2,1,1),
%figure,plot(real(y)),xlabel('t(单位:s)'),title('调制信号');
% s_fft_result=abs(fft(y(1:Ntao)));
%subplot(2,1,2),
%figure,plot((0:fs/Ntao:fs/2-fs/Ntao),abs(s_fft_result(1:Ntao/2))),xlabel('频率(单位:Hz)'),title('码内信号频谱');