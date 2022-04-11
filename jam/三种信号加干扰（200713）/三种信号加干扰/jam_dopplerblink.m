function [ sig_jam,t_jam ] = jam_dopplerblink( fd,Td,R0,echo,fs,Pj,flagT )
%========================================================================%
%功能：产生多普勒闪烁欺骗干扰                                              %
%输入：
%   fd          多普勒闪烁基准频率   
%   Td          多普勒闪烁周期
%   R0          随机数组
%   echo        雷达照射信号
%   fs       	采样频率
%   Pj          干扰功率
%   flagT       干扰信号起始时间    

%输出：
%   sig_jam     欺骗干扰信号
%   t_jam       欺骗干扰信号时间轴
%========================================================================%
fd1=[fd,5*fd,10*fd,30*fd,50*fd];%多普勒闪烁频率
Aj=sqrt(Pj);                    %干扰信号幅度
Ntd=floor(flagT/(2*Td));        %0-flagT之间包含了多少个2Td
index=floor(5*R0(Ntd+1)+0.5)+1;     %R0(Ntd)是指每过一个周期换一个随机数从而得到一个多普勒频率
                                %乘以5是因为R0是在0-1间的随机数，5*R0变为0-5之间的随机数
curtime=flagT-Ntd*2*Td;         %flagT在0-2Td中的位置
N=length(echo);
t=(0:N-1)/fs;
if index==6
    index=index-1;
end

% 在0-Td上闪烁多普勒频率为正的，在Td-2Td上闪烁多普勒频率为负的
if curtime<=Td
    sig_jam=Aj*echo.*exp(1i*(2*pi*fd1(index)*t));
else 
    sig_jam=Aj*echo.*exp(1i*(2*pi*(-1)*fd1(index)*t));
end
t_jam=t;
end

