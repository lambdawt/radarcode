function [ sig_jam,t_jam ] = jam_AGC( CurrentT,Pj,Period,D,radio,echo,fs )
%========================================================================%
%功能：AGC干扰                                              %
%输入：   
%   CurrentT       	当前时间
%   Pj              干扰功率
%   Period          通断周期
%   D               通断工作比
%   radio           幅度分配比例
%   echo            回波信号
%   fs              采样频率
%输出：
%   sig_jam         AGC干扰信号
%   t_jam           输出信号的时间轴
%========================================================================%
currt=mod(CurrentT,Period);
Am=radio*sqrt(Pj);
Amm=(1-radio)*sqrt(Pj);
N=length(echo);
if currt>Period*D
    sj=Am*echo;
else
    sj=Amm*echo;
end

sig_jam=sj;
t_jam=(0:N-1)/fs;

end

