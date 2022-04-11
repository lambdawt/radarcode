function [s_mtd]=mtd(pc_result,M1,num_jilei,num_tongdao)
%%  输入
%   pc_result     脉冲压缩后信号
%   M1            脉压后信号总点数
%   num_jilei     脉冲积累数
%   num_tongdao   通道数
%%  输出
%   s_mtd         动目标检测后信号
%%  功能         动目标检测
%%
s_pc=pc_result;
for i=1:M1
    s_temp=s_pc(i,1:num_jilei);
       win=hamming(length(s_temp));
       s_temp=s_temp.*win';
    s_mtd(i,:)=fft(s_temp,num_tongdao);
end