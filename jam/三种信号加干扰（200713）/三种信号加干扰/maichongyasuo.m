function [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame,match_filter_fft,tau,D,ts)
%%  输入
%   s_echo_mf 回波信号
%   M         脉冲采样点数（一个PRI内）
%   frame     仿真脉冲数
%   match_filter_fft    匹配滤波系数
%   tau       脉宽
%   D         脉冲压缩比
%   ts        采样周期
%%  输出
%   pc_result    脉压结果
%   pc_result1   降采样结果
%   M1           脉压后信号总点数
%%  功能   匹配滤波
%%
for i=1:frame
    s_echo_fft_result(:,i)=fft(s_echo_mf(:,i),M);
    pc_result_fft(:,i)=match_filter_fft.*s_echo_fft_result(:,i);
    pc_result1(:,i)=ifft(pc_result_fft(:,i),M);
    pc_result(:,i)=downsample(ifft(pc_result_fft(:,i),M),floor(tau/D/ts),floor(tau/D/ts/2));
end
M1=length(pc_result1(:,i));
pc_result1=pc_result1.';
% M=N;
% figure,plot(0:ts*c/2:(length(abs(pc_result(:,1)))*ts-ts)*c/2,(pc_result(:,1))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('第一个回波的脉冲压缩');
% M1=length(pc_result(:,i));%原始

% s_pc_result=reshape(pc_result,1,length(pc_result(:,i))*frame);%将矩阵转为1行M*frame列
% figure,subplot(221),plot(0:ts:(frame*M-1)*ts,abs(s_echo_1)),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('真目标、假目标、噪声');axis([-inf inf 0 max(abs(s_pc_result))])
% subplot(222),plot(0:ts:length(abs(s_pc_result))*ts-ts,abs(s_pc_result)),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('脉冲压缩处理后的包络');