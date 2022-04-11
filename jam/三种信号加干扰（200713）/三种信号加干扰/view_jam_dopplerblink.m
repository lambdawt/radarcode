function [ ] = view_jam_dopplerblink( echo,sig_jam,t_jam,fs )
%VIEW_JAM_COMBSPECTRUM 此处显示有关此函数的摘要
%   此处显示详细说明
N=length(t_jam);
figure;
% subplot(211);
% plot(t_jam,sig_jam);grid;xlabel('时间(s)');ylabel('幅度(V)');title('多普勒闪烁欺骗干扰信号时域波形');
Fecho=fftshift(fft(echo));
Fn=fftshift(fft(sig_jam));
K=N;k=floor(-K/2+0.5:K/2-0.5);dfs=fs/K;
plot(k*dfs,abs(Fecho),k*dfs,abs(Fn));grid;legend('雷达信号','干扰信号');xlabel('频率(Hz)');ylabel('幅度(V/Hz)');title('受到多普勒闪烁欺骗干扰的雷达回波信号频谱');

%plot(k*dfs,abs(Fn));grid;xlabel('频率(hz)');ylabel('幅度(V)');title('多普勒闪烁欺骗干扰信号幅度谱');
end

