function [ ] = view_jam_passive( sig_noise,t_noise,fs)
%VIEW_JAM_COMBSPECTRUM 此处显示有关此函数的摘要
%   此处显示详细说明
N=length(t_noise);
figure;
subplot(211);
plot(t_noise,real(sig_noise));grid;xlabel('时间(s)');ylabel('幅度(V)');title('箔条干扰信号时域波形');
subplot(212);
Fn=fftshift(fft(sig_noise));K=N;k=floor(-K/2+0.5:K/2-0.5);dfs=fs/K;
plot(k*dfs,abs(Fn));grid;xlabel('频率(Hz)');ylabel('幅度(V/Hz)');title('箔条干扰信号幅度谱')

end