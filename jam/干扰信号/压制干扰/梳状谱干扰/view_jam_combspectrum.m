function [ ] = view_jam_combspectrum( sig_noise,t_noise,fs )
%VIEW_JAM_COMBSPECTRUM 此处显示有关此函数的摘要
%   此处显示详细说明
h = waitbar(.33,'干扰产生中，请稍等！');
pause(.5)

N=length(t_noise);
figure;
subplot(311);
plot(t_noise,real(sig_noise));grid;xlabel('时间(s)');ylabel('幅度(V)');title('梳状谱噪声信号时域波形');axis([0 0.3e-4 -5 5]);
subplot(312);
Fn=fftshift(fft(sig_noise));K=N;k=floor(-K/2+0.5:K/2-0.5);dfs=fs/K;
plot(k*dfs,abs(Fn));grid;xlabel('频率(Hz)');ylabel('幅度(V/Hz)');title('梳状谱噪声信号幅度谱');
subplot(313);
window=hamming(100);noverlap=20;n=N;
yrj=pwelch(sig_noise,window,noverlap,n,fs);
Pdb=10*log10(yrj);
frj=0:fs/length(Pdb):(length(Pdb)-1)*(fs/length(Pdb));
 plot(frj-fs/2,fftshift(Pdb));grid;xlabel('频率(Hz)');ylabel('功率谱(dB/Hz)');title('梳状谱噪声信号功率谱密度');
 
waitbar(.67,h,'干扰产生中，请稍等！');
pause(1)
 
 close(h);
end

