function [ ] = view_SmartNoise( vRandSig,fs,tau,R)
%VIEW_RANDSIGNAL 此处显示有关此函数的摘要
%   此处显示详细说明
c=3e8;
N=length(vRandSig);
t=0:1/fs:(N-1)/fs;

h = waitbar(.33,'干扰产生中，请稍等！');
pause(.5)

figure;
subplot(411);
plot(t,real(vRandSig));grid;xlabel('时间(s)');ylabel('幅度(V)');title('灵巧噪声时域波形');axis([0 0.3e-4 -5 5]);

subplot(412);
nStart=ceil((2*R/c)*fs);
nEnd=nStart+ceil(tau*fs);
[y,x]=ksdensity(real(vRandSig(nStart:nEnd)));bar(x,y);grid;
xlabel('幅度(V)');ylabel('概率密度');title('灵巧噪声一个脉宽内的概率密度函数');

subplot(413);
[y,x]=ksdensity(real(vRandSig));bar(x,y);grid;
xlabel('幅度(V)');ylabel('概率密度');title('整个时域的灵巧噪声概率密度函数');

subplot(414);
window=hamming(100);noverlap=20;n=N;
vEstimationSpectrum=pwelch(vRandSig,window,noverlap,n,fs);
Pdb=10*log10(vEstimationSpectrum);
fx=0:fs/length(Pdb):(length(Pdb)-1)*(fs/length(Pdb));
 plot(fx-fs/2,fftshift(Pdb));grid;xlabel('频率(Hz)');ylabel('功率谱(dBW/Hz)');title('灵巧噪声功率谱密度');
 
 waitbar(.67,h,'干扰产生中，请稍等！');
pause(1)

close(h);
end

