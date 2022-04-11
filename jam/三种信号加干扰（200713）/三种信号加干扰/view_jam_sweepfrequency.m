function [ ] = view_jam_sweepfrequency( sig_noise,t_noise,fs )
%VIEW_JAM_COMBSPECTRUM �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
h = waitbar(.33,'���Ų����У����Եȣ�');
pause(.5)

N=length(t_noise);
figure;
subplot(311);
plot(t_noise,real(sig_noise));grid;xlabel('ʱ��(s)');ylabel('����(V)');title('ɨƵ�����ź�ʱ����');axis([0 0.3e-4 -5 5]);
subplot(312);
window=hamming(128);noverlap=64;n=N;
yrj=pwelch(sig_noise,window,noverlap,n,fs);
Pdb=10*log10(yrj);
frj=0:fs/length(Pdb):(length(Pdb)-1)*(fs/length(Pdb));
plot(frj-fs/2,fftshift(Pdb));grid;xlabel('Ƶ��(Hz)');ylabel('������(dB/Hz)');title('ɨƵ�����źŹ������ܶ�');
%Fn=fftshift(fft(sig_noise));K=N;k=floor(-K/2+0.5:K/2-0.5);dfs=fs/K;
%plot(k*dfs,abs(Fn));grid;xlabel('Ƶ��(Hz)');ylabel('����(V)');title('ɨƵ�����źŷ�����');
subplot(313);
window=hamming(128);noverlap=64;n=N;
[~,F,T,P]=spectrogram(sig_noise,window,noverlap,128,fs);
surf(T,F,10*log10(abs(P)),'EdgeColor','none');axis tight;
view(0,90);
grid;xlabel('ʱ�䣨s��');ylabel('Ƶ��(Hz)');

waitbar(.67,h,'���Ų����У����Եȣ�');
pause(1)

close(h);
%yrj=pwelch(sig_noise,window,noverlap,n,fs);
%Pdb=10*log10(yrj);
%frj=0:fs/length(Pdb):(length(Pdb)-1)*(fs/length(Pdb));
% plot(frj-fs/2,fftshift(Pdb));grid;xlabel('Ƶ��(Hz)');ylabel('������(dB/Hz)');title('ɨƵ�����źŹ������ܶ�');
end

