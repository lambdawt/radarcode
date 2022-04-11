function [ ] = view_jam_passive( sig_noise,t_noise,fs)
%VIEW_JAM_COMBSPECTRUM �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
N=length(t_noise);
figure;
subplot(211);
plot(t_noise,real(sig_noise));grid;xlabel('ʱ��(s)');ylabel('����(V)');title('���������ź�ʱ����');
subplot(212);
Fn=fftshift(fft(sig_noise));K=N;k=floor(-K/2+0.5:K/2-0.5);dfs=fs/K;
plot(k*dfs,abs(Fn));grid;xlabel('Ƶ��(Hz)');ylabel('����(V/Hz)');title('���������źŷ�����')

end