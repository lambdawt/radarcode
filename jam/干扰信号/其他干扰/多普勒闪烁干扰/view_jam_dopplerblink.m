function [ ] = view_jam_dopplerblink( echo,sig_jam,t_jam,fs )
%VIEW_JAM_COMBSPECTRUM �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
N=length(t_jam);
figure;
% subplot(211);
% plot(t_jam,sig_jam);grid;xlabel('ʱ��(s)');ylabel('����(V)');title('��������˸��ƭ�����ź�ʱ����');
Fecho=fftshift(fft(echo));
Fn=fftshift(fft(sig_jam));
K=N;k=floor(-K/2+0.5:K/2-0.5);dfs=fs/K;
plot(k*dfs,abs(Fecho),k*dfs,abs(Fn));grid;legend('�״��ź�','�����ź�');xlabel('Ƶ��(Hz)');ylabel('����(V/Hz)');title('�ܵ���������˸��ƭ���ŵ��״�ز��ź�Ƶ��');

%plot(k*dfs,abs(Fn));grid;xlabel('Ƶ��(hz)');ylabel('����(V)');title('��������˸��ƭ�����źŷ�����');
end

