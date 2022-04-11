function [ ] = view_SmartNoise( vRandSig,fs,tau,R)
%VIEW_RANDSIGNAL �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
c=3e8;
N=length(vRandSig);
t=0:1/fs:(N-1)/fs;

h = waitbar(.33,'���Ų����У����Եȣ�');
pause(.5)

figure;
subplot(411);
plot(t,real(vRandSig));grid;xlabel('ʱ��(s)');ylabel('����(V)');title('��������ʱ����');axis([0 0.3e-4 -5 5]);

subplot(412);
nStart=ceil((2*R/c)*fs);
nEnd=nStart+ceil(tau*fs);
[y,x]=ksdensity(real(vRandSig(nStart:nEnd)));bar(x,y);grid;
xlabel('����(V)');ylabel('�����ܶ�');title('��������һ�������ڵĸ����ܶȺ���');

subplot(413);
[y,x]=ksdensity(real(vRandSig));bar(x,y);grid;
xlabel('����(V)');ylabel('�����ܶ�');title('����ʱ����������������ܶȺ���');

subplot(414);
window=hamming(100);noverlap=20;n=N;
vEstimationSpectrum=pwelch(vRandSig,window,noverlap,n,fs);
Pdb=10*log10(vEstimationSpectrum);
fx=0:fs/length(Pdb):(length(Pdb)-1)*(fs/length(Pdb));
 plot(fx-fs/2,fftshift(Pdb));grid;xlabel('Ƶ��(Hz)');ylabel('������(dBW/Hz)');title('���������������ܶ�');
 
 waitbar(.67,h,'���Ų����У����Եȣ�');
pause(1)

close(h);
end

