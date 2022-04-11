function [noise_sp] = view_shepinzaosheng(fs,Bj,fj,frame,Prj,tr)
%========================================================================%
%���ܣ�������Ƶ���������ź�                                              %
%���룺
%   fs          ����Ƶ��   
%   Bj          �����źŴ���
%   fj          �����ź�����Ƶ��
%   frame       ������
%   Vj          ���ŷ���
%   tr          �����ظ�����
%�����
%   noise_sp    ��Ƶ�����ź�
%========================================================================%
h = waitbar(0, '���Ų����У����Եȣ�');
pause(.5)
t=0:1/fs:frame*tr-1/fs;N=length(t);
vRand=randn(1,2*N);
u=vRand(1:N)+1i*vRand(N+1:2*N);
Bn=Bj;

%Ƶ���ȡ
Fu=fftshift(fft(u));
n1=floor((fs-Bn)*N/2/fs);n2=floor((fs+Bn)*N/2/fs);
Fuj=[zeros(1,n1-1),Fu(n1:n2),zeros(1,N-n2)];
%��Ƶ��ʱ��
uj=ifft(ifftshift(Fuj));
Vj=dBmtoVwithoutr(Prj);
uj=PowerWeighted(Vj,uj);
waitbar(.33,h, '���Ų����У����Եȣ�');
pause(.1)

figure(1);
subplot(311);
plot(t,real(u));grid;xlabel('ʱ��(s)');ylabel('����(V)');
title('��������ʱ����');

subplot(312);
% [y,x]=hist(real(u),100);y=y./N./mean(diff(x));bar(x,y,1);grid;xlabel('����(V)');ylabel('�����ܶ�');
% title('�������������ܶȺ���');
[y,x]=ksdensity(real(u));bar(x,y);grid;xlabel('����(V)');ylabel('�����ܶ�');title('�������������ܶȺ���');
%������

window=hamming(100);noverlap=20;n=length(t);
yrj=pwelch(u,window,noverlap,n,fs);
Pdb=10*log10(yrj);
frj=0:fs/length(Pdb):(length(Pdb)-1)*(fs/length(Pdb));
subplot(313)
plot(frj-fs/2,fftshift(Pdb));grid;xlabel('Ƶ��(Hz)');ylabel('������(dBW/Hz)');title('���������������ܶ�');
waitbar(.67,h, '���Ų����У����Եȣ�');
pause(1)

%��������
phi=2*pi*rand(1,1);
sj=uj.*exp(1i*(2*pi*fj*t+phi));
figure(2);
subplot(311);
plot(t,real(sj));grid;xlabel('ʱ��(s)');ylabel('����(V)');
title('��Ƶ�����ź�ʱ����');
subplot(312);
% [y,x]=hist(real(sj),100);y=y./N./mean(diff(x));bar(x,y,1);grid;xlabel('����(V)');ylabel('�����ܶ�');
% title('��Ƶ�����źŸ����ܶȺ���');
[y,x]=ksdensity(real(sj));bar(x,y);grid;xlabel('����(V)');ylabel('�����ܶ�');title('��Ƶ�����źŸ����ܶȺ���');
subplot(313);
window=hamming(100);noverlap=20;n=N;
yrj=pwelch(sj,window,noverlap,n,fs);
Pdb=10*log10(yrj);
frj=0:fs/length(Pdb):(length(Pdb)-1)*(fs/length(Pdb));
 plot(frj-fs/2,fftshift(Pdb));grid;xlabel('Ƶ��(Hz)');ylabel('������(dBW/Hz)');title('��Ƶ�����źŹ������ܶ�');

noise_sp=sj;
close(h);
end

