function [noise_tx] =view_zaoshengtiaoxiang(fs,Bj,fj,Prj,Tr,frame)
%========================================================================%
%���ܣ�����������������ź�                                              %
%���룺
%   fs          ����Ƶ��  
%   Bj          ���Ŵ���
%   fj          �����ź�����Ƶ��
%   Prj         ���Ź��ʣ���λdBm��
%   Tr          �����ظ�����
%   frame       ������
%�����
%   noise_tx    ���������ź�
%========================================================================%
% % Ԥ���������
% fs=160e6;
% Bj=20e6;
% fj=10e6;
% frame=64;
% Prj=40;%dBmW
% Tr=40e-6;


h = waitbar(0,'���Ų����У����Եȣ�');
pause(.5)

t=0:1/fs:frame*Tr-1/fs;
N=length(t);
%������������
u=2*pi*rand(1,N)-pi;
waitbar(.33,h,'���Ų����У����Եȣ�');
pause(1)

figure(1);
subplot(311);
plot(t,real(u));grid;xlabel('ʱ��(s)');ylabel('����(V)');
title('��������ʱ����');

subplot(312);
[y,x]=ksdensity(u);bar(x,y);grid;xlabel('����(V)');ylabel('�����ܶ�');
title('�������������ܶȺ���');

subplot(313);
window=hamming(100);noverlap=20;n=N;
yrj=pwelch(u,window,noverlap,n,fs);
Pdb=10*log10(yrj);
frj=0:fs/length(Pdb):(length(Pdb)-1)*(fs/length(Pdb));
plot(frj-fs/2,fftshift(Pdb));grid;xlabel('Ƶ��(Hz)');ylabel('������(dB/Hz)');title('���������������ܶ�');
 
waitbar(.67,h,'���Ų����У����Եȣ�');
pause(1)
 %��������
 phi=2*pi*rand(1,1);
 sj=zeros(1,N);
 for i=1:N
 sj(i)=exp(1i*(2*pi*fj*t(i)+u(floor(t(i)*Bj)+1)+phi));
 end
 Vj=dBmtoV(Prj);
 sj=PowerWeighted(Vj,sj);
figure(2);
subplot(311);
plot(t,real(sj));grid;xlabel('ʱ��(s)');ylabel('����(V)');
title('���������ź�ʱ����');

subplot(312);
[y,x]=ksdensity(real(sj));bar(x,y);grid;xlabel('����(V)');ylabel('�����ܶ�');
title('���������źŸ����ܶȺ���');

subplot(313);
window=hamming(100);noverlap=20;n=N;
yrj=pwelch(sj,window,noverlap,n,fs);
Pdb=10*log10(yrj);
frj=0:fs/length(Pdb):(length(Pdb)-1)*(fs/length(Pdb));
 plot(frj-fs/2,fftshift(Pdb));grid;xlabel('Ƶ��(hz)');ylabel('������(db/hz)');title('�������๦�����ܶ�');

 noise_tx=sj;
  close(h);
end

