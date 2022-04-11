function [noise_tp] =view_zaoshengtiaopin(fs,Kfm,Pn,Prj,Bn,fj,frame,Tr)
%========================================================================%
%���ܣ�����������Ƶ�����ź�                                              %
%���룺
%   fs          ����Ƶ��  
%   Kfm         ��Ƶб��
%   Prj          ������Ƶ�źŹ���
%   Bj          ������Ƶ�źŵĵ�Ч���Ŵ���
%   Bn          ������������
%   fj          �����ź�����Ƶ��
%   frame       ������
%   Tr          �����ظ�����
%�����
%   noise_tp    ������Ƶ�ź�
%========================================================================%
% % Ԥ���������
% fs=160e6;%HZ
% Bj=20e6;%HZ
% fj=10e6;%HZ
% frame=64;
% Prj=40;%dBmW
% Tr=40e-6;%s
% Kfm=4e6;
% Bn=10e6;%HZ
% Pn=(Bj/2/(2.5*Kfm))^2;

t=0:1/fs:frame*Tr-1/fs;
N=length(t);
%��������������խ����˹��������
vGsRand=randn(1,2*N);
u=vGsRand(1:N)+1i*vGsRand(N+1:2*N);
Pn0=sqrt(fs*Pn/Bn);
u=Pn0*u;



%%%Ƶ���ȡ
Fu=fftshift(fft(u));
n1=floor((fs-Bn)*N/2/fs);n2=floor((fs+Bn)*N/2/fs);%
Fuj=[zeros(1,n1-1),Fu(n1:n2),zeros(1,N-n2)];
uj=ifft(ifftshift(Fuj));
figure(1);
subplot(311);
plot(t,real(uj));grid;xlabel('ʱ��(s)');ylabel('����(V)');
title('��������ʱ����');
subplot(312);
[y,x]=ksdensity(real(uj));bar(x,y);grid;xlabel('����(V)');ylabel('�����ܶ�');
title('�������������ܶȺ���');
subplot(313);
window=hamming(100);noverlap=20;n=N;
yrj=pwelch(uj,window,noverlap,n,fs);
Pdb=10*log10(yrj);
frj=0:fs/length(Pdb):(length(Pdb)-1)*(fs/length(Pdb));
 plot(frj-fs/2,fftshift(Pdb));grid;xlabel('Ƶ��(Hz)');ylabel('������(dB/Hz)');title('���������������ܶ�');
 
% waitbar(.67,h,'���Ų����У����Եȣ�');
% pause(1)

%Ƶ�ʵ���
phi=2*pi*rand(1,1);
%%%������������
Uj=zeros(1,N);
for i=1:N
    if i==1
        Uj(i)=real(uj(i))/fs;
    else
    Uj(i)=Uj(i-1)+real(uj(i))/fs; 
    end
end
sj=exp(1i*(2*pi*fj*t+2*pi*Kfm*Uj+phi));
Vj=dBmtoV(Prj);
sj=PowerWeighted(Vj,sj);
figure(2);
subplot(311);
plot(t,real(sj));grid;xlabel('ʱ��(s)');ylabel('����(V)');
title('������Ƶ�ź�ʱ����');
subplot(312);
[y,x]=ksdensity(real(sj));bar(x,y);grid;xlabel('����(V)');ylabel('�����ܶ�');
title('������Ƶ�źŸ����ܶȺ���');
subplot(313);
window=hamming(100);noverlap=20;n=N;
yrj=pwelch(sj,window,noverlap,n,fs);
Pdb=10*log10(yrj);
frj=0:fs/length(Pdb):(length(Pdb)-1)*(fs/length(Pdb));
 plot(frj-fs/2,fftshift(Pdb));grid;xlabel('Ƶ��(hz)');ylabel('������(db/hz)');title('������Ƶ�������ܶ�');
 noise_tp=sj;
%  close(h);
end

