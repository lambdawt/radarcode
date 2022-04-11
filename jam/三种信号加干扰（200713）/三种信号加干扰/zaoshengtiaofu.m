function [noise_tf] =zaoshengtiaofu(fs,Bj,fj,Prj,Tr,frame)
%========================================================================%
%���ܣ������������������ź�                                              %
%���룺
%   fs          ����Ƶ��  
%   Bj          ���Ŵ���
%   fj          �����ź�����Ƶ��
%   Prj         ���Ź���
%   Tr          �����ظ�����
%   frame       ������
%�����
%   noise_tf    ������Ƶ�ź�
%========================================================================%
t=0:1/fs:frame*Tr-1/fs;
N=length(t);
Bn=Bj;
vGsRandnum=randn(1,2*N);
u=vGsRandnum(1:N)+1i*vGsRandnum(N+1:2*N);
Prj=dBmtoW(Prj);
Pn=Prj*2/19;
Pn0=sqrt(fs*Pn/Bj);
u=Pn0*u;
%%%Ƶ���ȡ
Fu=fftshift(fft(u));
n1=floor((fs-Bn)*N/2/fs);n2=floor((fs+Bn)*N/2/fs);%
Fuj=[zeros(1,n1-1),Fu(n1:n2),zeros(1,N-n2)];
uj=ifft(ifftshift(Fuj));

%���ȵ���
 phi=2*pi*rand(1,1);
 U0=sqrt(Prj-Pn/2);
 sj=(U0+real(uj)).*exp(1i*(2*pi*fj*t+phi));

 noise_tf=sj;
end

