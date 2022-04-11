function [noise_sp] = shepinzaosheng(fs,Bj,fj,frame,Prj,Tr)
%========================================================================%
%���ܣ�������Ƶ���������ź�                                              %
%���룺
%   fs          ����Ƶ��   
%   Bj          �����źŴ���
%   fj          �����ź�����Ƶ��
%   frame       ������
%   Vj          ���ŷ���
%   Tr          �����ظ�����
%�����
%   noise_sp    ��Ƶ�����ź�
%========================================================================%
t=0:1/fs:frame*Tr-1/fs;N=length(t);
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



%��������
phi=2*pi*rand(1,1);
sj=uj.*exp(1i*(2*pi*fj*t+phi));


noise_sp=sj;
end

