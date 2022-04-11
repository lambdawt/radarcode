function [ vSmartNoiseSig ] = jam_smartnoise( vRadarSig,Pn,Prj,Bn,Kfm,fs )
%JAM_SMARTNOISE ������������
%   �˴���ʾ��ϸ˵��
vRadarSig=vRadarSig/max(vRadarSig);
N=length(vRadarSig);
t=0:1/fs:(N-1)/fs;
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
sj=vRadarSig.*exp(1i*(2*pi*Kfm*Uj+phi));
Vj=dBmtoV(Prj);
sj=Vj*sj;

vSmartNoiseSig=sj;

end

