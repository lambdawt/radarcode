function [ sig_noise,t_noise ] = jam_sweepfrequency( fs,Bj,fj,frame,Prj,Tr,T_fr,Time_begin,K_sweep )
%========================================================================%
%���ܣ�����ɨƵ�����ź�                                              %
%���룺
%   fs          ����Ƶ��   
%   Bj          �����źŴ���
%   fj          �����ź�����Ƶ��
%   frame       ������
%   Prj         ���Ź���
%   Tr          �����ظ�����
%   T_fr        ɨƵ����    �������ʱ���᳤�ȷ�Χ�ڣ��������Դﵽ���Ƶ�ʣ�
%   Time_begin  ���ſ�ʼʱ�䣨����ǰʱ�̣��źŵ���ʼʱ�䣬����Ϊ����ֵ��
%   K_sweep     ����Ƶ�ʱ仯б�ʣ�K_sweep*T_frӦ����fj��һ���������ϣ�����仯�����ԣ�
%�����
%   sig_noise    ɨƵ�����ź�
%   t_noise      ����źŵ�ʱ����
%========================================================================%

t=0:1/fs:frame*Tr-1/fs;N=length(t);
Bn=Bj;
u=sqrt(Prj)*randn(1,N);
%Ƶ���ȡ
Fu=fftshift(fft(u));
n1=floor((fs-Bn)*N/2/fs);n2=floor((fs+Bn)*N/2/fs);
Fuj=[zeros(1,n1-1),Fu(n1:n2),zeros(1,N-n2)];
%��Ƶ��ʱ��
uj=ifft(ifftshift(Fuj));

%��������
F_fr=1/T_fr;
t_begin=Time_begin-floor(Time_begin*F_fr)*T_fr;      %t_begin �ǵ�ǰʱ����һ��ɨƵ�����е�λ��
t_f=t_begin+(0:N-1)/fs;                              %��ɨƵ�����е�t_beginλ���𣬺���ʱ�䳤��N�ķ�Χ��������Ƶ�ʽ������Ա仯
tj=Time_begin+t;
sj=uj.*exp(1i*(2*pi*(fj+K_sweep*mod(t_f,T_fr)).*tj));%��mod������ʱ�����ǿ�ʼʱ���Ǵ�t_begin��ʼ�����ڼ��ϳ���N*ts��ʱ���п��ܳ���ɨƵ���ڣ���ʱҪ�ѳ����Ĳ����ۻش�ɨƵ���ڵ�0ʱ�̿�ʼ��
                                                     %ɨƵ����Сֵ��fj�����ֵ��:��ʱ�䳬��һ��ɨƵ���ڣ�fj+K_sweep*T_fr����ʱ��С��һ��ɨƵ���ڣ�;
sig_noise=sj;
t_noise=tj;
end

