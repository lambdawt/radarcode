function [noise_tx] =zaoshengtiaoxiang(fs,Bj,fj,Prj,Tr,frame)
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
t=0:1/fs:frame*Tr-1/fs;
N=length(t);
%������������
u=2*pi*rand(1,N)-pi;


 %��������
 phi=2*pi*rand(1,1);
 sj=zeros(1,N);
 for i=1:N
 sj(i)=exp(1i*(2*pi*fj*t(i)+u(floor(t(i)*Bj)+1)+phi));
 end
 Vj=dBmtoV(Prj);
 sj=PowerWeighted(Vj,sj);

 noise_tx=sj;
end

