function [ sig_jam,t_jam ] = jam_dopplerblink( fd,Td,R0,echo,fs,Pj,flagT )
%========================================================================%
%���ܣ�������������˸��ƭ����                                              %
%���룺
%   fd          ��������˸��׼Ƶ��   
%   Td          ��������˸����
%   R0          �������
%   echo        �״������ź�
%   fs       	����Ƶ��
%   Pj          ���Ź���
%   flagT       �����ź���ʼʱ��    

%�����
%   sig_jam     ��ƭ�����ź�
%   t_jam       ��ƭ�����ź�ʱ����
%========================================================================%
fd1=[fd,5*fd,10*fd,30*fd,50*fd];%��������˸Ƶ��
Aj=sqrt(Pj);                    %�����źŷ���
Ntd=floor(flagT/(2*Td));        %0-flagT֮������˶��ٸ�2Td
index=floor(5*R0(Ntd+1)+0.5)+1;     %R0(Ntd)��ָÿ��һ�����ڻ�һ��������Ӷ��õ�һ��������Ƶ��
                                %����5����ΪR0����0-1����������5*R0��Ϊ0-5֮��������
curtime=flagT-Ntd*2*Td;         %flagT��0-2Td�е�λ��
N=length(echo);
t=(0:N-1)/fs;
if index==6
    index=index-1;
end

% ��0-Td����˸������Ƶ��Ϊ���ģ���Td-2Td����˸������Ƶ��Ϊ����
if curtime<=Td
    sig_jam=Aj*echo.*exp(1i*(2*pi*fd1(index)*t));
else 
    sig_jam=Aj*echo.*exp(1i*(2*pi*(-1)*fd1(index)*t));
end
t_jam=t;
end

