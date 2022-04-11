function [ sig_jam,t_jam ] = jam_AGC( CurrentT,Pj,Period,D,radio,echo,fs )
%========================================================================%
%���ܣ�AGC����                                              %
%���룺   
%   CurrentT       	��ǰʱ��
%   Pj              ���Ź���
%   Period          ͨ������
%   D               ͨ�Ϲ�����
%   radio           ���ȷ������
%   echo            �ز��ź�
%   fs              ����Ƶ��
%�����
%   sig_jam         AGC�����ź�
%   t_jam           ����źŵ�ʱ����
%========================================================================%
currt=mod(CurrentT,Period);
Am=radio*sqrt(Pj);
Amm=(1-radio)*sqrt(Pj);
N=length(echo);
if currt>Period*D
    sj=Am*echo;
else
    sj=Amm*echo;
end

sig_jam=sj;
t_jam=(0:N-1)/fs;

end

