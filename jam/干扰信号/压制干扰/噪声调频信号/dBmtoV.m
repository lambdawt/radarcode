function [ Volt ] = dBmtoV( Vdbm )
%DBMTOMV ����(dBmW)ת��Ϊ��ѹ(V)
%   �Ƚ�������dBmW�����W
%   �ٸ��ݹ�ʽP=U^2/R ��õ�ѹU
%   ����RΪ50ŷ
Volt=sqrt(50*10^((Vdbm-30)/10));
end

