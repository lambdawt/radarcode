function [ Volt ] = dBmtoVwithoutR( Vdbm )
%DBMTOMV ����(dBmW)ת��Ϊ��ѹ(V)
%   �����ǵ���
Volt=sqrt(10^((Vdbm-30)/10));
end

