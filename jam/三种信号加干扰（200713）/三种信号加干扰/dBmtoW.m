function [ Pw ] = dBmtoW( Pdbm )
%DBMTOV ��λת��
%   ��dBmWת��ΪW
Pw=(10^(Pdbm/10))/1e3;
end

