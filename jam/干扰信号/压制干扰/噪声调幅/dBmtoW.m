function [ Pw ] = dBmtoW( Pdbm )
%DBMTOV 单位转换
%   将dBmW转换为W
Pw=(10^(Pdbm/10))/1e3;
end

