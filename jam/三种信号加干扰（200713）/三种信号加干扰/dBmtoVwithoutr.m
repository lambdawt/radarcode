function [ Volt ] = dBmtoVwithoutR( Vdbm )
%DBMTOMV 功率(dBmW)转换为电压(V)
%   不考虑电阻
Volt=sqrt(10^((Vdbm-30)/10));
end

