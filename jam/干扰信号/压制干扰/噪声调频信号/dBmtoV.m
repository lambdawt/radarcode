function [ Volt ] = dBmtoV( Vdbm )
%DBMTOMV 功率(dBmW)转换为电压(V)
%   先将功率由dBmW换算成W
%   再根据公式P=U^2/R 求得电压U
%   其中R为50欧
Volt=sqrt(50*10^((Vdbm-30)/10));
end

