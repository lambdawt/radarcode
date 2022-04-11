function [ PassivePara ] = paraset_passivejaming( tf,sf,vl,fv,ts,bt,al,sref,smax )
%========================================================================%
%功能：箔条干扰参数设置                                              %
%输入：
%   tf              第K发箔条弹的发射时刻   
%   sf              射程
%   vl              各种干扰弹的射速
%   fv              各种干扰弹的下落速度
%   ts              各种干扰弹的展开时间(从开始展开到完全展开所需要的时间)
%   bt              垂直射角，以海面为准 
%   al              水平射角，以航首方向为准
%   sref            反射面积
%   smax            各种干扰弹的反射面积  
%输出：
%   PassivePara     箔条干扰参数结构体
%========================================================================%
PassivePara.FireTime=tf;
PassivePara.FireScore=sf;
PassivePara.LauchVel=vl;
PassivePara.FallVel=fv;
PassivePara.SpreadTime=ts;
PassivePara.Bt=bt;
PassivePara.Al=al;
PassivePara.Sref=sref;
PassivePara.Smax=smax;

end

