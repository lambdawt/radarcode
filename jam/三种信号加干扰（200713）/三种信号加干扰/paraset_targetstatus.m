function [ TargetStatus ] = paraset_targetstatus( px,py,pz,vx,vy,vz,ax,ay,az,phi )
%========================================================================%
%功能：目标状态参数设置                                              %
%输入：
%   px,py,pz         目标位置   
%   vx,vy,vz         目标速度
%   ax,ay,az         目标加速度
%   phi              滚动角
%输出：
%   TargetStatus     目标状态参数结构体
%========================================================================%
TargetStatus.TarPosition.m_X=px;
TargetStatus.TarPosition.m_Y=py;
TargetStatus.TarPosition.m_Z=pz;
TargetStatus.TarVelocity.Vx=vx;
TargetStatus.TarVelocity.Vy=vy;
TargetStatus.TarVelocity.Vz=vz;
TargetStatus.TarAcceleretion.m_Ax=ax;
TargetStatus.TarAcceleretion.m_Ay=ay;
TargetStatus.TarAcceleretion.m_Az=az;
TargetStatus.TarAngle.m_Phi=phi;

end

