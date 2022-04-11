function [ Vr ] = GetVmt( TargetPosition,InputVelocity )
%========================================================================%
%功能：计算径向速度                                            %
%输入：
%   TargetPosition        目标位置参数结构体   
%   InputVelocity         速度参数结构体
%输出：
%   Vr                    径向速度                  
%========================================================================%
	
	L=sqrt(TargetPosition.m_X^2+TargetPosition.m_Y^2+TargetPosition.m_Z^2);%构造立方体最长对角线
	Vx=InputVelocity.Vx*TargetPosition.m_X/L;
	Vy=InputVelocity.Vy*TargetPosition.m_Y/L;
	Vz=InputVelocity.Vz*TargetPosition.m_Z/L;
	Vr=-Vx-Vy-Vz;

end

