function [ Vr ] = GetVmt( TargetPosition,InputVelocity )
%========================================================================%
%���ܣ����㾶���ٶ�                                            %
%���룺
%   TargetPosition        Ŀ��λ�ò����ṹ��   
%   InputVelocity         �ٶȲ����ṹ��
%�����
%   Vr                    �����ٶ�                  
%========================================================================%
	
	L=sqrt(TargetPosition.m_X^2+TargetPosition.m_Y^2+TargetPosition.m_Z^2);%������������Խ���
	Vx=InputVelocity.Vx*TargetPosition.m_X/L;
	Vy=InputVelocity.Vy*TargetPosition.m_Y/L;
	Vz=InputVelocity.Vz*TargetPosition.m_Z/L;
	Vr=-Vx-Vy-Vz;

end

