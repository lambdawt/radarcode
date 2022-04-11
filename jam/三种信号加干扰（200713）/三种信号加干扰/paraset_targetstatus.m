function [ TargetStatus ] = paraset_targetstatus( px,py,pz,vx,vy,vz,ax,ay,az,phi )
%========================================================================%
%���ܣ�Ŀ��״̬��������                                              %
%���룺
%   px,py,pz         Ŀ��λ��   
%   vx,vy,vz         Ŀ���ٶ�
%   ax,ay,az         Ŀ����ٶ�
%   phi              ������
%�����
%   TargetStatus     Ŀ��״̬�����ṹ��
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

