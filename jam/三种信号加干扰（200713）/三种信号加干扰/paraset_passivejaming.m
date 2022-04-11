function [ PassivePara ] = paraset_passivejaming( tf,sf,vl,fv,ts,bt,al,sref,smax )
%========================================================================%
%���ܣ��������Ų�������                                              %
%���룺
%   tf              ��K���������ķ���ʱ��   
%   sf              ���
%   vl              ���ָ��ŵ�������
%   fv              ���ָ��ŵ��������ٶ�
%   ts              ���ָ��ŵ���չ��ʱ��(�ӿ�ʼչ������ȫչ������Ҫ��ʱ��)
%   bt              ��ֱ��ǣ��Ժ���Ϊ׼ 
%   al              ˮƽ��ǣ��Ժ��׷���Ϊ׼
%   sref            �������
%   smax            ���ָ��ŵ��ķ������  
%�����
%   PassivePara     �������Ų����ṹ��
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

