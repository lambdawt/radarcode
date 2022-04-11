function [ sig_jam,t_jam ] = jam_passive( echo,fs,fc,CurrentT,TargetStatus,WindV,PassivePara )
%========================================================================%
%���ܣ���������                                              %
%���룺
%   echo            �ز��ź�   
%   fs          	����Ƶ��
%   fc              �ز��ź���Ƶ
%   CurrentT       	��ǰʱ��
%   TargetStatus    Ŀ��״̬���ṹ�壩
%   WindV          	����
%   PassivePara     ��������

%�����
%   sig_jam         ���������ź�
%   t_jam           ����źŵ�ʱ����
%========================================================================%
global c;
	T1=PassivePara.FireTime+PassivePara.FireScore/PassivePara.LauchVel;%����T1Ϊ������ʼչ��ʱ��

	T2=T1+PassivePara.SpreadTime;                                      %����T2Ϊ������ȫչ��ʱ��
	
	T3=T2+PassivePara.FireScore*sin(PassivePara.Bt)/PassivePara.FallVel;
 
	PassTime=CurrentT-PassivePara.FireTime;                            %����������ʱ��ͻ������͵�ʱ����

	Tn=CurrentT-T1;                                                    %�Ӳ�����ʼչ�������ڵ�ʱ��

	Tm=CurrentT-T2;                                                     %��������˶�ʱ��
    
	%�ȼ��㲭������ʱĿ���λ�ã����ڻ����ṩ���Ƕ�Ӧʱ���λ�ã�ͬʱ�ṩ��Ŀ����ٶȣ����ٶȣ��������ȵȣ��ɸ���ʱ����������λ��
	CurTarX=TargetStatus.TarPosition.m_X-TargetStatus.TarVelocity.Vx*PassTime-0.5*TargetStatus.TarAcceleretion.m_Ax*PassTime*PassTime;

	CurTarY=TargetStatus.TarPosition.m_Y-TargetStatus.TarVelocity.Vy*PassTime-0.5*TargetStatus.TarAcceleretion.m_Ay*PassTime*PassTime;

	CurTarZ=TargetStatus.TarPosition.m_Z-TargetStatus.TarVelocity.Vz*PassTime-0.5*TargetStatus.TarAcceleretion.m_Az*PassTime*PassTime;

	StartX=CurTarX+PassivePara.FireScore*cos(PassivePara.Bt)*sin(PassivePara.Al+TargetStatus.TarAngle.m_Phi);

	EndX=StartX+PassivePara.FireScore*cos(PassivePara.Bt)*sin(PassivePara.Al+TargetStatus.TarAngle.m_Phi);

	StartY=CurTarY+PassivePara.FireScore*cos(PassivePara.Bt)*cos(PassivePara.Al+TargetStatus.TarAngle.m_Phi);

	EndY=StartY+PassivePara.FireScore*cos(PassivePara.Bt)*cos(PassivePara.Al+TargetStatus.TarAngle.m_Phi);

	Z=CurTarZ+PassivePara.FireScore*sin(PassivePara.Bt);
	
    if (CurrentT>T1&&CurrentT<T2)

		%������������չ����״̬�����㾶���ٶ�ʱֻ�������� 
		Chaff_X=(StartX+EndX)/2;        %����X����
		Chaff_Y=(StartY+EndY)/2;        %����Y����
		Chaff_Z=Z;                      %����Z����
		Xmt=sqrt(Chaff_X^2+Chaff_Y^2+Chaff_Z^2);  %�������״�/�����ľ���.�״�λ��Ϊ����ԭ�㣬���㲭�������״��λ�ã��Ӷ���������ź�����״��źŵ��ӳ�ʱ��
		PassivePara.Sref=PassivePara.Smax*(1-exp(-2.3*Tn/PassivePara.SpreadTime));
		%���㾶���ٶ�,�ֽ����ٷֽ⵽x,y,z���������ڷֱ�ͶӰ����������
		InputVelocity.Vx=PassivePara.LauchVel*cos(PassivePara.Bt)*sin(PassivePara.Al+TargetStatus.TarAngle.m_Phi);
		InputVelocity.Vy=PassivePara.LauchVel*cos(PassivePara.Bt)*cos(PassivePara.Al+TargetStatus.TarAngle.m_Phi);
		InputVelocity.Vz=PassivePara.LauchVel*sin(PassivePara.Bt);
		Vmt=GetVmt(TargetStatus.TarPosition,InputVelocity);%���㾶���ٶȣ����������Ƶ��
        
    elseif (CurrentT>T2&&CurrentT<T3)
	
		Chaff_X=(StartX+EndX)/2+WindV.Vx*Tm;
		Chaff_Y=(StartY+EndY)/2+WindV.Vy*Tm;
		Chaff_Z=Z+PassivePara.FallVel*Tm;
		Xmt=sqrt(Chaff_X^2+Chaff_Y^2+Chaff_Z^2);   %�������״�ľ���
		PassivePara.Sref=PassivePara.Smax*(1-exp(-2.3*Tn/PassivePara.SpreadTime));
		
		InputVelocity.Vx=WindV.Vx;
		InputVelocity.Vy=WindV.Vy;
		InputVelocity.Vz=PassivePara.FallVel;
		Vmt=GetVmt(TargetStatus.TarPosition,InputVelocity);
	
    else
		PassivePara.Sref=0;  
		Xmt=0;   %��ʱ�����Ѿ����������ã����ò������״�ľ���Ϊ0���Ӷ�����������ж�
    end
	sj=zeros(1,length(echo));
    t=(0:length(echo)-1)/fs;
    if (Xmt~=0)   
        c=3e8;
		DelayTime=2*Xmt/c;       %�������״����Ծ�����ɵ��ź�ʱ���ӳ�
		DelayNum=DelayTime*360e-3/1e-6;   %�ӳٵ���
		WaveLength=c/fc;         %�״��źŲ���
		fd=2*Vmt/WaveLength;     %������Ƶ��
        t=t+DelayNum;
		sj=20*PassivePara.Sref*echo.*exp(1j*2*pi*fd*t);
    end
sig_jam=sj;
t_jam=t;
end

