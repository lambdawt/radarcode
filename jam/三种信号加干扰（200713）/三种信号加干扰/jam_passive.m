function [ sig_jam,t_jam ] = jam_passive( echo,fs,fc,CurrentT,TargetStatus,WindV,PassivePara )
%========================================================================%
%功能：箔条干扰                                              %
%输入：
%   echo            回波信号   
%   fs          	采样频率
%   fc              回波信号载频
%   CurrentT       	当前时间
%   TargetStatus    目标状态（结构体）
%   WindV          	风速
%   PassivePara     箔条参数

%输出：
%   sig_jam         箔条干扰信号
%   t_jam           输出信号的时间轴
%========================================================================%
global c;
	T1=PassivePara.FireTime+PassivePara.FireScore/PassivePara.LauchVel;%定义T1为箔条开始展开时间

	T2=T1+PassivePara.SpreadTime;                                      %定义T2为箔条完全展开时间
	
	T3=T2+PassivePara.FireScore*sin(PassivePara.Bt)/PassivePara.FallVel;
 
	PassTime=CurrentT-PassivePara.FireTime;                            %箔条弹发射时间和环境发送的时间间隔

	Tn=CurrentT-T1;                                                    %从箔条开始展开到现在的时间

	Tm=CurrentT-T2;                                                     %箔条随风运动时长
    
	%先计算箔条发射时目标的位置，由于环境提供的是对应时间的位置，同时提供了目标的速度，加速度，方向俯仰等等，可根据时间计算所需的位置
	CurTarX=TargetStatus.TarPosition.m_X-TargetStatus.TarVelocity.Vx*PassTime-0.5*TargetStatus.TarAcceleretion.m_Ax*PassTime*PassTime;

	CurTarY=TargetStatus.TarPosition.m_Y-TargetStatus.TarVelocity.Vy*PassTime-0.5*TargetStatus.TarAcceleretion.m_Ay*PassTime*PassTime;

	CurTarZ=TargetStatus.TarPosition.m_Z-TargetStatus.TarVelocity.Vz*PassTime-0.5*TargetStatus.TarAcceleretion.m_Az*PassTime*PassTime;

	StartX=CurTarX+PassivePara.FireScore*cos(PassivePara.Bt)*sin(PassivePara.Al+TargetStatus.TarAngle.m_Phi);

	EndX=StartX+PassivePara.FireScore*cos(PassivePara.Bt)*sin(PassivePara.Al+TargetStatus.TarAngle.m_Phi);

	StartY=CurTarY+PassivePara.FireScore*cos(PassivePara.Bt)*cos(PassivePara.Al+TargetStatus.TarAngle.m_Phi);

	EndY=StartY+PassivePara.FireScore*cos(PassivePara.Bt)*cos(PassivePara.Al+TargetStatus.TarAngle.m_Phi);

	Z=CurTarZ+PassivePara.FireScore*sin(PassivePara.Bt);
	
    if (CurrentT>T1&&CurrentT<T2)

		%箔条处于正在展开的状态，计算径向速度时只考虑射速 
		Chaff_X=(StartX+EndX)/2;        %箔条X坐标
		Chaff_Y=(StartY+EndY)/2;        %箔条Y坐标
		Chaff_Z=Z;                      %箔条Z坐标
		Xmt=sqrt(Chaff_X^2+Chaff_Y^2+Chaff_Z^2);  %箔条与雷达/导弹的距离.雷达位置为坐标原点，计算箔条距离雷达的位置，从而计算干扰信号相对雷达信号的延迟时间
		PassivePara.Sref=PassivePara.Smax*(1-exp(-2.3*Tn/PassivePara.SpreadTime));
		%计算径向速度,现将射速分解到x,y,z三个方向，在分别投影到径向后相加
		InputVelocity.Vx=PassivePara.LauchVel*cos(PassivePara.Bt)*sin(PassivePara.Al+TargetStatus.TarAngle.m_Phi);
		InputVelocity.Vy=PassivePara.LauchVel*cos(PassivePara.Bt)*cos(PassivePara.Al+TargetStatus.TarAngle.m_Phi);
		InputVelocity.Vz=PassivePara.LauchVel*sin(PassivePara.Bt);
		Vmt=GetVmt(TargetStatus.TarPosition,InputVelocity);%计算径向速度，计算多普勒频移
        
    elseif (CurrentT>T2&&CurrentT<T3)
	
		Chaff_X=(StartX+EndX)/2+WindV.Vx*Tm;
		Chaff_Y=(StartY+EndY)/2+WindV.Vy*Tm;
		Chaff_Z=Z+PassivePara.FallVel*Tm;
		Xmt=sqrt(Chaff_X^2+Chaff_Y^2+Chaff_Z^2);   %箔条到雷达的距离
		PassivePara.Sref=PassivePara.Smax*(1-exp(-2.3*Tn/PassivePara.SpreadTime));
		
		InputVelocity.Vx=WindV.Vx;
		InputVelocity.Vy=WindV.Vy;
		InputVelocity.Vz=PassivePara.FallVel;
		Vmt=GetVmt(TargetStatus.TarPosition,InputVelocity);
	
    else
		PassivePara.Sref=0;  
		Xmt=0;   %此时箔条已经不再有作用，设置箔条到雷达的距离为0，从而跳出下面的判断
    end
	sj=zeros(1,length(echo));
    t=(0:length(echo)-1)/fs;
    if (Xmt~=0)   
        c=3e8;
		DelayTime=2*Xmt/c;       %箔条与雷达的相对距离造成的信号时间延迟
		DelayNum=DelayTime*360e-3/1e-6;   %延迟点数
		WaveLength=c/fc;         %雷达信号波长
		fd=2*Vmt/WaveLength;     %多普勒频移
        t=t+DelayNum;
		sj=20*PassivePara.Sref*echo.*exp(1j*2*pi*fd*t);
    end
sig_jam=sj;
t_jam=t;
end

