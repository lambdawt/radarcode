function [ OutSig ] = PowerWeighted( Vin,mSig )
%POWERWEIGHTED ���ʼ�Ȩ
%   ���ݵ�ǰ�źŹ��ʣ������빦��Pt����ɱ����˵��ź���
%   �������źŵ�ƽ������ΪPin
vSig=mSig;
EnergeSig=sum(real(vSig).^2+imag(vSig).^2);
Vaverage=sqrt(EnergeSig/length(vSig));
vSig=(Vin/Vaverage)*vSig;
OutSig=vSig;
end

