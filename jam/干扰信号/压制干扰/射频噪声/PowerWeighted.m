function [ OutSig ] = PowerWeighted( Vin,mSig )
%POWERWEIGHTED 功率加权
%   根据当前信号功率，将输入功率Pt换算成比例乘到信号上
%   结果输出信号的平均功率为Pin
vSig=mSig;
EnergeSig=sum(real(vSig).^2+imag(vSig).^2);
Vaverage=sqrt(EnergeSig/length(vSig));
vSig=(Vin/Vaverage)*vSig;
OutSig=vSig;
end

