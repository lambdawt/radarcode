%雷达目标有效散射面积
function sigma=rcs(rcsk,sigma0)
%% 输入变量
%rcsk    目标的平均截面积
%sigma0        目标的起伏类型
%% 输出变量
%sigma       目标的散射面积
%% 程序
if rcsk==0 
    sigma=sigma0;  %rcsk =0时采用不起伏模型
elseif rcsk==1
    sigma=-sigma0*log(1-rand(1,1));  %rcsk =1时采用斯威林Ⅰ,Ⅱ起伏模型
elseif rcsk==2
    sigma=-(sigma0/2)*(log(1-rand(1,1))+log(1-rand(1,1))); %rcsk =2时采用斯威林Ⅲ,Ⅳ起伏模型
end