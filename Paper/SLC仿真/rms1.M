function  RMS=rms1(x)
%Calculation of RMS(Ratio of Mainlobe to Sidelobe).
t=abs(x);N=length(t);[m,pm]=max(t);s=0;ps=0;SPR=0;%s=t(1);ps=1;SPR=0;
% tt=20*log10(t/max(t));

for k=2:N-1
    if (t(k)>=t(k-1) & t(k)>=t(k+1) &t(k)~=m), SPR=t(k);PN=k; end %SPR即旁瓣高度
    if  (SPR>s), s=SPR;ps=PN; end%选出最大旁瓣值  s  ,和相应序号  ps
end
if  (s<t(N)),  s=t(N);ps=N; end
if (s==0)
    RMS=[];
    else
    RMS=[20*log10(m/s)]; %返回数组格式：[主副比 最大值 最大值序号 旁瓣最大值 旁瓣最大序号]
end
    %revised in Ocotober 2008
