function theta_measure = phasesd( S1 , S2 )
c=3e8;                                               
d=2;         %天线间距
fc=1e8;      %雷达工作频率
lambda=c/fc; %工作波长
ysum=S1+S2;
ydif=S2-S1;
k=imag(ydif/ysum);
u=atan(k)*lambda/(pi*d);
theta_measure=asin(u);
theta_measure=theta_measure*180/pi;

end

