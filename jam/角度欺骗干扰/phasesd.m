function theta_measure = phasesd( S1 , S2 )
c=3e8;                                               
d=2;         %���߼��
fc=1e8;      %�״﹤��Ƶ��
lambda=c/fc; %��������
ysum=S1+S2;
ydif=S2-S1;
k=imag(ydif/ysum);
u=atan(k)*lambda/(pi*d);
theta_measure=asin(u);
theta_measure=theta_measure*180/pi;

end

