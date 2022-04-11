function [s_echo_2,echo]=BKhuiboxinhao(y1,R,tr,ts,A,N,frame,fs,f_doppler,tau)
global c;
c=3e8;
for i=0:frame-1
    echo(i+1,:)=A.*[zeros(1,floor(2*R/c/ts)),y1,zeros(1,ceil((tr-2*R/c-tau)/ts))];%生成frame行回波信号
%    tm1=0:1/fs:2*R/c-1/fs;
%    tm2=0:1/fs:tr-2*R/c-tau-1/fs;
%     y1=rectpuls(tm1-R/c,(2*R/c));
%     y2=rectpuls(tm2-R/c-tau-(tr-2*R/c-tau)/2,tr-2*R/c-tau);
%     echo(i+1,:)=A.*[y1,y,y2];
end
for i=1:frame
    s_echo_1(1,(i-1)*N+1:i*N)=echo(i,:);%将frame行回波信号重复脉冲接成一行
   end
n=1:frame*N;%与frame个脉冲长度相同的采样序列
s_doppler=exp(1j*n*2*pi*f_doppler/fs);
s_echo_2=s_echo_1.*s_doppler;%回波信号叠加多普勒频率