function [s_echo_1,f0]=hunpin(s_echo_1,N,frame,f1,fs,f0)
s_echo_1=reshape(s_echo_1,N,frame);
n=0:N-1;
s_echo1=zeros(N,frame);
for i=1:frame
    s_echo_1(:,i)=s_echo_1(:,i).*exp(-1j*n*2*pi*(f1)/fs).';%(2*cos(n*f1/fs*2*pi)).'     %把回波信号的每一个都换成中频，e上的相乘相当于加减
end
s_echo_1=reshape(s_echo_1,1,N*frame);
f0=f0-f1;
