function [s_echo_2,echo]=JDhuiboxinhao(R,c,A,N,frame,fs,f_doppler,tm,f0,tau) 
for i=0:frame-1
%     echo(i+1,:)=A.*[zeros(1,fix(2*R/c/ts)),y,zeros(1,ceil((tr-2*R/c)/ts)-length(y))];%����frame�лز��ź�
    %echo(i+1,:)=A*rectpuls(tm-2*(R+v*tm)/c-tau/2,tau).*exp(j*2*pi*(((f0)).*(tm-2*(R+v*tm)/c)));
     echo(i+1,:)=A*rectpuls(tm-2*R/c-tau/2,tau).*exp(1j*2*pi*f0.*(tm-2*R/c));
end
for i=1:frame
    s_echo_1(1,(i-1)*N+1:i*N)=echo(i,:);%��frame�лز��ź��ظ�����ӳ�һ��
end
n=1:frame*N;%��frame�����峤����ͬ�Ĳ�������
s_doppler=exp(1j*n*2*pi*f_doppler/fs);
s_echo_2=s_echo_1.*s_doppler;%�ز��źŵ��Ӷ�����Ƶ��