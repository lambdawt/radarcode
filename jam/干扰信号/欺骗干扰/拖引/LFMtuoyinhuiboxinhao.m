function [s_echo_2,echo]=LFMtuoyinhuiboxinhao(y,R,tr,ts,c,A,N,frame,fs,f_doppler,tm,f0,B1,tau,k)
for i=0:frame-1
%     echo(i+1,:)=A.*[zeros(1,fix(2*R/c/ts)),y,zeros(1,ceil((tr-2*R/c)/ts)-length(y))];%����frame�лز��ź�
    echo(i+1,:)=A*rectpuls(tm-2*R/c-tau/2,tau).*exp(1j*2*pi*(((f0-B1/2)+k*(tm-2*R/c)/2).*(tm-2*R/c)));

end
for i=1:frame
    s_echo_1(1,(i-1)*N+1:i*N)=echo(i,:);%��frame�лز��ź��ظ�����ӳ�һ��
end
n=1:frame*N;%��frame�����峤����ͬ�Ĳ�������
s_doppler=exp(1j*n*2*pi*f_doppler/fs);
s_echo_2=s_echo_1.*s_doppler;%�ز��źŵ��Ӷ�����Ƶ��