function [s_ft,echo3]=LFMDeceptionJamming(D,~,R1,tr,~,c,Aj,N,frame,fs,f_doppler1,tm,f0,B1,tau,k,congmubiao,temp1)
s_ft=zeros(1,frame*N);%����1�У�frame*N�и�0���󣬱��ڲ��������Ŀ����
delt1=2*R1/c;
if temp1 == 3 %����
    for i1=1:congmubiao
    for i=0:frame-1
    %     echo(i+1,:)=A.*[zeros(1,fix(2*R/c/ts)),y,zeros(1,ceil((tr-2*R/c)/ts)-length(y))];%����frame�лز��ź�
        echo3(i+1,:)=Aj*rectpuls(tm-delt1-tau/2,tau).*exp(1j*2*pi*(((f0-B1/2)+k*(tm-delt1)/2).*(tm-delt1)));

    end
    for i=1:frame
        s_echo_1(1,(i-1)*N+1:i*N)=echo3(i,:);%��frame�лز��ź��ظ�����ӳ�һ��
    end
    n=1:frame*N;%��frame�����峤����ͬ�Ĳ�������
    s_doppler1=exp(1j*n*2*pi*f_doppler1/fs);
    s_echo_3=s_echo_1.*s_doppler1;%�ز��źŵ��Ӷ�����Ƶ��
    s_ft=s_ft+s_echo_3;
    f_doppler1=f_doppler1-4/(frame*tr);
    delt1=delt1+10*tau/D;
    end
elseif temp1 == 2 %����
    for i1=1:congmubiao
    for i=0:frame-1
    %     echo(i+1,:)=A.*[zeros(1,fix(2*R/c/ts)),y,zeros(1,ceil((tr-2*R/c)/ts)-length(y))];%����frame�лز��ź�
        echo3(i+1,:)=Aj*rectpuls(tm-delt1-tau/2,tau).*exp(1j*2*pi*(((f0-B1/2)+k*(tm-delt1)/2).*(tm-delt1)));

    end
    for i=1:frame
        s_echo_1(1,(i-1)*N+1:i*N)=echo3(i,:);%��frame�лز��ź��ظ�����ӳ�һ��
    end
    n=1:frame*N;%��frame�����峤����ͬ�Ĳ�������
    s_doppler1=exp(1j*n*2*pi*f_doppler1/fs);
    s_echo_3=s_echo_1.*s_doppler1;%�ز��źŵ��Ӷ�����Ƶ��
    s_ft=s_ft+s_echo_3;
    %f_doppler1=f_doppler1-4/(frame*tr);
    delt1=delt1+10*tau/D;
    end
elseif temp1 == 1 %�ٶ�
    for i1=1:congmubiao
    for i=0:frame-1
    %     echo(i+1,:)=A.*[zeros(1,fix(2*R/c/ts)),y,zeros(1,ceil((tr-2*R/c)/ts)-length(y))];%����frame�лز��ź�
        echo3(i+1,:)=Aj*rectpuls(tm-delt1-tau/2,tau).*exp(1j*2*pi*(((f0-B1/2)+k*(tm-delt1)/2).*(tm-delt1)));

    end
    for i=1:frame
        s_echo_1(1,(i-1)*N+1:i*N)=echo3(i,:);%��frame�лز��ź��ظ�����ӳ�һ��
    end
    n=1:frame*N;%��frame�����峤����ͬ�Ĳ�������
    s_doppler1=exp(1j*n*2*pi*f_doppler1/fs);
    s_echo_3=s_echo_1.*s_doppler1;%�ز��źŵ��Ӷ�����Ƶ��
    s_ft=s_ft+s_echo_3;
    f_doppler1=f_doppler1-4/(frame*tr);
    %delt1=delt1+tau/D;
    end
end

