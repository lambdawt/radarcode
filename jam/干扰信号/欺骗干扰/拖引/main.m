load('data_LFMParameter.mat')
load('data_target0Parameter.mat')
load('data_DeceptionJammingParameter.mat')
temp1=2;%1Ϊ�ٶȡ�2Ϊ���롢3Ϊ����
c=3e8;
lamta=c/fz_LFM;
Rj=2e3;
Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L_LFM/10)));
Aj=sqrt(Prj);
tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%һ�������ظ����ڲ�������
N=length(tm);%һ�������ظ����ڲ�����������
f_doppler=2*v/lamta;
f_doppler1=f_doppler;
k=B1_LFM/tau_LFM;%���Ե�Ƶ�źŵ���ϵ��
ts=1/fs_LFM;%�������
%���ɷ����ź�
Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L/10))); %Ŀ��ز��źŹ���
A=sqrt(Prs);
Te=1143;%�¶�
F=4;%����ϵ��
An=10*log10((1.382e-23)*Te*B_LFM*10^(F/10));%����ǿ��
[y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k); 
[~,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y,tr_LFM,ts,N);
%���ɻز��ź�
[s_echo_2,echo]=LFMtuoyinhuiboxinhao(y,R,tr_LFM,ts,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
%figure,plot(0:ts:(N-1)*ts,real(s_echo_2(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('�ز��ź�');
%figure,plot((0:fs/length(abs(fft(abs(fftshift(echo(1,:)))))):fs-fs/length(abs(fft(abs(fftshift(echo(1,:))))))),abs((fft(echo(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�ز��źŵ�Ƶ��');

%���ɸ����ź�
[s_ft,echo3]=LFMtuoyinganrao(R,frame_LFM,tf,Aj,tm,tau_LFM,f0_LFM,B1_LFM,N,f_doppler1,f_doppler,vf,lamta,k,ts,c,fs_LFM,temp1);

%��������
[s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
%Ŀ��ز��źš���Ŀ���źš�����������һ��������ջ�
s_echo_1=s_echo_2+s_noise+s_ft;
figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('�״�����ź�');
figure,plot((0:fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_LFM-fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�״�����źŵ�Ƶ��');
