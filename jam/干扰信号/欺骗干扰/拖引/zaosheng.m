function [s_noise]=zaosheng(frame,N,An,B,fs,ts)
%%  ���� 
%   frame   ����������
%   N       һ�������ظ����ڲ�����������
%   An      ����ǿ��
%   B       ����
%   fs      ����Ƶ��
%%  ���
%   s_noise  �����ź�
%%  ����  ���������ź�
%%
s_noise=wgn(1,frame*N,An);
%hold on,plot(0:ts:(frame*N-1)*ts,abs(s_noise),'g');
%figure,plot(0:fs/length(abs(fft(s_noise))):fs-fs/length(abs(fft(s_noise))),abs(fft(s_noise)));
window=chebwin(51,40);
[b,a]=fir1(50,2*B/fs,window);
s_noise=filter(b,a,s_noise);