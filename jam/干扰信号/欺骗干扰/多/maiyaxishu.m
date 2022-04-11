function [M,match_filter_fft]=maiyaxishu(f0,fs,y,tr,ts,N)
%%  ����
%   f0    ����Ƶ��
%   fs    ����Ƶ��
%   y     �����ź�
%   tr    �����ظ�����
%   ts    ��������
%   N     һ�������ظ����ڲ�����������
%%  ���
%   M     �������������һ��PRI�ڣ�
%   match_filter_fft     ƥ���˲���ϵ��
%%  ����   ��������ѹ��ϵ����ƥ���˲���Ƶ����Ӧ��
%%


n=0:N-1;%��Ƶ�ź�ʱ�����еĳ���
local_oscillator_i=cos(n*(-f0)/fs*2*pi);
local_oscillator_q=sin(n*(-f0)/fs*2*pi);
fbb_i=local_oscillator_i.*y;
fbb_q=local_oscillator_q.*y;
% window=chebwin(51,40);%�б�ѩ��
% [b,a]=fir1(50,2*B/fs,window);%��ͨ�˲�������ֹƵ��Ϊ2*B/fs
% fbb_i=[fbb_i,zeros(1,25)];
% fbb_q=[fbb_q,zeros(1,25)];
% fbb_i=filter(b,a,fbb_i);%һά�����˲�����aΪ��ĸ��bΪ����
% fbb_q=filter(b,a,fbb_q);
% fbb_i=fbb_i(26:end);
% fbb_q=fbb_q(26:end);
fbb=fbb_i+1j*fbb_q;
% figure,plot(ft,fbb_i),xlabel('t(��λ����)'),ylabel('y(��λ����)'),title('������I·�ź�');
% figure,plot(ft,fbb_q),xlabel('t(��λ����)'),ylabel('y(��λ����)'),title('������Q·�ź�');
% figure,plot(0:fs/length(abs(fft(fbb))):fs-fs/length(abs(fft(fbb))),abs((fft(fbb)))),xlabel('f(��λ��Hz)'),ylabel('y'),title('������źŵ�Ƶ��');
%%%%%%%%%%%%%%%%%%%����������ѹϵ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=fix(tr/ts);%2^nextpow2(tr/ts);%MΪ�ӽ�PRI���ȵ�2��K�η�
match_filter=fliplr(conj(fbb));%conjΪ���fliprΪ�����תѹ�����������sqrt(D)����FFTҪ�����ʵ���ȣ�Ӧ������2����n
% win=hamming(O)';
% match_filter=match_filter.*win;
match_filter_fft_1=fft(match_filter,M);
match_filter_fft=match_filter_fft_1.';
% figure,plot(real(match_filter_fft)),title('����ѹ��ϵ��(ʵ����');
% figure,plot(imag(match_filter_fft)),title('����ѹ��ϵ��(�鲿��');