function [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame,match_filter_fft,tau,D,ts)
%%  ����
%   s_echo_mf �ز��ź�
%   M         �������������һ��PRI�ڣ�
%   frame     ����������
%   match_filter_fft    ƥ���˲�ϵ��
%   tau       ����
%   D         ����ѹ����
%   ts        ��������
%%  ���
%   pc_result    ��ѹ���
%   pc_result1   ���������
%   M1           ��ѹ���ź��ܵ���
%%  ����   ƥ���˲�
%%
for i=1:frame
    s_echo_fft_result(:,i)=fft(s_echo_mf(:,i),M);
    pc_result_fft(:,i)=match_filter_fft.*s_echo_fft_result(:,i);
    pc_result1(:,i)=ifft(pc_result_fft(:,i),M);
    pc_result(:,i)=downsample(ifft(pc_result_fft(:,i),M),floor(tau/D/ts),floor(tau/D/ts/2));
end
M1=length(pc_result1(:,i));
pc_result1=pc_result1.';
% M=N;
% figure,plot(0:ts*c/2:(length(abs(pc_result(:,1)))*ts-ts)*c/2,(pc_result(:,1))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('��һ���ز�������ѹ��');
% M1=length(pc_result(:,i));%ԭʼ

% s_pc_result=reshape(pc_result,1,length(pc_result(:,i))*frame);%������תΪ1��M*frame��
% figure,subplot(221),plot(0:ts:(frame*M-1)*ts,abs(s_echo_1)),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('��Ŀ�ꡢ��Ŀ�ꡢ����');axis([-inf inf 0 max(abs(s_pc_result))])
% subplot(222),plot(0:ts:length(abs(s_pc_result))*ts-ts,abs(s_pc_result)),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('����ѹ�������İ���');