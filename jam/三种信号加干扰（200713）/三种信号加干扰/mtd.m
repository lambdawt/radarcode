function [s_mtd]=mtd(pc_result,M1,num_jilei,num_tongdao)
%%  ����
%   pc_result     ����ѹ�����ź�
%   M1            ��ѹ���ź��ܵ���
%   num_jilei     ���������
%   num_tongdao   ͨ����
%%  ���
%   s_mtd         ��Ŀ������ź�
%%  ����         ��Ŀ����
%%
s_pc=pc_result;
for i=1:M1
    s_temp=s_pc(i,1:num_jilei);
       win=hamming(length(s_temp));
       s_temp=s_temp.*win';
    s_mtd(i,:)=fft(s_temp,num_tongdao);
end