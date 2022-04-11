function [s_mtd,s_mtd1]=MTD_tuoyin(pc_result,M1,num_jilei,num_tongdao)
s_pc=pc_result;
for i=1:M1
    s_temp=s_pc(i,1:num_jilei);
       win=hamming(length(s_temp));
       s_temp=s_temp.*win';
    s_mtd(i,:)=fft(s_temp,num_tongdao);
end

for i=1:M1
    s_temp=s_pc(i,15:46);
       win=hamming(length(s_temp));
       s_temp=s_temp.*win';
    s_mtd1(i,:)=fft(s_temp,num_tongdao);
end