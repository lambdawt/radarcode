function [cfar_k_result]=hengxujing(M1,Pfa,s_mtd,num_cankao,num_tongdao,num_baohu,ts,c,tau,D)
%%  输入
%   M1            脉压后信号总点数
%   Pfa           虚警率
%   s_mtd         动目标检测后结果
%   num_cankao    参考单元数
%   num_tongdao   通道数
%   num_baohu     保护单元数
%   ts            采样周期
%   c             光速
%   tau           脉宽
%   D             脉冲压缩比
%%  输出
%%  功能
%%
N=M1;
Kc=sqrt(-2*log(Pfa));
h=figure('Name','采用恒虚警处理结果');
cfarmax=max(max(abs(s_mtd)));
for i=1:num_tongdao
        v=i;
        s_mtd_1=abs(s_mtd(:,i));
        s_pc1=s_mtd_1.';
    %     N=4096;
    cfar_k_result(1,1)=mean(s_pc1(1,2:num_cankao+1));
    if s_pc1(1,1)>=sqrt(2/pi)*mean(s_pc1(1,2:num_cankao+1))*Kc
        cfar_k_result(1,1)=s_pc1(1,1);
    else
        cfar_k_result(1,1) = 0;
    end
    
    for i=2:num_cankao
        noise_mean=(mean(s_pc1(1,1:i-1))+mean(s_pc1(1,i+1:i+num_cankao)))/2;
        cfar_k_result(1,i)=noise_mean;
        if s_pc1(1,i)>=sqrt(2/pi)*noise_mean*Kc
            cfar_k_result(1,i)=s_pc1(1,i);
        else
            cfar_k_result(1,i)=0;
        end

    end
    
    for i=num_cankao+1:N-num_cankao-1
        noise_mean=max(mean(s_pc1(1,i-num_cankao:i-num_baohu)),mean(s_pc1(1,i+num_baohu:i+num_cankao)));
        cfar_k_result(1,i)=noise_mean;
        if s_pc1(1,i)>=sqrt(2/pi)*noise_mean*Kc
            cfar_k_result(1,i)=s_pc1(1,i);
        else
            cfar_k_result(1,i)=0;
        end
    end
    for i=N-num_cankao:N-1
        noise_mean=mean(s_pc1(1,i-num_cankao:i-1))+mean(s_pc1(1,i+1:N));
        cfar_k_result(1,i)=s_pc1(1,i)/noise_mean;
        if s_pc1(1,i)>=sqrt(2/pi)*noise_mean*Kc
            cfar_k_result(1,i)=s_pc1(1,i);
        else
            cfar_k_result(1,i)=0;
        end
    end
    cfar_k_result(1,N)=(mean(s_pc1(1,N-num_cankao:N-1)));
        if s_pc1(1,N)>=sqrt(2/pi)*mean(s_pc1(1,N-num_cankao:N-1))
            cfar_k_result(1,N)=s_pc1(1,N);
        else
            cfar_k_result(1,N)=0;
        end
%     subplot(4,num_tongdao/4,v),plot(0:ts*(tau/D/ts)*c/2:(length(abs(cfar_k_result))*ts*(tau/D/ts)-ts)*c/2,abs(cfar_k_result)),xlabel('距离，单位：米'),ylabel('y(单位：伏)');axis([0 12000 0 cfarmax]); 
    subplot(4,num_tongdao/4,v),plot(0:ts*c/2:(length(abs(cfar_k_result))*ts-ts)*c/2,abs(cfar_k_result)),xlabel('距离，单位：米'),ylabel('y(单位：伏)');axis([0 12000 0 cfarmax]); 
   
%set(gca,'xtick',[0:1000:18000]);
    grid on;
end