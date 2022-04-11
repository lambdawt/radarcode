function [cfar_k_result]=hengxujing(M1,Pfa,s_mtd,num_cankao,num_tongdao,num_baohu,ts,c,tau,D)
%%  ����
%   M1            ��ѹ���ź��ܵ���
%   Pfa           �龯��
%   s_mtd         ��Ŀ�������
%   num_cankao    �ο���Ԫ��
%   num_tongdao   ͨ����
%   num_baohu     ������Ԫ��
%   ts            ��������
%   c             ����
%   tau           ����
%   D             ����ѹ����
%%  ���
%%  ����
%%
N=M1;
Kc=sqrt(-2*log(Pfa));
h=figure('Name','���ú��龯������');
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
%     subplot(4,num_tongdao/4,v),plot(0:ts*(tau/D/ts)*c/2:(length(abs(cfar_k_result))*ts*(tau/D/ts)-ts)*c/2,abs(cfar_k_result)),xlabel('���룬��λ����'),ylabel('y(��λ����)');axis([0 12000 0 cfarmax]); 
    subplot(4,num_tongdao/4,v),plot(0:ts*c/2:(length(abs(cfar_k_result))*ts-ts)*c/2,abs(cfar_k_result)),xlabel('���룬��λ����'),ylabel('y(��λ����)');axis([0 12000 0 cfarmax]); 
   
%set(gca,'xtick',[0:1000:18000]);
    grid on;
end