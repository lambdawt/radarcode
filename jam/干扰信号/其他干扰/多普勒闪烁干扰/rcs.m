%�״�Ŀ����Чɢ�����
function sigma=rcs(rcsk,sigma0)
%% �������
%rcsk    Ŀ���ƽ�������
%sigma0        Ŀ����������
%% �������
%sigma       Ŀ���ɢ�����
%% ����
if rcsk==0 
    sigma=sigma0;  %rcsk =0ʱ���ò����ģ��
elseif rcsk==1
    sigma=-sigma0*log(1-rand(1,1));  %rcsk =1ʱ����˹���֢�,�����ģ��
elseif rcsk==2
    sigma=-(sigma0/2)*(log(1-rand(1,1))+log(1-rand(1,1))); %rcsk =2ʱ����˹���֢�,�����ģ��
end