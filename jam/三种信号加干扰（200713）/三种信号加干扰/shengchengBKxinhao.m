function [y,y1,D]=shengchengBKxinhao(tau,fs,f0,~,number1,code,Pt,tr,ts)
%% ���� 
%   ����       tau
%   ����Ƶ��   fs
%   ����Ƶ��   f0
%   modelģʽ 1�Ϳ��� 2�Զ���
%   number1   ��������ܳ��ȣ��Ϳ���ģʽ��
%   code      �����Ӧ����2���4�������У��Զ�ģʽ��
%   Pt        �״﷢�书��
%   tr        �����ظ�����
%   ts        ��������
%% ���
%   y         ��λ�����źţ�PRI�ڣ�     
%   y1        ��λ�����źţ�ԭʼ�����źţ�
%   D         ��ѹϵ��
%% ���� ������λ�����ź�
%   ����ѡ��Ϳ�����Զ���
%   ����ѡ������������
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rmin=sqrt((Kj*Pt*sigma)/(4*pi*(Pj)*rj))
% f1=-f0+f1;
%Prs=((Pt*(10^((Gt/10)))*(10^((Gr/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L/10))); %Ŀ��ز��źŹ���
%Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L/10)));%�����źŹ���
%A=sqrt(Prs);%�ز��źŷ���
%Aj=sqrt(Prj);%�����źŷ���
%An=10*log10((1.382e-23)*Te*B*10^(F/10));%����ǿ��
% deltt1=2*(R1-R)/c;%��Ŀ�������Ŀ���ӳ�ʱ��
%fr=1/tr;%�����ظ�Ƶ��
%ts=1/fs; 
%f_doppler=2*v/lamta;%��Ŀ�������Ƶ��
%f_doppler1=2*v1/lamta;%��Ŀ�������Ƶ��
%tm=0:1/fs:tr-1/fs;%һ�������ظ����ڲ�������
%O=tau/ts;%һ�������������
%ft=linspace(-tau/2,tau/2,O);%һ�������������
%N=length(tm);%һ�������ظ����ڲ�����������
%k=B/tau;%B/fs*2*pi/max(ft);  %����ϵ��
%y=exp(j*2*pi*((f0+k*ft/2).*ft));%��������������ͬ�����Ե�Ƶ�ź�

%number1=length��Ԫ����  ,number2= 2��λ [1,-1][0,pi]   4��λ [0,2*pi]/[0,1,2,3]phase[0,pi/2,pi,3pi/2]
%code=[���롭��]
%0�Ϳ���or�Զ��� 
%0.1�Ϳ���  0.1.1ѡ��λ��7λ 0.1.2ѡ��λ��13λ 
%0.2�Զ���  0.2.1ѡ������� 0.2.2 ѡ��������  0.2.*.1��������Ԫ���� 0.2.*.2���������       
%model=2;
%number2=4;
% % if model==1
% %     code=[1,1,1,-1,-1,-1,1];
% % else if model==2
% %         code=[1,1,1,1,1,-1,-1,1,1,-1,1,-1,1];%13λ�Ϳ���
% %     else
% %         code=code;
% %     end
% % end
        

% code=[1,1,1,1,1,-1,-1,1,1,-1,1,-1,1];%13λ�Ϳ���%code=
%code=[1,1,1,-1,-1,-1,1];%7λ�Ϳ���
%code=[0,1,2,3,0,1,2,3,0,1,2,3,0];%����ָ�������ı�
D=number1;
ncode=length(code);%number1
tao=tau/ncode; %number1=13
t_tao=0:1/fs:tao-1/fs;
Ntao=length(t_tao);
pha=0;
t=0:1/fs:ncode*tao-1/fs;%number1=13
y1=zeros(1,length(t));
for i=1:ncode
   % if code(i)==1
        pha=code(i)*pi/2;
    %else pha=0;
    %end
    y1(1,(i-1)*Ntao+1:i*Ntao)=exp(1j*(2*pi*f0*t_tao+pha));
end
y=sqrt(Pt).*[y1,zeros(1,fix((tr-tau)/ts))];


%figure(1),subplot(2,1,1),
%figure,plot(real(y)),xlabel('t(��λ:s)'),title('�����ź�');
% s_fft_result=abs(fft(y(1:Ntao)));
%subplot(2,1,2),
%figure,plot((0:fs/Ntao:fs/2-fs/Ntao),abs(s_fft_result(1:Ntao/2))),xlabel('Ƶ��(��λ:Hz)'),title('�����ź�Ƶ��');