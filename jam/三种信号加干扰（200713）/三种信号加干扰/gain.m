function G = gain(radar,AzAngle,ElAngle)
%========================================================================%
%���ܣ����㷽��ͼ����                                                     %
%���룺�״����radar��Ŀ�����״������ӳ���λ��AzAngle��������ElAngle       %
%���������ͼ����Gainsum,GainAz,GainEl                                   %
%========================================================================%
% Gain = 40;                                          %����������棬��ʱ
% antennaMode = radar.antennaMode;                    %����ͼ����
% switch antennaMode
%     case 'SINC'
%         antenaGain = abs(sinc(-pi:0.001:pi).^2);
%         s = trapz(antenaGain)*2*pi/length(antenaGain); %integral over antena gain
%         antenaGain = antenaGain/s;
%         GainAz = antenaGain(fix(pi/0.001+AzAngle/180*pi/0.001));
%         GainEl = antenaGain(fix(pi/0.001+ElAngle/180*pi/0.001));
%         G = Gain * GainAz * GainEl;      
% end
%% ������ʼ��

Gain = 10^(radar.Gt/10);                                     %��������
antennaAzAngle = radar.antennaAzAngle;                      %���߷�λָ��
antennaELAngle = radar.antennaELAngle;                      %���߸���ָ��
AdjacentAngle = radar.AdjacentAngle;                         %������б�ӽ�
beamwidth = radar.antennabeamwidth;
antennaMode = radar.strfunc{radar.nfunc};
% %% ��ͼ�ò���
% % AzAngle = -(AdjacentAngle+beamwidth):0.1:(AdjacentAngle+beamwidth);
% % ElAngle = -(AdjacentAngle+beamwidth):0.1:(AdjacentAngle+beamwidth);
% AzAngle = -3:0.1:3;
% ElAngle = -3:0.1:3;
%% ������
        antennaLU = subAntennaPattern(antennaMode,beamwidth,AzAngle,ElAngle,...
            antennaAzAngle-AdjacentAngle,antennaELAngle+AdjacentAngle); %����������
        antennaLD = subAntennaPattern(antennaMode,beamwidth,AzAngle,ElAngle,...
            antennaAzAngle-AdjacentAngle,antennaELAngle-AdjacentAngle); %����������
        antennaRU = subAntennaPattern(antennaMode,beamwidth,AzAngle,ElAngle,...
            antennaAzAngle+AdjacentAngle,antennaELAngle+AdjacentAngle); %����������
        antennaRD = subAntennaPattern(antennaMode,beamwidth,AzAngle,ElAngle,...
            antennaAzAngle+AdjacentAngle,antennaELAngle-AdjacentAngle); %����������
        G.Gainsum =(antennaRU + antennaRD + antennaLU + antennaLD);        %�Ͳ���
        G.GainAz = Gain*(antennaRU + antennaRD - antennaLU - antennaLD);         %��λ���
        G.GainEl = Gain*(antennaRU - antennaRD + antennaLU - antennaLD);         %�������
        G.NDSAz =G.GainAz./G.Gainsum;             %��λ��һ����б��normalization difference slop 
        G.NDSEl =G.GainEl./G.Gainsum;             %������һ����б��normalization difference slop 

%% ��ͼ
% figure(31);
% mesh(AzAngle*pi/180,ElAngle*pi/180,antennaLU);
% hold on;
% mesh(AzAngle*pi/180,ElAngle*pi/180,antennaLD);
% mesh(AzAngle*pi/180,ElAngle*pi/180,antennaRU);
% mesh(AzAngle*pi/180,ElAngle*pi/180,antennaRD);
% xlabel('\theta��^����' );
% ylabel('\phi��^����' );
% zlabel('���ȣ�^����' );
% title('�ĸ�������Ӧ');
% grid on;

% figure(32);
% subplot(411);
% plot(AzAngle*pi/180,antennaLU(:,round(length(AzAngle)/2)),'r');
% hold on;
% plot(AzAngle*pi/180,antennaLU(round(length(ElAngle)/2),:),'b');
% xlabel('\theta��^����' );
% ylabel('���ȣ�^����' );
% title('����������Ӧ');
% grid on;
% 
% subplot(412);
% plot(AzAngle*pi/180,G.Gainsum(:,round(length(AzAngle)/2)));
% grid on;
% xlabel('\theta��^����' );
% title('�Ͳ�����Ӧ');
% 
% subplot(413);
% plot(AzAngle*pi/180,G.GainAz(:,round(length(AzAngle)/2)),'r');
% % hold on;
% % plot(ElAngle*pi/180,G.GainEl(round(length(ElAngle)/2),:),'b');
% grid on;
% xlabel('\theta��^����' );
% title('�����Ӧ');
% 
% subplot(414); 
% plot(AzAngle*pi/180,G.NDSAz(:,round(length(AzAngle)/2)),'r');
% grid on;
% xlabel('\theta��^����' );
% title('���/�Ͳ���');

%�����ӷ���ͼ
function subGain = subAntennaPattern(antennaMode,beamwidth,AzAngle,ElAngle,antennaAz,antennaEl)
%========================================================================%
%���ܣ������ӷ���ͼ����                                                   %
%���룺����ͼ����antennaMode���������beamwidth���Ƕ�angle                 %
%������ӷ���ͼ����subGain                                                        %
%========================================================================%
%�ǶȻ��ɻ���
AzAngle = AzAngle/180*pi;
ElAngle = ElAngle/180*pi;
antennaAz = antennaAz/180*pi;
antennaEl = antennaEl/180*pi;
beamwidth = beamwidth/180*pi;
switch antennaMode
    case 'SINC'                                         % ���˺���
       %F1 = abs(sinc(AzAngle)); 
        for m=1:length(AzAngle)
             F1(m) = sinc(2*(AzAngle(m)-antennaAz)/beamwidth); 
             for n = 1:length(ElAngle)               
                F2(n) = sinc(2*(ElAngle(n)-antennaEl)/beamwidth);
                subGain(m,n) = F1(m).*F2(n);  
            end
        end
    case '��˹'                                        % ��˹���� 
        F1 = exp(-1.4*(AzAngle-antennaAz).^2/beamwidth^2);
        F2 = exp(-1.4*(ElAngle-antennaEl).^2/beamwidth^2);
        subGain = F1.*F2;
    case '����'                                          % ���Һ���
        F1 = cos((AzAngle-antennaAz)*pi/(2*beamwidth));
        F2 = cos((ElAngle-antennaEl)*pi/(2*beamwidth));
        subGain = F1.*F2;
end
