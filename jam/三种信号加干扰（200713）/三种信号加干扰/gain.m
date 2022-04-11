function G = gain(radar,AzAngle,ElAngle)
%========================================================================%
%功能：计算方向图增益                                                     %
%输入：雷达参数radar，目标在雷达天线视场方位角AzAngle，俯仰角ElAngle       %
%输出：方向图增益Gainsum,GainAz,GainEl                                   %
%========================================================================%
% Gain = 40;                                          %天线最大增益，临时
% antennaMode = radar.antennaMode;                    %方向图类型
% switch antennaMode
%     case 'SINC'
%         antenaGain = abs(sinc(-pi:0.001:pi).^2);
%         s = trapz(antenaGain)*2*pi/length(antenaGain); %integral over antena gain
%         antenaGain = antenaGain/s;
%         GainAz = antenaGain(fix(pi/0.001+AzAngle/180*pi/0.001));
%         GainEl = antenaGain(fix(pi/0.001+ElAngle/180*pi/0.001));
%         G = Gain * GainAz * GainEl;      
% end
%% 参数初始化

Gain = 10^(radar.Gt/10);                                     %天线增益
antennaAzAngle = radar.antennaAzAngle;                      %天线方位指向
antennaELAngle = radar.antennaELAngle;                      %天线俯仰指向
AdjacentAngle = radar.AdjacentAngle;                         %子天线斜视角
beamwidth = radar.antennabeamwidth;
antennaMode = radar.strfunc{radar.nfunc};
% %% 画图用参数
% % AzAngle = -(AdjacentAngle+beamwidth):0.1:(AdjacentAngle+beamwidth);
% % ElAngle = -(AdjacentAngle+beamwidth):0.1:(AdjacentAngle+beamwidth);
% AzAngle = -3:0.1:3;
% ElAngle = -3:0.1:3;
%% 计算用
        antennaLU = subAntennaPattern(antennaMode,beamwidth,AzAngle,ElAngle,...
            antennaAzAngle-AdjacentAngle,antennaELAngle+AdjacentAngle); %左上子天线
        antennaLD = subAntennaPattern(antennaMode,beamwidth,AzAngle,ElAngle,...
            antennaAzAngle-AdjacentAngle,antennaELAngle-AdjacentAngle); %左下子天线
        antennaRU = subAntennaPattern(antennaMode,beamwidth,AzAngle,ElAngle,...
            antennaAzAngle+AdjacentAngle,antennaELAngle+AdjacentAngle); %右上子天线
        antennaRD = subAntennaPattern(antennaMode,beamwidth,AzAngle,ElAngle,...
            antennaAzAngle+AdjacentAngle,antennaELAngle-AdjacentAngle); %右下子天线
        G.Gainsum =(antennaRU + antennaRD + antennaLU + antennaLD);        %和波束
        G.GainAz = Gain*(antennaRU + antennaRD - antennaLU - antennaLD);         %方位差波束
        G.GainEl = Gain*(antennaRU - antennaRD + antennaLU - antennaLD);         %俯仰差波束
        G.NDSAz =G.GainAz./G.Gainsum;             %方位归一化差斜率normalization difference slop 
        G.NDSEl =G.GainEl./G.Gainsum;             %俯仰归一化差斜率normalization difference slop 

%% 画图
% figure(31);
% mesh(AzAngle*pi/180,ElAngle*pi/180,antennaLU);
% hold on;
% mesh(AzAngle*pi/180,ElAngle*pi/180,antennaLD);
% mesh(AzAngle*pi/180,ElAngle*pi/180,antennaRU);
% mesh(AzAngle*pi/180,ElAngle*pi/180,antennaRD);
% xlabel('\theta（^。）' );
% ylabel('\phi（^。）' );
% zlabel('幅度（^。）' );
% title('四个波束响应');
% grid on;

% figure(32);
% subplot(411);
% plot(AzAngle*pi/180,antennaLU(:,round(length(AzAngle)/2)),'r');
% hold on;
% plot(AzAngle*pi/180,antennaLU(round(length(ElAngle)/2),:),'b');
% xlabel('\theta（^。）' );
% ylabel('幅度（^。）' );
% title('两个波束响应');
% grid on;
% 
% subplot(412);
% plot(AzAngle*pi/180,G.Gainsum(:,round(length(AzAngle)/2)));
% grid on;
% xlabel('\theta（^。）' );
% title('和波束响应');
% 
% subplot(413);
% plot(AzAngle*pi/180,G.GainAz(:,round(length(AzAngle)/2)),'r');
% % hold on;
% % plot(ElAngle*pi/180,G.GainEl(round(length(ElAngle)/2),:),'b');
% grid on;
% xlabel('\theta（^。）' );
% title('差波束响应');
% 
% subplot(414); 
% plot(AzAngle*pi/180,G.NDSAz(:,round(length(AzAngle)/2)),'r');
% grid on;
% xlabel('\theta（^。）' );
% title('差波束/和波束');

%计算子方向图
function subGain = subAntennaPattern(antennaMode,beamwidth,AzAngle,ElAngle,antennaAz,antennaEl)
%========================================================================%
%功能：计算子方向图增益                                                   %
%输入：方向图类型antennaMode，波束宽度beamwidth，角度angle                 %
%输出：子方向图增益subGain                                                        %
%========================================================================%
%角度化成弧度
AzAngle = AzAngle/180*pi;
ElAngle = ElAngle/180*pi;
antennaAz = antennaAz/180*pi;
antennaEl = antennaEl/180*pi;
beamwidth = beamwidth/180*pi;
switch antennaMode
    case 'SINC'                                         % 辛克函数
       %F1 = abs(sinc(AzAngle)); 
        for m=1:length(AzAngle)
             F1(m) = sinc(2*(AzAngle(m)-antennaAz)/beamwidth); 
             for n = 1:length(ElAngle)               
                F2(n) = sinc(2*(ElAngle(n)-antennaEl)/beamwidth);
                subGain(m,n) = F1(m).*F2(n);  
            end
        end
    case '高斯'                                        % 高斯函数 
        F1 = exp(-1.4*(AzAngle-antennaAz).^2/beamwidth^2);
        F2 = exp(-1.4*(ElAngle-antennaEl).^2/beamwidth^2);
        subGain = F1.*F2;
    case '余弦'                                          % 余弦函数
        F1 = cos((AzAngle-antennaAz)*pi/(2*beamwidth));
        F2 = cos((ElAngle-antennaEl)*pi/(2*beamwidth));
        subGain = F1.*F2;
end
