%070228��Ӳ�λ���ݱ���
clear all
close all
clc
%�������,��λ����
BeamWidth = 3;  % ʾ����3.6
theta05 = round(sin(BeamWidth*pi/180)*1000)/1000;
%�״����ÿ���
AZ = [-60 60];              %��λ�Ƿ�Χ(�״�������ϵ)����λ���� ʾ����[-15 15]
EL = [-10 50];                %�����Ƿ�Χ(�״�������ϵ)����λ���� ʾ����[0 30]
%������ǣ���λ����
thetaT = atan(-(cos(EL(1)*pi/180) - cos(EL(2)*pi/180))/(sin(EL(1)*pi/180) - sin(EL(2)*pi/180))*cos(AZ(2)*pi/180))*180/pi;
%��߽�
az = AZ(1);
el = EL(1):EL(2);
alpha1 = cos(el*pi/180)*sin(az*pi/180);
beta1 = sin(el*pi/180)*cos(thetaT*pi/180)-cos(el*pi/180)*cos(az*pi/180)*sin(thetaT*pi/180);
%�ұ߽�
az = AZ(2);
el = EL(1):EL(2);
alpha2 = cos(el*pi/180)*sin(az*pi/180);
beta2 = sin(el*pi/180)*cos(thetaT*pi/180)-cos(el*pi/180)*cos(az*pi/180)*sin(thetaT*pi/180);
%�±߽�
az = AZ(1):AZ(2);
el = EL(1);
alpha3 = cos(el*pi/180)*sin(az*pi/180);
beta3 = sin(el*pi/180)*cos(thetaT*pi/180)-cos(el*pi/180)*cos(az*pi/180)*sin(thetaT*pi/180);
%�ϱ߽�
az = AZ(1):AZ(2);
el = EL(2);
alpha4 = cos(el*pi/180)*sin(az*pi/180);
beta4 = sin(el*pi/180)*cos(thetaT*pi/180)-cos(el*pi/180)*cos(az*pi/180)*sin(thetaT*pi/180);
lowup = max(beta3);
lowest = min(beta3);
rightest = max(alpha3) + theta05/2;
leftest = min(alpha3) - theta05/2;
upest = max(beta4);
uplow = min(beta4);
figure(1),plot(alpha1,beta1,'r');axis([ leftest-0.1 rightest+0.05 lowest-0.05 upest+0.1]);hold on;grid on;
title('\fontsize{10.5}\fontname{����}�״ﲨ�������ҿռ�ı���');
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 10.5);
xlabel('\fontname{����}��λ��\fontname{Times New Roman}/rad', 'FontSize', 10.5);
ylabel('\fontname{����}������\fontname{Times New Roman}/rad', 'FontSize', 10.5);
plot(alpha2,beta2,'r');
plot(alpha3,beta3,'r');
plot(alpha4,beta4,'r');
% pause(0.1)
%==========================================================================
BoweiNumber = 0;
Mx = 2*round(rightest/theta05)-1;
My = round((upest - lowest)/(0.866*theta05));
arr_x0 = zeros(My,Mx);  %��λ���꣨��λ�ǣ�
arr_y0 = zeros(My,Mx);  %��λ���꣨�����ǣ�
arr_NboweiRow = zeros(My,1); %ÿһ�еĲ�λ��Ŀ
%---------------������------------------------
m = 0;
for y0 = lowest+theta05/2: 2*0.866*theta05 : upest
    m = m + 1;
    n = 0;
    NboweiRow = 0;
    for x0 = 0:theta05:rightest
        n = n + 1;
        x = (-theta05/2:0.001:theta05/2) + x0;
        y1 = sqrt((theta05/2)^2 - (x - x0).^2) + y0;        
        y2 = -sqrt((theta05/2)^2 - (x - x0).^2) + y0;
        % alpha2,beta2:�ұ߽�    
        flag = 0;
        [minValue minSite] = min(abs(beta2 - y0));
        if alpha2(minSite) + theta05/2 >= x0
               flag = 1;
        end
        if y0 < lowup
            % alpha3,beta3:�±߽�       
            [minValue minSite] = min(abs(alpha3 - x0));
            if beta3(minSite) > y0
                flag = 0;
            end
        end
        if y0 > uplow
            % alpha4,beta4:�ϱ߽�       
            [minValue minSite] = min(abs(alpha4 - x0));
            if beta4(minSite) + theta05/2 <= y0
                flag = 0;
            end
        end
        if flag == 1
            NboweiRow = NboweiRow + 1; 
            plot(x,real(y1));
            plot(x,real(y2));            
            pause(0.01)
            BoweiNumber = BoweiNumber + 1;
            arr_x0((m-1)*2+1,round(Mx/2)+n) = x0;
            arr_y0((m-1)*2+1,round(Mx/2)+n) = y0;
        end
    end    
    arr_NboweiRow( (m-1)*2 + 1 ) = NboweiRow;
end
m = 0;
for y0 = lowest+theta05/2 : 2*0.866*theta05 : upest
    m = m + 1;
    n = 0;
    NboweiRow = 0;
    for x0 = -theta05:-theta05:leftest
        n = n + 1;
        x = (-theta05/2:0.001:theta05/2) + x0;
        y1 = sqrt((theta05/2)^2 - (x - x0).^2) + y0;        
        y2 = -sqrt((theta05/2)^2 - (x - x0).^2) + y0;
        % alpha1,beta1:��߽�   
        flag = 0;
        [minValue minSite] = min(abs(beta1 - y0));
        if alpha1(minSite) - theta05/2 <= x0
            flag = 1;
        end        
        if y0 < lowup
            % alpha3,beta3:�±߽�       
            [minValue minSite] = min(abs(alpha3 - x0));
            if beta3(minSite) > y0
                flag = 0;
            end
        end
       if y0 > uplow
            % alpha4,beta4:�ϱ߽�       
            [minValue minSite] = min(abs(alpha4 - x0));
            if beta4(minSite) + theta05/2 <= y0
                flag = 0;
            end
        end
        if flag == 1
            NboweiRow = NboweiRow + 1; 
            plot(x,real(y1));
            plot(x,real(y2));
            pause(0.01)
            BoweiNumber = BoweiNumber + 1;
            arr_x0((m-1)*2+1,round(Mx/2)-n+1) = x0;
            arr_y0((m-1)*2+1,round(Mx/2)-n+1) = y0;
        end
    end
    arr_NboweiRow( (m-1)*2 + 1 ) = arr_NboweiRow( (m-1)*2 + 1 ) + NboweiRow;
end
%----------------------ż����-------------------------------
m = 0;
for y0 = lowest+theta05/2+0.866*theta05 : 2*0.866*theta05 : upest
    m = m + 1;
    n = 0;
    NboweiRow = 0;
    for x0 = theta05/2:theta05:rightest
        n = n + 1;
        x = (-theta05/2:0.001:theta05/2) + x0;
        y1 = sqrt((theta05/2)^2 - (x - x0).^2) + y0;        
        y2 = -sqrt((theta05/2)^2 - (x - x0).^2) + y0;
        % alpha2,beta2:�ұ߽�   
        flag = 0;
        [minValue minSite] = min(abs(beta2 - y0));
        if alpha2(minSite)  + theta05/2 >= x0
            flag = 1;
        end          
        if y0 < lowup
            % alpha3,beta3:�±߽�       
            [minValue minSite] = min(abs(alpha3 - x0));
            if beta3(minSite) > y0
                flag = 0;
            end
        end
        if y0 > uplow
            % alpha4,beta4:�ϱ߽�       
            [minValue minSite] = min(abs(alpha4 - x0));
            if beta4(minSite) + theta05/2 <= y0
                flag = 0;
            end
        end
        if flag == 1
            NboweiRow = NboweiRow + 1; 
            plot(x,real(y1));
            plot(x,real(y2));
            pause(0.01)
            BoweiNumber = BoweiNumber + 1;
            arr_x0(m*2,round(Mx/2)+n+1) = x0;
            arr_y0(m*2,round(Mx/2)+n+1) = y0;
        end
    end
    arr_NboweiRow( m*2 ) =  NboweiRow;
end
m = 0;
for y0 = lowest+theta05/2+0.866*theta05 : 2*0.866*theta05 : upest
    m = m + 1;
    n = 0;
    NboweiRow = 0;
    for x0 = -theta05/2:-theta05:leftest
        n = n + 1;
        x = (-theta05/2:0.001:theta05/2) + x0;
        y1 = sqrt((theta05/2)^2 - (x - x0).^2) + y0;        
        y2 = -sqrt((theta05/2)^2 - (x - x0).^2) + y0;       
        % alpha1,beta1:��߽�  
        flag = 0;
        [minValue minSite] = min(abs(beta1 - y0));
        if alpha1(minSite) - theta05/2 <= x0
            flag = 1;
        end  
        if y0 < lowup
            % alpha3,beta3:�±߽�       
            [minValue minSite] = min(abs(alpha3 - x0));
            if beta3(minSite) > y0
                flag = 0;
            end
        end
        if y0 > uplow
            % alpha4,beta4:�ϱ߽�       
            [minValue minSite] = min(abs(alpha4 - x0));
            if beta4(minSite) + theta05/2 <= y0
                flag = 0;
            end
        end        
        if flag == 1
            NboweiRow = NboweiRow + 1;
            plot(x,real(y1));
            plot(x,real(y2));
            pause(0.01)
            BoweiNumber = BoweiNumber + 1;
            arr_x0(m*2,round(Mx/2)-n+1) = x0;
            arr_y0(m*2,round(Mx/2)-n+1) = y0;
        end        
    end
    arr_NboweiRow( m*2 ) = arr_NboweiRow( m*2 ) + NboweiRow;
end
%-----------------------------------------------------
% figure(2),
% [b,a] = size(arr_x0);
% for mm = 1:b;
%     for nn = 1:a
%         plot(arr_x0(mm,nn),arr_y0(mm,nn),'o','LineWidth',2,...
%                 'MarkerEdgeColor','y',...
%                 'MarkerFaceColor','g',...
%                 'MarkerSize',10);
%             axis([leftest rightest lowest upest]);hold on;            
%     end
% end
% title('��λʾ��ͼ');
BoweiNumber = BoweiNumber
%-----------------------------------------------------
%�״�վ����ϵ�µĲ����ֲ�
theta = thetaT*pi/180;
[row,collum] = size(arr_x0);
arr_alpha = zeros(row,collum);
arr_beta = zeros(row,collum);
for mm = 1:row
    for nn = 1:collum
        Tx = arr_x0(mm,nn);
        Ty = arr_y0(mm,nn);
        temp = sqrt(1 - Tx^2 - Ty^2);
        arr_beta(mm,nn) = asin(Ty*cos(theta) + temp*sin(theta));
        arr_alpha(mm,nn) = asin(Tx/cos(arr_beta(mm,nn)));
    end
end
arr_beta = arr_beta*180/pi;         %��λ��λ��
arr_alpha = arr_alpha*180/pi ;      %��λ������
arr_NboweiRow;                      %ÿһ�еĲ�λ����
return;
%==========================================================================