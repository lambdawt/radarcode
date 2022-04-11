function [s_echo_mf]=jianbo(s_echo_1,N,frame,f0,fs)
N2=10;            %10位AD
Vmax=max(max(abs(s_echo_1)));
ad1=Vmax/2.^N2;
ad2=round(s_echo_1.*2^N2/Vmax);  %取整
s_echo_1=ad1.*ad2;

s_echo=reshape(s_echo_1,N,frame);
s_echo=s_echo.';
% f0=f1;
n=0:N-1;
% local_oscillator_i=2*cos(n*(-f0)/fs*2*pi);
% local_oscillator_q=2*sin(n*(-f0)/fs*2*pi);
local_oscillator_i=cos(n*(-f0)/fs*2*pi);
local_oscillator_q=sin(n*(-f0)/fs*2*pi);
% window=chebwin(51,40);
% [b,a]=fir1(50,2*B/fs,window);
s_echo_i=zeros(frame,N);
s_echo_q=zeros(frame,N);
for i=1:frame
    s_echo_i(i,:)=local_oscillator_i.*s_echo(i,:);
%     s1(i,:)=[s_echo_i(i,:),zeros(1,25)];
% %     s1=[s1,zeros(1,25)];
%     s_echo_i(i,:)=[s_echo_i(:,i),zeros(1,25)];
    s_echo_q(i,:)=local_oscillator_q.*s_echo(i,:);
%    s2(i,:)=[s_echo_q(i,:),zeros(1,25)];
%     s2=[s2,zeros(1,25)];
%     s_echo_q(:,i)=local_oscillator_q.*s_echo(:,i);
%     s2=s_echo_q(:,i).';
%     s2=[s2,zeros(1,25)];
%     s_echo_q(:,i)=[s_echo_q(:,i),zeros(1,25)'];
%     s_echo_q(:,i)=s_echo_q(:,i)
%     s1(i,:)=filter(b,a,s1(i,:));
%     s2(i,:)=filter(b,a,s2(i,:));
%     s_echo_i(i,:)=s1(i,26:end);
%     s_echo_q(i,:)=s2(i,26:end);
% N2=10;            %10位AD
% Vmax=max(max(abs(s_echo_i(:,i))));%3.5e-6;          %电压5v
% ad1=Vmax/2.^N2;
% ad2=round(s_echo_i(:,i).*2^N2/Vmax);  %取整
% s_echo_i(:,i)=ad1.*ad2;
%     s_echo_q(:,i)=filter(b,a,s2.');
%     s_echo_q(:,i)=s_echo_q(26:end,i)
% N2=10;            %10位AD
% Vmax=max(max(abs(s_echo_q(:,i))));          %电压5v
% ad1=Vmax/2.^N2;
% ad2=round(s_echo_q(:,i).*2^N2/Vmax);  %取整
% s_echo_q(:,i)=ad1.*ad2;
    s_echo_mf(i,:)=s_echo_i(i,:)+1j*s_echo_q(i,:);
end
s_echo_mf=s_echo_mf.';

% figure,plot(0:ts:length(s_echo_mf(:,1))*ts-ts,s_echo_mf(:,1)),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('调制后第一个回波');
% figure,plot((0:fs/length(abs(fft(s_echo_mf(:,1)))):fs-fs/length(abs(fft(s_echo_mf(:,1))))),abs(fft(s_echo_mf(:,1)))),xlabel('f(单位：Hz)'),ylabel('y'),title('调制后回波信号的频域');
s_echo_result=reshape(s_echo_mf,1,frame*N);

