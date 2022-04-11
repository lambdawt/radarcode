function varargout = DeceptionJammingParameter(varargin)
% DECEPTIONJAMMINGPARAMETER MATLAB code for DeceptionJammingParameter.fig
%      DECEPTIONJAMMINGPARAMETER, by itself, creates a new DECEPTIONJAMMINGPARAMETER or raises the existing
%      singleton*.
%
%      H = DECEPTIONJAMMINGPARAMETER returns the handle to a new DECEPTIONJAMMINGPARAMETER or the handle to
%      the existing singleton*.
%
%      DECEPTIONJAMMINGPARAMETER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DECEPTIONJAMMINGPARAMETER.M with the given input arguments.
%
%      DECEPTIONJAMMINGPARAMETER('Property','Value',...) creates a new DECEPTIONJAMMINGPARAMETER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DeceptionJammingParameter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DeceptionJammingParameter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DeceptionJammingParameter

% Last Modified by GUIDE v2.5 02-Nov-2016 22:19:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DeceptionJammingParameter_OpeningFcn, ...
                   'gui_OutputFcn',  @DeceptionJammingParameter_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before DeceptionJammingParameter is made visible.
function DeceptionJammingParameter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DeceptionJammingParameter (see VARARGIN)

% Choose default command line output for DeceptionJammingParameter
handles.output = hObject;

global sigma0;
sigma0 = 2;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DeceptionJammingParameter wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DeceptionJammingParameter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global RadarS;
global sigma;
global rcsk;
rcsk = 2;
global sigma0;
sigma = rcs(rcsk,sigma0);
global string1;
global string2;
global temp1;
global temp2;
c=3e8;
Gj=str2num(get(handles.Gj,'string'));
Gjr=str2num(get(handles.Gjr,'string'));
Pj=str2num(get(handles.Pj,'string'));
tf=str2num(get(handles.tf,'string'));
vf=str2num(get(handles.vf,'string'));
L=3;%综合损耗
load data/data_target0Parameter
if RadarS == 1
string1=get(handles.n1,'string');
string2=get(handles.n2,'string');
temp1=get(handles.n1,'Value');
temp2=get(handles.n2,'Value');
load data/data_LFMParameter;
switch string1{temp1}
    case '速度'
        switch string2{temp2}
            case {'单' ,'多','密集' }
               R1=str2num(get(handles.R1,'string'));
               v1=str2num(get(handles.v1,'string'));
               congmubiao=str2num(get(handles.congmubiao,'string'));
               fr=1/tr_LFM;%脉冲重复频
                lamta=c/fz_LFM;%波长
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%一个脉冲重复周期采样序列
                N=length(tm);%一个脉冲重复周期采样点数长度
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%噪声强度
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %#ok<*NODEF> %目标回波信号功率
                A=sqrt(Prs);%回波信号幅度
                f_doppler=2*v/lamta;%真目标多普勒频率%线性调频信号调制系数
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
                Aj=sqrt(Prj);
                %速度
                f_doppler1=2*v1/lamta;
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [~,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y,tr_LFM,ts,N);
                [s_echo_2,echo]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                [s_ft,echo3]=LFMDeceptionJamming(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao,temp1);
                s_echo_1=s_echo_2+s_noise+s_ft;%%%%干扰信号
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('雷达接收信号');
%                 figure,plot(linspace(1.2e-5,1.3e-5,N),real(s_echo_1(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('雷达接收信号');
                figure,plot((0:fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_LFM-fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('雷达接收信号的频谱');
                save data/data_DeceptionJammingParameter R1 v1 congmubiao Gj Gjr Pj   L
            case '拖引'
%                  load data/data_DeceptionJammingParameter
                tf=10;
                c=3e8;
                lamta=c/fz_LFM;
                Rj=2e3;
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L_LFM/10)));
                Aj=sqrt(Prj);
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%一个脉冲重复周期采样序列
                N=length(tm);%一个脉冲重复周期采样点数长度
                f_doppler=2*v/lamta;
                f_doppler1=f_doppler;
                k=B1_LFM/tau_LFM;%线性调频信号调制系数
                ts=1/fs_LFM;%采样间隔
                %生成发射信号
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L/10))); %目标回波信号功率
                A=sqrt(Prs);
                Te=1143;%温度
                F=4;%噪声系数
                An=10*log10((1.382e-23)*Te*B_LFM*10^(F/10));%噪声强度
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k); 
                [~,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y,tr_LFM,ts,N);
                %生成回波信号
                [s_echo_2,echo]=LFMtuoyinhuiboxinhao(y,R,tr_LFM,ts,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                %figure,plot(0:ts:(N-1)*ts,real(s_echo_2(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('回波信号');
                %figure,plot((0:fs/length(abs(fft(abs(fftshift(echo(1,:)))))):fs-fs/length(abs(fft(abs(fftshift(echo(1,:))))))),abs((fft(echo(1,:))))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('回波信号的频谱');

                %生成干扰信号
                [s_ft,echo3]=LFMtuoyinganrao(R,frame_LFM,tf,Aj,tm,tau_LFM,f0_LFM,B1_LFM,N,f_doppler1,f_doppler,vf,lamta,k,ts,c,fs_LFM,temp1);

                %生成噪声
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                %目标回波信号、假目标信号、噪声叠加在一起送入接收机
                s_echo_1=s_echo_2+s_noise+s_ft;
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('雷达接收信号');
                figure,plot((0:fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_LFM-fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('雷达接收信号的频谱');
                save data/data_DeceptionJammingParameter  Gj Gjr Pj tf vf L
        end
    case '距离'
         switch string2{temp2}
            case {'单' ,'多','密集' }
                R1=str2num(get(handles.R1,'string'));
                v1=str2num(get(handles.v1,'string'));
                congmubiao=str2num(get(handles.congmubiao,'string'));
                fr=1/tr_LFM;%脉冲重复频
                lamta=c/fz_LFM;%波长
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%一个脉冲重复周期采样序列
                N=length(tm);%一个脉冲重复周期采样点数长度
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%噪声强度
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %目标回波信号功率
                A=sqrt(Prs);%回波信号幅度
                f_doppler=2*v/lamta;%真目标多普勒频率%线性调频信号调制系数
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
                Aj=sqrt(Prj);
                %速度
                f_doppler1=2*v1/lamta;
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [~,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y,tr_LFM,ts,N);
                [s_echo_2,echo]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                [s_ft,echo3]=LFMDeceptionJamming(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao,temp1);
                s_echo_1=s_echo_2+s_noise+s_ft;%%%%干扰信号
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('雷达接收信号');
                figure,plot((0:fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_LFM-fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('雷达接收信号的频谱');
                save data/data_DeceptionJammingParameter R1 v1 congmubiao Gj Gjr Pj tf vf L
           
            case '拖引'
                load data/data_target0Parameter
                c=3e8;
                lamta=c/fz_LFM;
                Rj=2e3;
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L_LFM/10)));
                Aj=sqrt(Prj);
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%一个脉冲重复周期采样序列
                N=length(tm);%一个脉冲重复周期采样点数长度
                f_doppler=2*v/lamta;
                f_doppler1=f_doppler;
                k=B1_LFM/tau_LFM;%线性调频信号调制系数
                ts=1/fs_LFM;%采样间隔
                %生成发射信号
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L/10))); %目标回波信号功率
                A=sqrt(Prs);
                Te=1143;%温度
                F=4;%噪声系数
                An=10*log10((1.382e-23)*Te*B_LFM*10^(F/10));%噪声强度
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k); 
                [~,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y,tr_LFM,ts,N);
                %生成回波信号
                [s_echo_2,echo]=LFMtuoyinhuiboxinhao(y,R,tr_LFM,ts,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                %figure,plot(0:ts:(N-1)*ts,real(s_echo_2(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('回波信号');
                %figure,plot((0:fs/length(abs(fft(abs(fftshift(echo(1,:)))))):fs-fs/length(abs(fft(abs(fftshift(echo(1,:))))))),abs((fft(echo(1,:))))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('回波信号的频谱');

                %生成干扰信号
                [s_ft,echo3]=LFMtuoyinganrao(R,frame_LFM,tf,Aj,tm,tau_LFM,f0_LFM,B1_LFM,N,f_doppler1,f_doppler,vf,lamta,k,ts,c,fs_LFM,temp1);

                %生成噪声
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                %目标回波信号、假目标信号、噪声叠加在一起送入接收机
                s_echo_1=s_echo_2+s_noise+s_ft;
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('雷达接收信号');
                figure,plot((0:fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_LFM-fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('雷达接收信号的频谱');
                save data/data_DeceptionJammingParameter  Gj Gjr Pj tf vf L
         end
    case '联合'
         switch string2{temp2}
            case {'单' ,'多','密集' }
                 R1=str2num(get(handles.R1,'string'));
                v1=str2num(get(handles.v1,'string'));
                congmubiao=str2num(get(handles.congmubiao,'string'));
                fr=1/tr_LFM;%脉冲重复频
                lamta=c/fz_LFM;%波长
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%一个脉冲重复周期采样序列
                N=length(tm);%一个脉冲重复周期采样点数长度
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%噪声强度
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %目标回波信号功率
                A=sqrt(Prs);%回波信号幅度
                f_doppler=2*v/lamta;%真目标多普勒频率%线性调频信号调制系数
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
                Aj=sqrt(Prj);
                %速度
                f_doppler1=2*v1/lamta;
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [~,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y,tr_LFM,ts,N);
                [s_echo_2,echo]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                [s_ft,echo3]=LFMDeceptionJamming(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao,temp1);
                s_echo_1=s_echo_2+s_noise+s_ft;%%%%干扰信号
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('雷达接收信号');
                figure,plot((0:fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_LFM-fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('雷达接收信号的频谱');
                save data/data_DeceptionJammingParameter R1 v1 congmubiao Gj Gjr Pj tf vf L
            
            case '拖引'
                load data/data_target0Parameter
                c=3e8;
                lamta=c/fz_LFM;
                Rj=2e3;
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L_LFM/10)));
                Aj=sqrt(Prj);
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%一个脉冲重复周期采样序列
                N=length(tm);%一个脉冲重复周期采样点数长度
                f_doppler=2*v/lamta;
                f_doppler1=f_doppler;
                k=B1_LFM/tau_LFM;%线性调频信号调制系数
                ts=1/fs_LFM;%采样间隔
                %生成发射信号
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L/10))); %目标回波信号功率
                A=sqrt(Prs);
                Te=1143;%温度
                F=4;%噪声系数
                An=10*log10((1.382e-23)*Te*B_LFM*10^(F/10));%噪声强度
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k); 
                [~,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y,tr_LFM,ts,N);
                %生成回波信号
                [s_echo_2,echo]=LFMtuoyinhuiboxinhao(y,R,tr_LFM,ts,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                %figure,plot(0:ts:(N-1)*ts,real(s_echo_2(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('回波信号');
                %figure,plot((0:fs/length(abs(fft(abs(fftshift(echo(1,:)))))):fs-fs/length(abs(fft(abs(fftshift(echo(1,:))))))),abs((fft(echo(1,:))))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('回波信号的频谱');

                %生成干扰信号
                [s_ft,echo3]=LFMtuoyinganrao(R,frame_LFM,tf,Aj,tm,tau_LFM,f0_LFM,B1_LFM,N,f_doppler1,f_doppler,vf,lamta,k,ts,c,fs_LFM,temp1);

                %生成噪声
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                %目标回波信号、假目标信号、噪声叠加在一起送入接收机
                s_echo_1=s_echo_2+s_noise+s_ft;
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('雷达接收信号');
                figure,plot((0:fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_LFM-fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('雷达接收信号的频谱');
                save data/data_DeceptionJammingParameter  Gj Gjr Pj tf vf L
         end
         
end
elseif RadarS==2
                        
      string1=get(handles.n1,'string');
      string2=get(handles.n2,'string');
      temp1=get(handles.n1,'Value');
      temp2=get(handles.n2,'Value');
      load data/data_BFParameter;
      global code; %#ok<*TLEV>
     switch string1{temp1}
    case '速度'
        switch string2{temp2}
        case {'单' ,'多','密集' }
        R1=str2num(get(handles.R1,'string'));
        v1=str2num(get(handles.v1,'string'));
        congmubiao=str2num(get(handles.congmubiao,'string'));
        tm_B=0:1/fs_B:tr_B-1/fs_B;%一个脉冲重复周期采样序列
        N=length(tm_B);%一个脉冲重复周期采样点数长度
        number1=length(code);
        ts=1/fs_B;
        lamta=c/fz_B;
        fr=1/tr_B;
        An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%噪声强度
        Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %目标回波信号功率
        A=sqrt(Prs);%回波信号幅度
        f_doppler=2*v/lamta;%真目标多普勒频率
        Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
        Aj=sqrt(Prj);
        %速度
        f_doppler1=2*v1/lamta; 
        [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
        [s_echo_2,echo]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
        [~,match_filter_fft]=maiyaxishu(f0_B,fs_B,y,tr_B,ts,N);
        [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
        [s_ft,echo3]=BKDeceptionJamming(D,y,R1,tr_B,ts,c,Aj,N,frame_B,fs_B,f_doppler1,tm_B,f0_B,tau_B,congmubiao,y1,temp1);
        s_echo_1=s_echo_2+s_noise+s_ft;%%%%干扰信号
        figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('雷达接收信号');
        figure,plot((0:fs_B/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_B-fs_B/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('雷达接收信号的频谱');
        save data/data_DeceptionJammingParameter R1 v1 congmubiao Gj Gjr Pj tf vf L
        case '拖引'
            c=3e8;
            Rj=2e3;
            number1=length(code);
            ts=1/fs_B;
            lamta=c/fz_B;
            tm=0:1/fs_B:tr_B-1/fs_B;%一个脉冲重复周期采样序列
            N=length(tm);%一个脉冲重复周期采样点数长度
            f_doppler=2*v/lamta;%真目标多普勒频率
            f_doppler1= f_doppler;%真目标多普勒频率
            Te=1143;%温度
            An=10*log10((1.382e-23)*Te*B_B*10^(F_B/10));%噪声强度
            Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %目标回波信号功率
            A=sqrt(Prs);%回波信号幅度
            Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L/10)));
            Aj=sqrt(Prj);
            [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
%             figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(单位：秒)'), ylabel('y(单位：伏)'),title('发射信号');
%             figure,plot((0:fs_B/length(abs(fft(abs(fftshift(y))))):fs_B-fs_B/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('发射信号的频谱');
%             %%%%%%%%%%%1.1.2生成理想脉冲压缩系数
            [~,match_filter_fft]=maiyaxishu(f0_B,fs_B,y,tr_B,ts,N);
             %%%%%%%%%%%%%1.3生成回波信号%%%%%%%%%%%%
            [s_echo_2,echo]=BKtuoyinhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tm,tau_B); 
            %生成噪声
            [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
            %干扰信号
            [s_ft,echo3]=BKtuoyinganrao(R,frame_B,tf,Aj,tau_B,N,f_doppler1,f_doppler,vf,lamta,ts,c,fs_B,y1,tr_B,temp1);
            %目标回波信号、假目标信号、噪声叠加在一起送入接收机
            s_echo_1=s_echo_2+s_noise+s_ft;
            figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('雷达接收信号');
            figure,plot((0:fs_B/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_B-fs_B/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('雷达接收信号的频谱');
            save data/data_DeceptionJammingParameter  Gj Gjr Pj tf vf L
        end
  case '距离'
        switch string2{temp2}
        case {'单' ,'多','密集' }
        R1=str2num(get(handles.R1,'string'));
        v1=str2num(get(handles.v1,'string'));
        congmubiao=str2num(get(handles.congmubiao,'string'));
        load data/data_BFParameter;
        tm_B=0:1/fs_B:tr_B-1/fs_B;%一个脉冲重复周期采样序列
        N=length(tm_B);%一个脉冲重复周期采样点数长度
        number1=length(code);
        ts=1/fs_B;
        lamta=c/fz_B;
        fr=1/tr_B;
        An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%噪声强度
        Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %目标回波信号功率
        A=sqrt(Prs);%回波信号幅度
        f_doppler=2*v/lamta;%真目标多普勒频率
        Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
        Aj=sqrt(Prj);
        %速度
        f_doppler1=2*v1/lamta; 
        [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
        [s_echo_2,echo]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
        [~,match_filter_fft]=maiyaxishu(f0_B,fs_B,y,tr_B,ts,N);
        [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
        [s_ft,echo3]=BKDeceptionJamming(D,y,R1,tr_B,ts,c,Aj,N,frame_B,fs_B,f_doppler1,tm_B,f0_B,tau_B,congmubiao,y1,temp1);
        s_echo_1=s_echo_2+s_noise+s_ft;%%%%干扰信号
        figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('雷达接收信号');
        figure,plot((0:fs_B/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_B-fs_B/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('雷达接收信号的频谱');
        save data/data_DeceptionJammingParameter R1 v1 congmubiao Gj Gjr Pj tf vf L
        case '拖引'
           c=3e8;
            Rj=2e3;
            number1=length(code);
            ts=1/fs_B;
            lamta=c/fz_B;
            tm=0:1/fs_B:tr_B-1/fs_B;%一个脉冲重复周期采样序列
            N=length(tm);%一个脉冲重复周期采样点数长度
            f_doppler=2*v/lamta;%真目标多普勒频率
            f_doppler1= f_doppler;%真目标多普勒频率
            Te=1143;%温度
            An=10*log10((1.382e-23)*Te*B_B*10^(F_B/10));%噪声强度
            Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %目标回波信号功率
            A=sqrt(Prs);%回波信号幅度
            Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L/10)));
            Aj=sqrt(Prj);
            [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
%             figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(单位：秒)'), ylabel('y(单位：伏)'),title('发射信号');
%             figure,plot((0:fs_B/length(abs(fft(abs(fftshift(y))))):fs_B-fs_B/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('发射信号的频谱');
%             %%%%%%%%%%%1.1.2生成理想脉冲压缩系数
            [~,match_filter_fft]=maiyaxishu(f0_B,fs_B,y,tr_B,ts,N);
             %%%%%%%%%%%%%1.3生成回波信号%%%%%%%%%%%%
            [s_echo_2,echo]=BKtuoyinhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tm,tau_B); 
            %生成噪声
            [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
            %干扰信号
            [s_ft,echo3]=BKtuoyinganrao(R,frame_B,tf,Aj,tau_B,N,f_doppler1,f_doppler,vf,lamta,ts,c,fs_B,y1,tr_B,temp1);
            %目标回波信号、假目标信号、噪声叠加在一起送入接收机
            s_echo_1=s_echo_2+s_noise+s_ft;
            figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('雷达接收信号');
            figure,plot((0:fs_B/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_B-fs_B/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('雷达接收信号的频谱');
            save data/data_DeceptionJammingParameter  Gj Gjr Pj tf vf L
        end
       case '联合'
        switch string2{temp2}
        case {'单' ,'多','密集' }
        R1=str2num(get(handles.R1,'string'));
        v1=str2num(get(handles.v1,'string'));
        congmubiao=str2num(get(handles.congmubiao,'string'));
        load data/data_BFParameter;
        tm_B=0:1/fs_B:tr_B-1/fs_B;%一个脉冲重复周期采样序列
        N=length(tm_B);%一个脉冲重复周期采样点数长度
        number1=length(code);
        ts=1/fs_B;
        lamta=c/fz_B;
        fr=1/tr_B;
        An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%噪声强度
        Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %目标回波信号功率
        A=sqrt(Prs);%回波信号幅度
        f_doppler=2*v/lamta;%真目标多普勒频率
        Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
        Aj=sqrt(Prj);
        %速度
        f_doppler1=2*v1/lamta; 
        [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
        [s_echo_2,echo]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
        [~,match_filter_fft]=maiyaxishu(f0_B,fs_B,y,tr_B,ts,N);
        [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
        [s_ft,echo3]=BKDeceptionJamming(D,y,R1,tr_B,ts,c,Aj,N,frame_B,fs_B,f_doppler1,tm_B,f0_B,tau_B,congmubiao,y1,temp1);
        s_echo_1=s_echo_2+s_noise+s_ft;%%%%干扰信号
        figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('雷达接收信号');
        figure,plot((0:fs_B/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_B-fs_B/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('雷达接收信号的频谱');
        save data/data_DeceptionJammingParameter R1 v1 congmubiao Gj Gjr Pj tf vf L
      case '拖引'
           c=3e8;
            Rj=2e3;
            number1=length(code);
            ts=1/fs_B;
            lamta=c/fz_B;
            tm=0:1/fs_B:tr_B-1/fs_B;%一个脉冲重复周期采样序列
            N=length(tm);%一个脉冲重复周期采样点数长度
            f_doppler=2*v/lamta;%真目标多普勒频率
            f_doppler1= f_doppler;%真目标多普勒频率
            Te=1143;%温度
            An=10*log10((1.382e-23)*Te*B_B*10^(F_B/10));%噪声强度
            Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %目标回波信号功率
            A=sqrt(Prs);%回波信号幅度
            Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L/10)));
            Aj=sqrt(Prj);
            [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
%             figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(单位：秒)'), ylabel('y(单位：伏)'),title('发射信号');
%             figure,plot((0:fs_B/length(abs(fft(abs(fftshift(y))))):fs_B-fs_B/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('发射信号的频谱');
%             %%%%%%%%%%%1.1.2生成理想脉冲压缩系数
            [~,match_filter_fft]=maiyaxishu(f0_B,fs_B,y,tr_B,ts,N);
             %%%%%%%%%%%%%1.3生成回波信号%%%%%%%%%%%%
            [s_echo_2,echo]=BKtuoyinhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tm,tau_B); 
            %生成噪声
            [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
            %干扰信号
            [s_ft,echo3]=BKtuoyinganrao(R,frame_B,tf,Aj,tau_B,N,f_doppler1,f_doppler,vf,lamta,ts,c,fs_B,y1,tr_B,temp1);
            %目标回波信号、假目标信号、噪声叠加在一起送入接收机
            s_echo_1=s_echo_2+s_noise+s_ft;
            figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('雷达接收信号');
            figure,plot((0:fs_B/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_B-fs_B/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('雷达接收信号的频谱');
            save data/data_DeceptionJammingParameter  Gj Gjr Pj tf vf L
        end
     end
elseif RadarS==3
     string1=get(handles.n1,'string');
      string2=get(handles.n2,'string');
      temp1=get(handles.n1,'Value');
      temp2=get(handles.n2,'Value');
      load data/data_JDParameter;
     switch string1{temp1}
     case '速度'
        switch string2{temp2}
        case {'单' ,'多','密集' }
        R1=str2num(get(handles.R1,'string'));
        v1=str2num(get(handles.v1,'string'));
        congmubiao=str2num(get(handles.congmubiao,'string'));
        fr=1/tr_JD;
        lamta=c/fz_JD;%波长
        tm=0:1/fs_JD:tr_JD-1/fs_JD;%一个脉冲重复周期采样序列
        N=length(tm);%一个脉冲重复周期采样点数长度
        An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%噪声强度
        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        ts=1/fs_JD; 
        Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %目标回波信号功率
        A=sqrt(Prs);%回波信号幅度
        f_doppler=2*v/lamta;%真目标多普勒频率
        Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
        Aj=sqrt(Prj);
        %速度
         f_doppler1=2*v1/lamta; 
        [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
        [~,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y,tr_JD,ts,N);
        [s_echo_2,echo]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        [s_ft,echo3]=JDDeceptionJamming(D,R1,tr_JD,c,Aj,N,frame_JD,fs_JD,f_doppler1,tm,f0_JD,tau_JD,congmubiao,temp1);
        s_echo_1=s_echo_2+s_noise+s_ft;%%%%干扰信号
        figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('雷达接收信号');
        figure,plot((0:fs_JD/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_JD-fs_JD/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('雷达接收信号的频谱');
        save data/data_DeceptionJammingParameter R1 v1 congmubiao Gj Gjr Pj tf vf L
        
        case '拖引'
             %生成发射信号
        Rj=2e3;
        c=3e8;
        lamta=c/fz_JD;%波长
        tm=0:1/fs_JD:tr_JD-1/fs_JD;%一个脉冲重复周期采样序列
        N=length(tm);%一个脉冲重复周期采样点数长度
        An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%噪声强度
        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        ts=1/fs_JD; 
        Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %目标回波信号功率
        A=sqrt(Prs);%回波信号幅度
        f_doppler=2*v/lamta;%真目标多普勒频率
   
        %速度
        f_doppler=2*v/lamta;
        f_doppler1=f_doppler;
        lamta=c/fz_JD;
        f_doppler=2*v/lamta;
        [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
        ts=1/fs_JD;%采样间隔
        tm=0:1/fs_JD:tr_JD-1/fs_JD;%一个脉冲重复周期采样序列
        N=length(tm);%一个脉冲重复周期采样点数长度
      
        Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L/10)));
        Aj=sqrt(Prj);
%         figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(单位：秒)'), ylabel('y(单位：伏)'),title('发射信号');
%         figure,plot((0:fs/length(abs(fft(abs(fftshift(y))))):fs-fs/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('发射信号的频谱');
        %生成理想脉冲压缩系数
        [~,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y,tr_JD,ts,N);

        %生成回波信号
        [s_echo_2,echo]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
             [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        %生成干扰信号
        [s_ft,echo3]=JDtuoyinganrao(R,frame_JD,tf,Aj,tm,tau_JD,f0_JD,N,f_doppler1,f_doppler,vf,lamta,ts,c,fs_JD,temp1);

        %目标回波信号、假目标信号、噪声叠加在一起送入接收机
        s_echo_1=s_echo_2+s_noise+s_ft;
        figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('雷达接收信号');
        figure,plot((0:fs_JD/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_JD-fs_JD/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('雷达接收信号的频谱');
        save data/data_DeceptionJammingParameter  Gj Gjr Pj tf vf L
        
        end
    case '距离'
        switch string2{temp2}
        case {'单' ,'多','密集' }
        R1=str2num(get(handles.R1,'string'));
        v1=str2num(get(handles.v1,'string'));
        congmubiao=str2num(get(handles.congmubiao,'string'));
        load data/data_JDParameter;
        fr=1/tr_JD;
        lamta=c/fz_JD;%波长
        tm=0:1/fs_JD:tr_JD-1/fs_JD;%一个脉冲重复周期采样序列
        N=length(tm);%一个脉冲重复周期采样点数长度
        An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%噪声强度
        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        ts=1/fs_JD; 
        Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %目标回波信号功率
        A=sqrt(Prs);%回波信号幅度
        f_doppler=2*v/lamta;%真目标多普勒频率
        Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
        Aj=sqrt(Prj);
        %速度
         f_doppler1=2*v1/lamta; 
        [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
        [~,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y,tr_JD,ts,N);
        [s_echo_2,echo]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        [s_ft,echo3]=JDDeceptionJamming(D,R1,tr_JD,c,Aj,N,frame_JD,fs_JD,f_doppler1,tm,f0_JD,tau_JD,congmubiao,temp1);
        s_echo_1=s_echo_2+s_noise+s_ft;%%%%干扰信号
        figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('雷达接收信号');
        figure,plot((0:fs_JD/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_JD-fs_JD/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('雷达接收信号的频谱');
        save data/data_DeceptionJammingParameter R1 v1 congmubiao Gj Gjr Pj tf vf L
        
        case '拖引'
             %生成发射信号
        Rj=2e3;
        c=3e8;
        lamta=c/fz_JD;%波长
        tm=0:1/fs_JD:tr_JD-1/fs_JD;%一个脉冲重复周期采样序列
        N=length(tm);%一个脉冲重复周期采样点数长度
        An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%噪声强度
        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        ts=1/fs_JD; 
        Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %目标回波信号功率
        A=sqrt(Prs);%回波信号幅度
        f_doppler=2*v/lamta;%真目标多普勒频率
   
        %速度
        f_doppler=2*v/lamta;
        f_doppler1=f_doppler;
        lamta=c/fz_JD;
        f_doppler=2*v/lamta;
        [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
        ts=1/fs_JD;%采样间隔
        tm=0:1/fs_JD:tr_JD-1/fs_JD;%一个脉冲重复周期采样序列
        N=length(tm);%一个脉冲重复周期采样点数长度
      
        Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L/10)));
        Aj=sqrt(Prj);
%         figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(单位：秒)'), ylabel('y(单位：伏)'),title('发射信号');
%         figure,plot((0:fs/length(abs(fft(abs(fftshift(y))))):fs-fs/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('发射信号的频谱');
        %生成理想脉冲压缩系数
        [~,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y,tr_JD,ts,N);

        %生成回波信号
        [s_echo_2,echo]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
             [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        %生成干扰信号
        [s_ft,echo3]=JDtuoyinganrao(R,frame_JD,tf,Aj,tm,tau_JD,f0_JD,N,f_doppler1,f_doppler,vf,lamta,ts,c,fs_JD,temp1);

        %目标回波信号、假目标信号、噪声叠加在一起送入接收机
        s_echo_1=s_echo_2+s_noise+s_ft;
        figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('雷达接收信号');
        figure,plot((0:fs_JD/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_JD-fs_JD/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('雷达接收信号的频谱');
        save data/data_DeceptionJammingParameter  Gj Gjr Pj tf vf L
        
        end
        case '联合'
        switch string2{temp2}
        case {'单' ,'多','密集' }
        R1=str2num(get(handles.R1,'string'));
        v1=str2num(get(handles.v1,'string'));
        congmubiao=str2num(get(handles.congmubiao,'string'));
        load data/data_JDParameter;
        fr=1/tr_JD;
        lamta=c/fz_JD;%波长
        tm=0:1/fs_JD:tr_JD-1/fs_JD;%一个脉冲重复周期采样序列
        N=length(tm);%一个脉冲重复周期采样点数长度
        An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%噪声强度
        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        ts=1/fs_JD; 
        Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %目标回波信号功率
        A=sqrt(Prs);%回波信号幅度
        f_doppler=2*v/lamta;%真目标多普勒频率
        Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
        Aj=sqrt(Prj);
        %速度
         f_doppler1=2*v1/lamta; 
        [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
        [~,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y,tr_JD,ts,N);
        [s_echo_2,echo]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        [s_ft,echo3]=JDDeceptionJamming(D,R1,tr_JD,c,Aj,N,frame_JD,fs_JD,f_doppler1,tm,f0_JD,tau_JD,congmubiao,temp1);
        s_echo_1=s_echo_2+s_noise+s_ft;%%%%干扰信号
        figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('雷达接收信号');
        figure,plot((0:fs_JD/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_JD-fs_JD/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('雷达接收信号的频谱');
        save data/data_DeceptionJammingParameter R1 v1 congmubiao Gj Gjr Pj tf vf L
        
        case '拖引'
             %生成发射信号
        Rj=2e3;
        c=3e8;
        lamta=c/fz_JD;%波长
        tm=0:1/fs_JD:tr_JD-1/fs_JD;%一个脉冲重复周期采样序列
        N=length(tm);%一个脉冲重复周期采样点数长度
        An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%噪声强度
        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        ts=1/fs_JD; 
        Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %目标回波信号功率
        A=sqrt(Prs);%回波信号幅度
        f_doppler=2*v/lamta;%真目标多普勒频率
   
        %速度
        f_doppler=2*v/lamta;
        f_doppler1=f_doppler;
        lamta=c/fz_JD;
        f_doppler=2*v/lamta;
        [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
        ts=1/fs_JD;%采样间隔
        tm=0:1/fs_JD:tr_JD-1/fs_JD;%一个脉冲重复周期采样序列
        N=length(tm);%一个脉冲重复周期采样点数长度
      
        Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L/10)));
        Aj=sqrt(Prj);
%         figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(单位：秒)'), ylabel('y(单位：伏)'),title('发射信号');
%         figure,plot((0:fs/length(abs(fft(abs(fftshift(y))))):fs-fs/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('发射信号的频谱');
        %生成理想脉冲压缩系数
        [~,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y,tr_JD,ts,N);

        %生成回波信号
        [s_echo_2,echo]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
             [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        %生成干扰信号
        [s_ft,echo3]=JDtuoyinganrao(R,frame_JD,tf,Aj,tm,tau_JD,f0_JD,N,f_doppler1,f_doppler,vf,lamta,ts,c,fs_JD,temp1);

        %目标回波信号、假目标信号、噪声叠加在一起送入接收机
        s_echo_1=s_echo_2+s_noise+s_ft;
        figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(单位：s)'), ylabel('y(单位：伏)'),title('雷达接收信号');
        figure,plot((0:fs_JD/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_JD-fs_JD/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('雷达接收信号的频谱');
        save data/data_DeceptionJammingParameter  Gj Gjr Pj tf vf L
        
        end
     end
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
string1=get(handles.n1,'string');
string2=get(handles.n2,'string');
temp1=get(handles.n1,'Value');
temp2=get(handles.n2,'Value');
Gj=str2num(get(handles.Gj,'string'));
Gjr=str2num(get(handles.Gjr,'string')); %#ok<*NASGU,*ST2NM>
Pj=str2num(get(handles.Pj,'string'));
tf=str2num(get(handles.tf,'string'));
vf=str2num(get(handles.vf,'string'));
L=3;%综合损耗
switch string1{temp1}
    case '速度'
        switch string2{temp2}
            case '单'
                R1=2100;
                v1=270;
                congmubiao=1;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                 
            case '多'
                R1=2100;
                v1=270;
                congmubiao=2;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                 
            case '密集'
                R1=2100;
                v1=270;
                congmubiao=10;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                 
            case '拖引'
                R1=2100;
                v1=270;
                tf=10;
                vf=10;
                set(handles.R1,'string','');
                set(handles.v1,'string','');
                set(handles.congmubiao,'string','');
                set(handles.R1,'Enable','off');
                set(handles.v1,'Enable','off');
                set(handles.congmubiao,'Enable','off');
                set(handles.vf,'Enable','on');
                set(handles.tf,'Enable','on');
                set(handles.tf,'string',tf);
                set(handles.vf,'string',vf);
                
        end
        
    case '距离'
         switch string2{temp2}
            case '单'
                R1=2500;
                v1=300;
                congmubiao=1;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                
            case '多'
                R1=2500;
                v1=300;
                congmubiao=2;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                
            case '密集'
                R1=2500;
                v1=300;
                congmubiao=10;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                
            case '拖引'
                R1=2500;
                v1=300;
                tf=10;
                vf=10;
                set(handles.R1,'string','');
                set(handles.v1,'string','');
                set(handles.congmubiao,'string','');
                set(handles.R1,'Enable','off');
                set(handles.v1,'Enable','off');
                set(handles.congmubiao,'Enable','off');
                set(handles.vf,'Enable','on');
                set(handles.tf,'Enable','on');
                set(handles.tf,'string',tf);
                set(handles.vf,'string',vf);
                
         end
    case '联合'
        switch string2{temp2}
            case '单'
                R1=2500;
                v1=270;
                congmubiao=1;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                
            case '多'
                R1=2500;
                v1=270;
                congmubiao=2;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                
            case '密集'
                R1=2500;
                v1=270;
                congmubiao=10;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                
            case '拖引'
                set(handles.R1,'string','');
                set(handles.v1,'string','');
                set(handles.congmubiao,'string','');
                set(handles.R1,'Enable','off');
                set(handles.v1,'Enable','off');
                set(handles.congmubiao,'Enable','off');
                set(handles.vf,'Enable','on');
                set(handles.tf,'Enable','on');
                tf=10;
                vf=10;
                R1=2500;
                v1=270;
                set(handles.tf,'string',tf);
                set(handles.vf,'string',vf);
                
        end
end


function Gj_Callback(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to Gj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Gj as text
%        str2double(get(hObject,'String')) returns contents of Gj as a double


% --- Executes during object creation, after setting all properties.
function Gj_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Gjr_Callback(hObject, eventdata, handles)
% hObject    handle to Gjr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Gjr as text
%        str2double(get(hObject,'String')) returns contents of Gjr as a double


% --- Executes during object creation, after setting all properties.
function Gjr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gjr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pj_Callback(hObject, eventdata, handles)
% hObject    handle to Pj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pj as text
%        str2double(get(hObject,'String')) returns contents of Pj as a double


% --- Executes during object creation, after setting all properties.
function Pj_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vf_Callback(hObject, eventdata, handles)
% hObject    handle to vf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vf as text
%        str2double(get(hObject,'String')) returns contents of vf as a double


% --- Executes during object creation, after setting all properties.
function vf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function v1_Callback(hObject, eventdata, handles)
% hObject    handle to v1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of v1 as text
%        str2double(get(hObject,'String')) returns contents of v1 as a double


% --- Executes during object creation, after setting all properties.
function v1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to v1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function R1_Callback(hObject, eventdata, handles)
% hObject    handle to R1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R1 as text
%        str2double(get(hObject,'String')) returns contents of R1 as a double


% --- Executes during object creation, after setting all properties.
function R1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function congmubiao_Callback(hObject, eventdata, handles)
% hObject    handle to congmubiao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of congmubiao as text
%        str2double(get(hObject,'String')) returns contents of congmubiao as a double


% --- Executes during object creation, after setting all properties.
function congmubiao_CreateFcn(hObject, eventdata, handles)
% hObject    handle to congmubiao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in n1.
function n1_Callback(hObject, eventdata, handles)
% hObject    handle to n1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns n1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from n1


% --- Executes during object creation, after setting all properties.
function n1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in n2.
function n2_Callback(hObject, eventdata, handles)
% hObject    handle to n2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns n2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from n2
string1=get(handles.n1,'string');
string2=get(handles.n2,'string');
temp1=get(handles.n1,'Value');
temp2=get(handles.n2,'Value');
Gj=str2num(get(handles.Gj,'string'));
Gjr=str2num(get(handles.Gjr,'string')); %#ok<*NASGU,*ST2NM>
Pj=str2num(get(handles.Pj,'string'));
tf=str2num(get(handles.tf,'string'));
vf=str2num(get(handles.vf,'string'));
L=3;%综合损耗
switch string1{temp1}
    case '速度'
        switch string2{temp2}
            case '单'
                R1=2100;
                v1=270;
                congmubiao=1;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                 
            case '多'
                R1=2100;
                v1=270;
                congmubiao=2;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                 
            case '密集'
                R1=2100;
                v1=270;
                congmubiao=10;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                 
            case '拖引'
                R1=2100;
                v1=270;
                tf=10;
                vf=10;
                set(handles.R1,'string','');
                set(handles.v1,'string','');
                set(handles.congmubiao,'string','');
                set(handles.R1,'Enable','off');
                set(handles.v1,'Enable','off');
                set(handles.congmubiao,'Enable','off');
                set(handles.vf,'Enable','on');
                set(handles.tf,'Enable','on');
                set(handles.tf,'string',tf);
                set(handles.vf,'string',vf);
                
        end
        
    case '距离'
         switch string2{temp2}
            case '单'
                R1=2500;
                v1=300;
                congmubiao=1;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                
            case '多'
                R1=2500;
                v1=300;
                congmubiao=2;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                
            case '密集'
                R1=2500;
                v1=300;
                congmubiao=10;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                
            case '拖引'
                R1=2500;
                v1=300;
                tf=10;
                vf=10;
                set(handles.R1,'string','');
                set(handles.v1,'string','');
                set(handles.congmubiao,'string','');
                set(handles.R1,'Enable','off');
                set(handles.v1,'Enable','off');
                set(handles.congmubiao,'Enable','off');
                set(handles.vf,'Enable','on');
                set(handles.tf,'Enable','on');
                set(handles.tf,'string',tf);
                set(handles.vf,'string',vf);
                
         end
    case '联合'
        switch string2{temp2}
            case '单'
                R1=2500;
                v1=270;
                congmubiao=1;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                
            case '多'
                R1=2500;
                v1=270;
                congmubiao=2;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                
            case '密集'
                R1=2500;
                v1=270;
                congmubiao=10;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                
            case '拖引'
                set(handles.R1,'string','');
                set(handles.v1,'string','');
                set(handles.congmubiao,'string','');
                set(handles.R1,'Enable','off');
                set(handles.v1,'Enable','off');
                set(handles.congmubiao,'Enable','off');
                set(handles.vf,'Enable','on');
                set(handles.tf,'Enable','on');
                tf=10;
                vf=10;
                R1=2500;
                v1=270;
                set(handles.tf,'string',tf);
                set(handles.vf,'string',vf);
                
        end
end

% --- Executes during object creation, after setting all properties.
function n2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tf_Callback(hObject, eventdata, handles)
% hObject    handle to tf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tf as text
%        str2double(get(hObject,'String')) returns contents of tf as a double


% --- Executes during object creation, after setting all properties.
function tf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
