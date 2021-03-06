function varargout = A_main(varargin)
% A_MAIN MATLAB code for A_main.fig
%      A_MAIN, by itself, creates a new A_MAIN or raises the existing
%      singleton*.
%
%      H = A_MAIN returns the handle to a new A_MAIN or the handle to
%      the existing singleton*.
%
%      A_MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in A_MAIN.M with the given input arguments.
%
%      A_MAIN('Property','Value',...) creates a new A_MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before A_main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to A_main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help A_main

% Last Modified by GUIDE v2.5 06-Jul-2020 19:36:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @A_main_OpeningFcn, ...
                   'gui_OutputFcn',  @A_main_OutputFcn, ...
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


% --- Executes just before A_main is made visible.
function A_main_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to A_main (see VARARGIN)

% Choose default command line output for A_main
handles.output = hObject;
global sigma0;
sigma0=2;
global c;
c=3*10^8;
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes A_main wait for user response (see UIRESUME)
% uiwait(handles.figure1); 


% --- Outputs from this function are returned to the command line.
function varargout = A_main_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on selection change in radarSelect.
function radarSelect_Callback(hObject, ~, handles)   %#ok<*DEFNU>
global RadarS;
global str;
str=get(hObject,'string');RadarS=get(hObject,'Value');
set(handles.radar,'string',str{RadarS});

% --- Executes during object creation, after setting all properties.
function radarSelect_CreateFcn(hObject, ~, ~)
% hObject    handle to radarSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global RadarS;
RadarS=1;
% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in targetSelect.
function targetSelect_Callback(hObject, ~, handles)
% hObject    handle to targetSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 global rcsk;
str_t=get(hObject,'string');
m=get(hObject,'Value');
set(handles.target,'string',str_t{m});
switch m
    case 1
        rcsk=0;
    case 2
        rcsk=1;
    case 3
        rcsk=1;
    case 4
        rcsk=2; 
    case 5
        rcsk=2;
end
save data/rcsk rcsk;
% Hints: contents = cellstr(get(hObject,'String')) returns targetSelect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from targetSelect


% --- Executes during object creation, after setting all properties.
function targetSelect_CreateFcn(hObject, ~, ~)
% hObject    handle to targetSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global rcsk;
rcsk=0;
% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in jammingSelect.
function jammingSelect_Callback(hObject, ~, handles)
% hObject    handle to jammingSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns jammingSelect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from jammingSelect
str_j=get(hObject,'string');p=get(hObject,'Value');
set(handles.jamming,'string',str_j{p});
if p==1
        set(handles.pushbutton7,'Enable','off');
else
        set(handles.pushbutton7,'Enable','on');
end
        
% --- Executes during object creation, after setting all properties.
function jammingSelect_CreateFcn(hObject, ~, ~)
% hObject    handle to jammingSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(~, ~, handles)%????????????
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str=get(handles.radarSelect,'string');
n=get(handles.radarSelect,'Value');
switch str{n}
    case '????????????'
         h=LFMParameterSetting; %#ok<*NASGU>
    case '????????????'
         h=BFParameterSetting;   
    case '????????????'
         h=JDParameterSetting; 
end 
        
% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(~, ~, ~)%????????????
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%????????????
% str=get(handles.targetSelect,'string');
% n=get(handles.targetSelect,'Value');
        h=target0ParameterSetting;
       

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(~, ~, handles)%????????????
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str=get(handles.jammingSelect,'string');
n=get(handles.jammingSelect,'Value');
switch str{n}
     
    case '??????????'
         h=DeceptionJammingParameter; 
    case '??????????'
         h=oppressSet;
% h=oppressivejammingParameterSeting;

    case '????'
         h=other;
end 

% --- Executes during object creation, after setting all properties.
function pushbutton5_CreateFcn(~, ~, ~)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on button press in generateSignal.
function generateSignal_Callback(~, ~, handles)% %????????
% hObject    handle to generateSignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% h=LFMParameterSetting,
%  BI=str2double(get(findall(h ,'tag','BI'),'string'));   %
str=get(handles.radarSelect,'string');
n=get(handles.radarSelect,'Value');
switch str{n}
    case '????????????'
       load data/data_LFMParameter %#ok<*LOAD>
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%????????
        % global R;
        % global v;
        % global sigma0;
        % global rcsk;  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%????????????
        ts=1/fs_LFM;
        k=B1_LFM/tau_LFM;                                 %????????????????????
        tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;  
        N=length(tm);
        [y,~]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
        figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(????????)'), ylabel('y(????????)'),title('LFM????????');
        figure,plot((0:fs_LFM/length(abs(fft(abs(fftshift(y))))):fs_LFM-fs_LFM/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('????f(??????Hz)'), ylabel('y(????????/Hz)'),title('LFM??????????????');

    case '????????????'
        global code ;
        load data/data_BFParameter
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%????????
        % global R;
        % global v;
        % global sigma0;
        % global rcsk;  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%????????????
        tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
        N=length(tm_B);%????????????????????????????
        %code=[1,1,1,-1,-1,-1,1];
        number1=length(code);
        ts=1/fs_B;
        [y,~,~]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
        figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(????????)'), ylabel('y(????????)'),title('????????');
        figure,plot((0:fs_B/length(abs(fft(abs(fftshift(y))))):fs_B-fs_B/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('????f(??????Hz)'), ylabel('y(????????/Hz)'),title('??????????????');
        % [y,A,Aj,An,ts,f_doppler,f_doppler1,O,N,fr]=shengchengLFMxinhao(B,Gt,Gr,lamta,F,Srmin,Pt,L,Pfa,R,v,R1,v1,sigma,Rmax,Rmin,Rf0,Gj,Gjr,rj,Kj,Te,Pj,tau,tr,f0,f1,fs,time_jam,frame,num_jilei,num_tongdao,num_cankao,num_baohu,number,deltf,deltt,Rj);
        % y_fft_result=(fft(y));
   
    case '????????????'
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%????????
        load data/data_JDParameter
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%????????
        % global R;
        % global v;
        % global sigma0;
        % global rcsk;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%????????????
        ts=1/fs_JD;
        tm=0:1/fs_JD:tr_JD-1/fs_JD;  
        N=length(tm);
        [y,~]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
        figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(????????)'), ylabel('y(????????)'),title('????????????????');
        figure,plot((0:fs_JD/length(abs(fft(abs(fftshift(y))))):fs_JD-fs_JD/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('????f(??????Hz)'), ylabel('y(????????/Hz)'),title('??????????????????????');

end 


% --- Executes on button press in echo_signal.
function echo_signal_Callback(~, ~, handles)%????????
% hObject    handle to echo_signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


str=get(handles.radarSelect,'string');
n=get(handles.radarSelect,'Value');
global c;
global code 
c=3e8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%????????


load data/data_target0Parameter
load data/antannaParameter
global sigma;
global rcsk;
sigma = rcs(rcsk,sigma0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if get(handles.jammingSelect,'value')==1
switch str{n}
    case '????????????'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LFM????????
       load data/data_LFMParameter
       %%%%%%%%%%%%LFM????????????    
        lamta=c/fz_LFM;%????
        ts=1/fs_LFM;
        k=B1_LFM/tau_LFM;                                 %????????????????????
        tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;  
        N=length(tm);
        f_doppler=2*v/lamta;%????????????????
        Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
        A=sqrt(Prs);%????????????
        [s_echo_2,echo]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
         
        [~,RTAz,RTEl] = xyz2radar(radar.x,radar.y,radar.z,...
                    target.x,target.y,target.z);                        %????????????
%         [FTR,FTAz,FTEl] = xyz2radar(radar.x1,radar.y1,radar.z1,...
%                     jammer.x3,jammer.y3,jammer.z3);                        %????????????
        Gt = 10^(radar.Gt/10);                                                           %????????????????
%         Gj = jammer.Gain;                                                          %??????????????????
        Gr = gain(radar,RTAz,RTEl);                                                %??????????????????
%         Grj = gain(radar,FTAz,FTEl)/Gj;/Gt                                               %????????????????????
        echo1.sum=s_echo_2*Gr.Gainsum;
        echo1.Az=s_echo_2*Gr.GainAz;
        echo1.Ez=s_echo_2*Gr.GainEl;
        
        figure,plot(0:ts:(N-1)*ts,real(echo1.sum(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('LFM????????');
        figure,plot((0:fs_LFM/length(abs(fft(abs(fftshift(echo(1,:)))))):fs_LFM-fs_LFM/length(abs(fft(abs(fftshift(echo(1,:))))))),abs((fft(echo(1,:))))),xlabel('????f(??????Hz)'), ylabel('y(????????/Hz)'),title('LFM??????????????');     
    
    case '????????????'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BK????????
        load data/data_BFParameter
        tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
        N=length(tm_B);%????????????????????????????
        number1=length(code);
        ts=1/fs_B;
        lamta=c/fz_B;
        Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
        A=sqrt(Prs);%????????????
        f_doppler=2*v/lamta;%????????????????
        [~,y1,~]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
        [s_echo_2,echo]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
        
        [~,RTAz,RTEl] = xyz2radar(radar.x,radar.y,radar.z,...
                    target.x,target.y,target.z);                        %????????????
%         [FTR,FTAz,FTEl] = xyz2radar(radar.x1,radar.y1,radar.z1,...
%                     jammer.x3,jammer.y3,jammer.z3);                        %????????????
        Gt = 10^(radar.Gt/10);                                                           %????????????????
%         Gj = jammer.Gain;                                                          %??????????????????
        Gr = gain(radar,RTAz,RTEl);                                                %??????????????????
%         Grj = gain(radar,FTAz,FTEl)/Gj;/Gt                                               %????????????????????
        echo1.sum=s_echo_2*Gr.Gainsum;
        echo1.Az=s_echo_2*Gr.GainAz;
        echo1.Ez=s_echo_2*Gr.GainEl;
        
        
        figure,plot(0:ts:(N-1)*ts,real( echo1.sum(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????????');
        figure,plot((0:fs_B/length(abs(fft(abs(fftshift(echo(1,:)))))):fs_B-fs_B/length(abs(fft(abs(fftshift(echo(1,:))))))),abs((fft(echo(1,:))))),xlabel('????f(??????Hz)'), ylabel('y(????????/Hz)'),title('??????????????????????');
 
    case '????????????'
       %%%%%%%%%%%%%%%%%%%%%%%%%
     load data/data_JDParameter
        %%%%%%%%%%%????????????????????
        ts=1/fs_JD;
        lamta=c/fz_JD;
        f_doppler=2*v/lamta;
        Prs=(Pt_JD*(10^((Gt_JD*0.1)))*(10^((Gr_JD*0.1)))*(lamta^2)*sigma)/((4*pi)^3*R^4*10^(L_JD*0.1)); %????????????????
        A=sqrt(Prs);%????????????
        tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
        N=length(tm);%????????????????????????????
        [~,~]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
        [s_echo_2,echo]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
        
        [~,RTAz,RTEl] = xyz2radar(radar.x,radar.y,radar.z,...
                    target.x,target.y,target.z);                        %????????????
%         [FTR,FTAz,FTEl] = xyz2radar(radar.x1,radar.y1,radar.z1,...
%                     jammer.x3,jammer.y3,jammer.z3);                        %????????????
        Gt = 10^(radar.Gt/10);                                                           %????????????????
%         Gj = jammer.Gain;                                                          %??????????????????
        Gr = gain(radar,RTAz,RTEl);                                                %??????????????????
%         Grj = gain(radar,FTAz,FTEl)/Gj;/Gt                                               %????????????????????
        echo1.sum=s_echo_2*Gr.Gainsum;
        echo1.Az=s_echo_2*Gr.GainAz;
        echo1.Ez=s_echo_2*Gr.GainEl;
        
        figure,plot(0:ts:(N-1)*ts,real(echo1.sum(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????????');
        figure,plot((0:fs_JD/length(abs(fft(abs(fftshift(echo(1,:)))))):fs_JD-fs_JD/length(abs(fft(abs(fftshift(echo(1,:))))))),abs((fft(echo(1,:))))),xlabel('????f(??????Hz)'), ylabel('y(????????/Hz)'),title('??????????????????????');    
end

% --- Executes on button press in matchfilter.
function matchfilter_Callback(~, ~, ~)%????????
% hObject    handle to matchfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%????????
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%??????????????????
global RadarS;
global c;
c=3e8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%????????
  global code ;
load data/data_target0Parameter
load data/antannaParameter;
  global sigma;
  global rcsk;
  sigma = rcs(rcsk,sigma0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if RadarS==1
        %%%%%%%%%%%%%%????????????
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LFM????????
        load data/data_LFMParameter
        lamta=c/fz_LFM;%????
        tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
        N=length(tm);%????????????????????????????
        An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
        ts=1/fs_LFM;
        k=B1_LFM/tau_LFM;   
        Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
        A=sqrt(Prs);%????????????
        f_doppler=2*v/lamta;%????????????????%????????????????????   
        [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
         [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
        [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k);
        
        [~,RTAz,RTEl] = xyz2radar(radar.x,radar.y,radar.z,...
                    target.x,target.y,target.z);                        %????????????
%         [FTR,FTAz,FTEl] = xyz2radar(radar.x1,radar.y1,radar.z1,...
%                     jammer.x3,jammer.y3,jammer.z3);                        %????????????
        Gt = 10^(radar.Gt/10);                                                           %????????????????
%         Gj = jammer.Gain;                                                          %??????????????????
        Gr = gain(radar,RTAz,RTEl);                                                %??????????????????
%         Grj = gain(radar,FTAz,FTEl)/Gj;/Gt                                               %????????????????????
        echo1.sum=s_echo_2*Gr.Gainsum;
        echo1.Az=s_echo_2*Gr.GainAz;
        echo1.Ez=s_echo_2*Gr.GainEl;
        
        
        
        s_echo_2=echo1.sum;
        [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
        s_echo_1=s_echo_2+s_noise;
        [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
        [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
        [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
        [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
        figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('LFM????????');

elseif RadarS==2  %%%%%%%%%%%%%????????????
        
%         [M,match_filter_fft]=maiyaxishu(f0,fs,y,tr,ts,N);
%         [s_echo_2,echo]=BKhuiboxinhao(y1,R,tr,ts,A,N,frame,fs,f_doppler,tau); 
%         [s_noise]=zaosheng(frame,N,An,B,fs);     
%         s_echo_1=s_echo_2+s_noise;
%         [s_echo_1]=gaofang(f0,B,fs,s_echo_1);
%         [s_echo_1,f0]=hunpin(s_echo_1,N,frame,f1,fs,f0);
%         [s_echo_mf]=jianbo(s_echo_1,N,frame,f0,fs);
%         [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame,match_filter_fft,tau,D,ts);
%         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BK????????
         load data/data_BFParameter
        
        tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
        N=length(tm_B);%????????????????????????????
        number1=length(code);
        ts=1/fs_B;
        lamta=c/fz_B;
        An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
        Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
        A=sqrt(Prs);%????????????
        f_doppler=2*v/lamta;%????????????????
        [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
        [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
        
         [~,RTAz,RTEl] = xyz2radar(radar.x,radar.y,radar.z,...
                    target.x,target.y,target.z);                        %????????????
%         [FTR,FTAz,FTEl] = xyz2radar(radar.x1,radar.y1,radar.z1,...
%                     jammer.x3,jammer.y3,jammer.z3);                        %????????????
        Gt = 10^(radar.Gt/10);                                                           %????????????????
%         Gj = jammer.Gain;                                                          %??????????????????
        Gr = gain(radar,RTAz,RTEl);                                                %??????????????????
%         Grj = gain(radar,FTAz,FTEl)/Gj;/Gt                                               %????????????????????
        echo1.sum=s_echo_2*Gr.Gainsum;
        echo1.Az=s_echo_2*Gr.GainAz;
        echo1.Ez=s_echo_2*Gr.GainEl;
        
        s_echo_2=echo1.sum;
        
        [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
        [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
        s_echo_1=s_echo_2+s_noise;
        [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
        [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
        [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
        [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
        figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????????????');

  
elseif RadarS==3
        %%%%%%%%%%%%%%%%%%%%%%%????????
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%JD????????
        load data/data_JDParameter
        lamta=c/fz_JD;%????
        tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
        N=length(tm);%????????????????????????????
        An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        ts=1/fs_JD; 
        Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
        A=sqrt(Prs);%????????????
        f_doppler=2*v/lamta;%????????????????%???????????????????? 
        [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
        [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
        [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD);
        
         [~,RTAz,RTEl] = xyz2radar(radar.x,radar.y,radar.z,...
                    target.x,target.y,target.z);                        %????????????
%         [FTR,FTAz,FTEl] = xyz2radar(radar.x1,radar.y1,radar.z1,...
%                     jammer.x3,jammer.y3,jammer.z3);                        %????????????
        Gt = 10^(radar.Gt/10);                                                           %????????????????
%         Gj = jammer.Gain;                                                          %??????????????????
        Gr = gain(radar,RTAz,RTEl);                                                %??????????????????
%         Grj = gain(radar,FTAz,FTEl)/Gj;/Gt                                               %????????????????????
        echo1.sum=s_echo_2*Gr.Gainsum;
        echo1.Az=s_echo_2*Gr.GainAz;
        echo1.Ez=s_echo_2*Gr.GainEl;
        
        
        
        s_echo_2=echo1.sum;
        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        s_echo_1=s_echo_2+s_noise;
        [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
        [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
        [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
        [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
        figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????????????');

        
end

% --- Executes on button press in mtd.
function mtd_Callback(~, ~, ~)%MTD
% hObject    handle to mtd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%??????????????????
global RadarS;
global c;
c=3e8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%????????


load data/data_target0Parameter
load data/antannaParameter;
  global sigma;
  global rcsk;
  sigma = rcs(rcsk,sigma0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if RadarS==1
        %%%%%%%%%%%%%%????????????
        load data/data_LFMParameter
        fr=1/tr_LFM;%??????????
        lamta=c/fz_LFM;%????
        tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
        N=length(tm);%????????????????????????????
        An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
        [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
        ts=1/fs_LFM;
        k=B1_LFM/tau_LFM;   
        Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
        A=sqrt(Prs);%????????????
        f_doppler=2*v/lamta;%????????????????%???????????????????? 
        [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
        [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
        [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
        
        
         [~,RTAz,RTEl] = xyz2radar(radar.x,radar.y,radar.z,...
                    target.x,target.y,target.z);                        %????????????
%         [FTR,FTAz,FTEl] = xyz2radar(radar.x1,radar.y1,radar.z1,...
%                     jammer.x3,jammer.y3,jammer.z3);                        %????????????
        Gt = 10^(radar.Gt/10);                                                           %????????????????
%         Gj = jammer.Gain;                                                          %??????????????????
        Gr = gain(radar,RTAz,RTEl);                                                %??????????????????
%         Grj = gain(radar,FTAz,FTEl)/Gj;/Gt                                               %????????????????????
        echo1.sum=s_echo_2*Gr.Gainsum;
        echo1.Az=s_echo_2*Gr.GainAz;
        echo1.Ez=s_echo_2*Gr.GainEl;
        
        
        
        s_echo_2=echo1.sum;
        [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
        s_echo_1=s_echo_2+s_noise;
        [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
        [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
        [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
%         [pc_result,~,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts); 
        [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
        [s_mtd]=mtd(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
        figure,mesh(1:fr/num_tongdao_LFM:fr,0:ts*(tau_LFM/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_LFM/D/ts)-ts*(tau_LFM/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('LFM????????MTD????');

elseif RadarS==2  %%%%%%%%%%%%%????????????
         global code ;
         load data/data_BFParameter
        tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
        N=length(tm_B);%????????????????????????????
        number1=length(code);
        ts=1/fs_B;
        lamta=c/fz_B;
        fr=1/tr_B;
        An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
        Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
        A=sqrt(Prs);%????????????
        f_doppler=2*v/lamta;%????????????????
        [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
        [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
        
         [~,RTAz,RTEl] = xyz2radar(radar.x,radar.y,radar.z,...
                    target.x,target.y,target.z);                        %????????????
%         [FTR,FTAz,FTEl] = xyz2radar(radar.x1,radar.y1,radar.z1,...
%                     jammer.x3,jammer.y3,jammer.z3);                        %????????????
        Gt = 10^(radar.Gt/10);                                                           %????????????????
%         Gj = jammer.Gain;                                                          %??????????????????
        Gr = gain(radar,RTAz,RTEl);                                                %??????????????????
%         Grj = gain(radar,FTAz,FTEl)/Gj;/Gt                                               %????????????????????
        echo1.sum=s_echo_2*Gr.Gainsum;
        echo1.Az=s_echo_2*Gr.GainAz;
        echo1.Ez=s_echo_2*Gr.GainEl;
        
        
        
        s_echo_2=echo1.sum;
        [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
        [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
        s_echo_1=s_echo_2+s_noise;
        [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
        [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
        [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
        [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
        [s_mtd]=mtd(pc_result1.',M1,num_jilei_B,num_tongdao_B);
        figure,mesh(1:fr/num_tongdao_B:fr,0:ts*(tau_B/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_B/D/ts)-ts*(tau_B/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('????????MTD????');
  
        
elseif RadarS==3
        %%%%%%%%%%%%%%%%%%%%%%%????????
         load data/data_JDParameter
         fr=1/tr_JD;
        lamta=c/fz_JD;%????
        tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
        N=length(tm);%????????????????????????????
        An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
        ts=1/fs_JD; 
        Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
        A=sqrt(Prs);%????????????
        f_doppler=2*v/lamta;%???????????????? 
        [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
        [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
        [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
        
        
         [~,RTAz,RTEl] = xyz2radar(radar.x,radar.y,radar.z,...
                    target.x,target.y,target.z);                        %????????????
%         [FTR,FTAz,FTEl] = xyz2radar(radar.x1,radar.y1,radar.z1,...
%                     jammer.x3,jammer.y3,jammer.z3);                        %????????????
        Gt = 10^(radar.Gt/10);                                                           %????????????????
%         Gj = jammer.Gain;                                                          %??????????????????
        Gr = gain(radar,RTAz,RTEl);                                                %??????????????????
%         Grj = gain(radar,FTAz,FTEl)/Gj;/Gt                                               %????????????????????
        echo1.sum=s_echo_2*Gr.Gainsum;
        echo1.Az=s_echo_2*Gr.GainAz;
        echo1.Ez=s_echo_2*Gr.GainEl;
        
        
        
        s_echo_2=echo1.sum;
        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        s_echo_1=s_echo_2+s_noise;
        [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
        [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
        [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
        [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
        [s_mtd]=mtd(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
        figure,mesh(1:fr/num_tongdao_JD:fr,0:ts*(tau_JD/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_JD/D/ts)-ts*(tau_JD/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('????????MTD????');

        
        %%%%%%%%%%%%%%%%%%%%%%%%%
end

% --- Executes on button press in hengxujing.
function hengxujing_Callback(~, ~, ~)%??????
% hObject    handle to hengxujing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%??????????????????
global RadarS;
global c;
c=3e8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%????????


load data/data_target0Parameter
load data/antannaParameter;
  global sigma;
  global rcsk;
  sigma = rcs(rcsk,sigma0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if RadarS==1
        %%%%%%%%%%%%%%????????????
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LFM????????
 load data/data_LFMParameter
        fr=1/tr_LFM;%????????????
      
        lamta=c/fz_LFM;%????
        tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
        N=length(tm);%????????????????????????????
        An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
        [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
        ts=1/fs_LFM;
        k=B1_LFM/tau_LFM;   
        Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
        A=sqrt(Prs);%????????????
        f_doppler=2*v/lamta;%????????????????%???????????????????? 
        [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
        [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
        [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
        
        
         [~,RTAz,RTEl] = xyz2radar(radar.x,radar.y,radar.z,...
                    target.x,target.y,target.z);                        %????????????
%         [FTR,FTAz,FTEl] = xyz2radar(radar.x1,radar.y1,radar.z1,...
%                     jammer.x3,jammer.y3,jammer.z3);                        %????????????
        Gt = 10^(radar.Gt/10);                                                           %????????????????
%         Gj = jammer.Gain;                                                          %??????????????????
        Gr = gain(radar,RTAz,RTEl);                                                %??????????????????
%         Grj = gain(radar,FTAz,FTEl)/Gj;/Gt                                               %????????????????????
        echo1.sum=s_echo_2*Gr.Gainsum;
        echo1.Az=s_echo_2*Gr.GainAz;
        echo1.Ez=s_echo_2*Gr.GainEl;
        
        
        
        s_echo_2=echo1.sum;
        [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
        s_echo_1=s_echo_2+s_noise;
        [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
        [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
        [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
        [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts); 
        [s_mtd]=mtd(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
        hengxujing(M1,Pfa_LFM,s_mtd,num_cankao_LFM,num_tongdao_LFM,num_baohu_LFM,ts,c,tau_LFM,D);
        %%%%%%%%%%%%%%%%%%%%%%%%%
elseif RadarS==2  %%%%%%%%%%%%%????????????
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BF????????
         global code 
 load data/data_BFParameter
         tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
        N=length(tm_B);%????????????????????????????
        number1=length(code);
        ts=1/fs_B;
        lamta=c/fz_B;
        fr=1/tr_B;
        An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
        Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
        A=sqrt(Prs);%????????????
        f_doppler=2*v/lamta;%????????????????
        [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
        [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
        
        
         [~,RTAz,RTEl] = xyz2radar(radar.x,radar.y,radar.z,...
                    target.x,target.y,target.z);                        %????????????
%         [FTR,FTAz,FTEl] = xyz2radar(radar.x1,radar.y1,radar.z1,...
%                     jammer.x3,jammer.y3,jammer.z3);                        %????????????
        Gt = 10^(radar.Gt/10);                                                           %????????????????
%         Gj = jammer.Gain;                                                          %??????????????????
        Gr = gain(radar,RTAz,RTEl);                                                %??????????????????
%         Grj = gain(radar,FTAz,FTEl)/Gj;/Gt                                               %????????????????????
        echo1.sum=s_echo_2*Gr.Gainsum;
        echo1.Az=s_echo_2*Gr.GainAz;
        echo1.Ez=s_echo_2*Gr.GainEl;
        
        
        
        s_echo_2=echo1.sum;
        [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
        [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
        s_echo_1=s_echo_2+s_noise;
        [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
        [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
        [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
        [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
        [s_mtd]=mtd(pc_result1.',M1,num_jilei_B,num_tongdao_B);
        hengxujing(M1,Pfa_B,s_mtd,num_cankao_B,num_tongdao_B,num_baohu_B,ts,c,tau_B,D);
elseif RadarS==3
        %%%%%%%%%%%%%%%%%%%%%%%????????
         load data/data_JDParameter

        fr=1/tr_JD;
        lamta=c/fz_JD;%????
        tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
        N=length(tm);%????????????????????????????
        An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        ts=1/fs_JD; 
        Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
        A=sqrt(Prs);%????????????
        f_doppler=2*v/lamta;%????????????????%???????????????????? 
        [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
        [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
        [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
        
        
         [~,RTAz,RTEl] = xyz2radar(radar.x,radar.y,radar.z,...
                    target.x,target.y,target.z);                        %????????????
%         [FTR,FTAz,FTEl] = xyz2radar(radar.x1,radar.y1,radar.z1,...
%                     jammer.x3,jammer.y3,jammer.z3);                        %????????????
        Gt = 10^(radar.Gt/10);                                                           %????????????????
%         Gj = jammer.Gain;                                                          %??????????????????
        Gr = gain(radar,RTAz,RTEl);                                                %??????????????????
%         Grj = gain(radar,FTAz,FTEl)/Gj;/Gt                                               %????????????????????
        echo1.sum=s_echo_2*Gr.Gainsum;
        echo1.Az=s_echo_2*Gr.GainAz;
        echo1.Ez=s_echo_2*Gr.GainEl;
        
        
        
        s_echo_2=echo1.sum;
        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        s_echo_1=s_echo_2+s_noise;
        [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
        [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
        [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
        [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
        [s_mtd]=mtd(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
        hengxujing(M1,Pfa_JD,s_mtd,num_cankao_JD,num_tongdao_JD,num_baohu_JD,ts,c,tau_JD,D);
        
end


% --- Executes on button press in zongtifangzhen.
function zongtifangzhen_Callback(~, ~, ~)%????????
% hObject    handle to zongtifangzhen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%??????????????????
global RadarS;
global c;
c=3e8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%????????
global sigma0;
global rcsk;

load data/data_target0Parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if RadarS==1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LFM????????
         load data/data_LFMParameter
        LFMmain(fz_LFM,B_LFM,B1_LFM,Gt_LFM,Gr_LFM,F_LFM,Pt_LFM,L_LFM,Te_LFM,tau_LFM,tr_LFM,f0_LFM,f1_LFM,fs_LFM,frame_LFM,num_jilei_LFM,num_tongdao_LFM,num_cankao_LFM,num_baohu_LFM,Pfa_LFM,R,v,rcsk,sigma0) ;
elseif RadarS==2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BF????????
        global code  %#ok<*TLEV>
        load data/data_BFParameter
        number1=length(code);
        BKmain(fz_B,B_B,Gt_B,Gr_B,F_B,Pt_B,L_B,Pfa_B,R,v,rcsk,sigma0,Te_B,tau_B,tr_B,f0_B,f1_B,fs_B,frame_B,num_jilei_B,num_tongdao_B,num_cankao_B,num_baohu_B,flag,number1,code)  ;
elseif RadarS==3     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%JD????????
 load data/data_JDParameter
        JDmain(fz_JD,B_JD,Gt_JD,Gr_JD,F_JD,Pt_JD,L_JD,Pfa_JD,R,v,rcsk,sigma0,Te_JD,tau_JD,tr_JD,f0_JD,f1_JD,fs_JD,frame_JD,num_jilei_JD,num_tongdao_JD,num_cankao_JD,num_baohu_JD) ;
end

function zongtifangzhen_CreateFcn(~, ~, ~)
% hObject    handle to zongtifangzhen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function pushbutton7_CreateFcn(hObject, ~, ~)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Enable','off');


   % --- Executes on button press in pushbutton19.
function pushbutton19_Callback(~, ~, handles)%??????????????
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global sigma;
global rcsk;
global sigma0;
sigma = rcs(rcsk,sigma0);
global temp1;
global string2;
global temp2;
global c;
global code;
global strOPSet;
global nOPSet;
global strO;
global nO;
str0=get(handles.jammingSelect,'string');
n0=get(handles.jammingSelect,'Value');

if n0==3
str=get(handles.radarSelect,'string');
n=get(handles.radarSelect,'Value');
switch str{n}
    case '????????????'
       load data/data_LFMParameter
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%????????
        load data/data_DeceptionJammingParameter
        load data/data_target0Parameter
        switch string2{temp2}
            case {'??' ,'??','????' }
                    %%%%%%%%%%%%%%%%????????
                     load data/data_LFMParameter
                    %%%%%%%%%%%%%%%%%%%%????????
                     load data/data_DeceptionJammingParameter
                     load data/data_target0Parameter
                    % global rcsk;  
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%????????????
                    fr=1/tr_LFM;%??????????
                    lamta=c/fz_LFM;%#ok<*NODEF> %????
                    tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                    N=length(tm);%????????????????????????????
                    An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                    ts=1/fs_LFM;
                    k=B1_LFM/tau_LFM;   
                    Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                    A=sqrt(Prs);%????????????
                    f_doppler=2*v/lamta;%????????????????%????????????????????
                    Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
                    Aj=sqrt(Prj);
                    %????
                    f_doppler1=2*v1/lamta;
                    [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                    [~,~]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                    [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                    [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                    [s_ft,echo3]=LFMDeceptionJamming(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao,temp1);
                    s_echo_1=s_echo_2+s_noise+s_ft;%%%%????????
                    figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                    figure,plot((0:fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_LFM-fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('????f(??????Hz)'), ylabel('y(????????/Hz)'),title('??????????????????');
               
            case '????'
                load data/data_LFMParameter
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%????????
                load data/data_DeceptionJammingParameter
                 load data/data_target0Parameter
                c=3e8;
                lamta=c/fz_LFM;
                Rj=2e3;
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L_LFM/10)));
                Aj=sqrt(Prj);
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                f_doppler=2*v/lamta;
                f_doppler1=f_doppler;
                k=B1_LFM/tau_LFM;%????????????????????
                ts=1/fs_LFM;%????????
                %????????????
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L/10))); %????????????????
                A=sqrt(Prs);
                Te=1143;%????
                F=4;%????????
                 tf=10;
                An=10*log10((1.382e-23)*Te*B_LFM*10^(F/10));%????????
                [y,~]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k); 
                [~,~]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                %????????????
                [s_echo_2,~]=LFMtuoyinhuiboxinhao(y,R,tr_LFM,ts,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                %figure,plot(0:ts:(N-1)*ts,real(s_echo_2(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????');
                %figure,plot((0:fs/length(abs(fft(abs(fftshift(echo(1,:)))))):fs-fs/length(abs(fft(abs(fftshift(echo(1,:))))))),abs((fft(echo(1,:))))),xlabel('????f(??????Hz)'), ylabel('y(????????)'),title('??????????????');

                %????????????
                [s_ft,echo3]=LFMtuoyinganrao(R,frame_LFM,tf,Aj,tm,tau_LFM,f0_LFM,B1_LFM,N,f_doppler1,f_doppler,vf,lamta,k,ts,c,fs_LFM,temp1);

                %????????
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                %??????????????????????????????????????????????????
                s_echo_1=s_echo_2+s_noise+s_ft;
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                figure,plot((0:fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_LFM-fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('????f(??????Hz)'), ylabel('y(????????/Hz)'),title('??????????????????');
                
        end
    case '????????????'
        c=3e8;
        load data/data_DeceptionJammingParameter;
        load data/data_BFParameter;
        load data/data_target0Parameter;
        switch string2{temp2}
            case {'??' ,'??','????' }
                    tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                    N=length(tm_B);%????????????????????????????
                    number1=length(code);
                    ts=1/fs_B;
                    lamta=c/fz_B;
                    fr=1/tr_B;
                    An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                    Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                    A=sqrt(Prs);%????????????
                    f_doppler=2*v/lamta;%????????????????
                    Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
                    Aj=sqrt(Prj);
                    %????
                    f_doppler1=2*v1/lamta; 
                    [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                    [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                    [~,~]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                    [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                    [s_ft,echo3]=BKDeceptionJamming(D,y,R1,tr_B,ts,c,Aj,N,frame_B,fs_B,f_doppler1,tm_B,f0_B,tau_B,congmubiao,y1,temp1);
                    s_echo_1=s_echo_2+s_noise+s_ft;%%%%????????
                    figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                    figure,plot((0:fs_B/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_B-fs_B/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('????f(??????Hz)'), ylabel('y(????????/Hz)'),title('??????????????????');
                    
            case '????'
                Rj=2e3;
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                tm=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm);%????????????????????????????
                f_doppler=2*v/lamta;%????????????????
                f_doppler1= f_doppler;%????????????????
                Te=1143;%????
                An=10*log10((1.382e-23)*Te*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L/10)));
                Aj=sqrt(Prj);
                [y,y1,~]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
    %             figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(????????)'), ylabel('y(????????)'),title('????????');
    %             figure,plot((0:fs_B/length(abs(fft(abs(fftshift(y))))):fs_B-fs_B/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('????f(??????Hz)'), ylabel('y(????????)'),title('??????????????');
    %             %%%%%%%%%%%1.1.2????????????????????
                [~,~]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                 %%%%%%%%%%%%%1.3????????????%%%%%%%%%%%%
                [s_echo_2,~]=BKtuoyinhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tm,tau_B); 
                %????????
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                %????????
                [s_ft,echo3]=BKtuoyinganrao(R,frame_B,tf,Aj,tau_B,N,f_doppler1,f_doppler,vf,lamta,ts,c,fs_B,y1,tr_B,temp1);
                %??????????????????????????????????????????????????
                s_echo_1=s_echo_2+s_noise+s_ft;
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                figure,plot((0:fs_B/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_B-fs_B/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('????f(??????Hz)'), ylabel('y(????????/Hz)'),title('??????????????????');

        end
    case '????????????'
        c=3e8;
        load data/data_JDParameter;
        load data/data_DeceptionJammingParameter;
        load data/data_target0Parameter;
        switch string2{temp2}
            case {'??' ,'??','????' }
                    fr=1/tr_JD;
                    lamta=c/fz_JD;%????
                    tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                    N=length(tm);%????????????????????????????
                    An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                    [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                    ts=1/fs_JD; 
                    Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                    A=sqrt(Prs);%????????????
                    f_doppler=2*v/lamta;%????????????????
                    Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
                    Aj=sqrt(Prj);
                    %????
                     f_doppler1=2*v1/lamta; 
                    [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                    [~,~]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                    [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                    [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                    [s_ft,echo3]=JDDeceptionJamming(D,R1,tr_JD,c,Aj,N,frame_JD,fs_JD,f_doppler1,tm,f0_JD,tau_JD,congmubiao,temp1);
                    s_echo_1=s_echo_2+s_noise+s_ft;%%%%????????
                    figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                    figure,plot((0:fs_JD/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_JD-fs_JD/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('????f(??????Hz)'), ylabel('y(????????/Hz)'),title('??????????????????');
   
            case '????'
                Rj=2e3;
                c=3e8;
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????

                %????
                f_doppler=2*v/lamta;
                f_doppler1=f_doppler;
                lamta=c/fz_JD;
                f_doppler=2*v/lamta;
                [y,~]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                ts=1/fs_JD;%????????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????

                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L/10)));
                Aj=sqrt(Prj);
                %figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(????????)'), ylabel('y(????????)'),title('????????');
                %figure,plot((0:fs/length(abs(fft(abs(fftshift(y))))):fs-fs/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('????f(??????Hz)'), ylabel('y(????????)'),title('??????????????');
                %????????????????????
                [~,~]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);

                %????????????
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                     [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                %????????????
                [s_ft,echo3]=JDtuoyinganrao(R,frame_JD,tf,Aj,tm,tau_JD,f0_JD,N,f_doppler1,f_doppler,vf,lamta,ts,c,fs_JD,temp1);

                %??????????????????????????????????????????????????
                s_echo_1=s_echo_2+s_noise+s_ft;
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                figure,plot((0:fs_JD/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_JD-fs_JD/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('????f(??????Hz)'), ylabel('y(????????/Hz)'),title('??????????????????');

        end
end
elseif n0==2

      switch strOPSet{nOPSet}
          
      case '????????????'
     
        load data/data_sp
              
        Kfm=4e6;tau=1e-6; c=3e8;
        % global str;
        % global RadarS;
        str=get(handles.radarSelect,'string');
        n=get(handles.radarSelect,'Value');
        [noise_sp]=shepinzaosheng(fs_sp,Bj_sp,fj_sp,frame_sp,Prj_sp,Tr_sp);
        switch str{n}   
            case '????????????'
        % Tr=40e-6;tau=1e-6;frame=64;fs=160e6;f0=10e6;B1=4e6;A=5;R=2e3;tm=0:1/fs:Tr-1/fs;k=B1/tau;N=length(tm);Pn=(B1/(2.5*Kfm))^2;Bn=B1/2;
                load data/data_LFMParameter;
                load data/data_target0Parameter
                lamta=c/fz_LFM;%???? 
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;                                 %????????????????????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;  
                N=length(tm);
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));
                f_doppler=2*v/lamta;%????????????????
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [vRadarSig]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                s_echo_1=vRadarSig+s_noise+noise_sp;%%%%????????
                t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                
          case '????????????'
               load data/data_BFParameter;
               load data/data_target0Parameter 
               tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                number1=length(code);
                N=length(tm_B);%????????????????????????????
                ts=1/fs_B;
                lamta=c/fz_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                 Pn=(B_B/(2.5*Kfm))^2;Bn=B_B/2;
                 [~,y1,~]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [vRadarSig]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                s_echo_1=vRadarSig+s_noise+noise_sp;%%%%????????
                t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                
           case '????????????'
               load data/data_JDParameter;
               load data/data_target0Parameter
                ts=1/fs_JD;
                lamta=c/fz_JD;
                f_doppler=2*v/lamta;
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));
                Prs=(Pt_JD*(10^((Gt_JD*0.1)))*(10^((Gr_JD*0.1)))*(lamta^2)*sigma)/((4*pi)^3*R^4*10^(L_JD*0.1)); %????????????????
                A=sqrt(Prs);%????????????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
        %       [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                Pn=(B_JD/(2.5*Kfm))^2;Bn=B_JD/2;
                [vRadarSig]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD);
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                s_echo_1=vRadarSig+s_noise+noise_sp;%%%%????????
                t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD;
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                
        end
     case '????????????'
        load data/data_tx   
        load data/data_sp     
        Kfm=4e6;tau=1e-6;  c=3e8;
        % global str;
        % global RadarS;
        str=get(handles.radarSelect,'string');
        n=get(handles.radarSelect,'Value');
        [noise_tx]=zaoshengtiaoxiang(fs_tx,Bj_tx,fj_tx,Prj_tx,Tr_tx,frame_tx);
         switch str{n}   
            case '????????????'
        % Tr=40e-6;tau=1e-6;frame=64;fs=160e6;f0=10e6;B1=4e6;A=5;R=2e3;tm=0:1/fs:Tr-1/fs;k=B1/tau;N=length(tm);Pn=(B1/(2.5*Kfm))^2;Bn=B1/2;
                load data/data_LFMParameter;
                load data/data_target0Parameter
                lamta=c/fz_LFM;%???? 
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;                                 %????????????????????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;  
                N=length(tm);
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));
                f_doppler=2*v/lamta;%????????????????
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [vRadarSig]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                s_echo_1=vRadarSig+s_noise+noise_tx;%%%%????????
                 t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                 s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                
            case '????????????'
                load data/data_BFParameter;
                load data/data_target0Parameter 
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                number1=length(code);
                N=length(tm_B);%????????????????????????????
                ts=1/fs_B;
                lamta=c/fz_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                 Pn=(B_B/(2.5*Kfm))^2;Bn=B_B/2;
                 [~,y1,~]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [vRadarSig]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                s_echo_1=vRadarSig+s_noise+noise_tx;%%%%????????
                 t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                 s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                
           case '????????????'
               load data/data_JDParameter;
               load data/data_target0Parameter
                ts=1/fs_JD;
                lamta=c/fz_JD;
                f_doppler=2*v/lamta;
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));
                Prs=(Pt_JD*(10^((Gt_JD*0.1)))*(10^((Gr_JD*0.1)))*(lamta^2)*sigma)/((4*pi)^3*R^4*10^(L_JD*0.1)); %????????????????
                A=sqrt(Prs);%????????????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
        %       [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                Pn=(B_JD/(2.5*Kfm))^2;Bn=B_JD/2;
                [vRadarSig]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD);
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                s_echo_1=vRadarSig+s_noise+noise_tx;%%%%????????
                 t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD;
                 s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
               
         end
    case '????????????'
     
        load data/data_sp
        Kfm=4e6;tau=1e-6; c=3e8;
        % global str;
        % global RadarS;
        str=get(handles.radarSelect,'string');
        n=get(handles.radarSelect,'Value');
       load data/data_tf
       [noise_tf] =zaoshengtiaofu(fs_tf,Bj_tf,fj_tf,Prj_tf,Tr_tf,frame_tf);
        switch str{n}   
            case '????????????'
        % Tr=40e-6;tau=1e-6;frame=64;fs=160e6;f0=10e6;B1=4e6;A=5;R=2e3;tm=0:1/fs:Tr-1/fs;k=B1/tau;N=length(tm);Pn=(B1/(2.5*Kfm))^2;Bn=B1/2;
                load data/data_LFMParameter;
                load data/data_target0Parameter
                lamta=c/fz_LFM;%???? 
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;                                 %????????????????????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;  
                N=length(tm);
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));
                f_doppler=2*v/lamta;%????????????????
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [vRadarSig]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                s_echo_1=vRadarSig+s_noise+noise_tf;%%%%????????
                 t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                 s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
               
            case '????????????'
               load data/data_BFParameter;
               load data/data_target0Parameter 
               tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                number1=length(code);
                N=length(tm_B);%????????????????????????????
                ts=1/fs_B;
                lamta=c/fz_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                 Pn=(B_B/(2.5*Kfm))^2;Bn=B_B/2;
                 [~,y1,~]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [vRadarSig]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                s_echo_1=vRadarSig+s_noise+noise_tf;%%%%????????
                 t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                 s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                
   
           case '????????????'
               load data/data_JDParameter;
               load data/data_target0Parameter
                ts=1/fs_JD;
                lamta=c/fz_JD;
                f_doppler=2*v/lamta;
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));
                Prs=(Pt_JD*(10^((Gt_JD*0.1)))*(10^((Gr_JD*0.1)))*(lamta^2)*sigma)/((4*pi)^3*R^4*10^(L_JD*0.1)); %????????????????
                A=sqrt(Prs);%????????????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
        %       [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                Pn=(B_JD/(2.5*Kfm))^2;Bn=B_JD/2;
                [vRadarSig]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD);
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                s_echo_1=vRadarSig+s_noise+noise_tf;%%%%????????
                 t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD; 
                 s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                

        end

    case '????????????'
    
        Kfm=4e6;tau=1e-6; c=3e8;
        % global str;
        % global RadarS;
        str=get(handles.radarSelect,'string');
        n=get(handles.radarSelect,'Value');
        load data/data_tp
        Pn=(Bj_tp/2/(2.5*Kfm))^2;
        [noise_tp]=zaoshengtiaopin(fs_tp,Kfm,Prj_tp,Pn,Bn,fj_tp,frame_tp,Tr_tp);
        switch str{n}   
            case '????????????'
        % Tr=40e-6;tau=1e-6;frame=64;fs=160e6;f0=10e6;B1=4e6;A=5;R=2e3;tm=0:1/fs:Tr-1/fs;k=B1/tau;N=length(tm);Pn=(B1/(2.5*Kfm))^2;Bn=B1/2;
                load data/data_LFMParameter;
                load data/data_target0Parameter
                lamta=c/fz_LFM;%???? 
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;                                 %????????????????????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;  
                N=length(tm);
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));
                f_doppler=2*v/lamta;%????????????????
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [vRadarSig]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                s_echo_1=vRadarSig+s_noise+noise_tp;%%%%????????
                 t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM;
                 s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                
            case '????????????'
               load data/data_BFParameter;
               load data/data_target0Parameter 
               tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                number1=length(code);
                N=length(tm_B);%????????????????????????????
                ts=1/fs_B;
                lamta=c/fz_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                 Pn=(B_B/(2.5*Kfm))^2;Bn=B_B/2;
                 [~,y1,~]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [vRadarSig]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                s_echo_1=vRadarSig+s_noise+noise_tp;%%%%????????
                 t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                 s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                   
           case '????????????'
               load data/data_JDParameter;
               load data/data_target0Parameter
                ts=1/fs_JD;
                lamta=c/fz_JD;
                f_doppler=2*v/lamta;
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));
                Prs=(Pt_JD*(10^((Gt_JD*0.1)))*(10^((Gr_JD*0.1)))*(lamta^2)*sigma)/((4*pi)^3*R^4*10^(L_JD*0.1)); %????????????????
                A=sqrt(Prs);%????????????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
        %       [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                Pn=(B_JD/(2.5*Kfm))^2;Bn=B_JD/2;
                [vRadarSig]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD);
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                s_echo_1=vRadarSig+s_noise+noise_tp;%%%%????????
                 t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD; 
                 s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                
        end

    case  '??????????'
      load data/data_shuzhuangpu
      fj=[1e7,4e7,7e7];
      [noise_szp,t_noise] = jam_combspectrum(fs_shuzhuangpu,Bj_shuzhuangpu,Ns_shuzhuangpu,fj,frame_shuzhuangpu,Prj_shuzhuangpu,Tr_shuzhuangpu);
%       view_jam_combspectrum( sig_noise,t_noise,fs_shuzhuangpu );
            c=3e8;
            str=get(handles.radarSelect,'string');
            n=get(handles.radarSelect,'Value');
    %         load data/data_tp
    %         Pn=(Bj_tp/2/(2.5*Kfm))^2;
    %         [noise_tp]=zaoshengtiaopin(fs_tp,Kfm,Prj_tp,Pn,Bn,fj_tp,frame_tp,Tr_tp)
        switch str{n}   
            case '????????????'
        % Tr=40e-6;tau=1e-6;frame=64;fs=160e6;f0=10e6;B1=4e6;A=5;R=2e3;tm=0:1/fs:Tr-1/fs;k=B1/tau;N=length(tm);Pn=(B1/(2.5*Kfm))^2;Bn=B1/2;
                load data/data_LFMParameter;
                load data/data_target0Parameter
                lamta=c/fz_LFM;%???? 
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;                                 %????????????????????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;  
                N=length(tm);
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));
                f_doppler=2*v/lamta;%????????????????
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [vRadarSig]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                sig_noise=vRadarSig+s_noise+noise_szp;%%%%????????
                t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                sig_noise=sig_noise.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                figure,plot(0:ts:(N-1)*ts,real(sig_noise(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');

            case '????????????'
               load data/data_BFParameter;
               load data/data_target0Parameter 
               tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                number1=length(code);
                N=length(tm_B);%????????????????????????????
                ts=1/fs_B;
                lamta=c/fz_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
%                 Pn=(B_B/(2.5*Kfm))^2;Bn=B_B/2;
                [~,y1,~]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [vRadarSig]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                sig_noise=vRadarSig+s_noise+noise_szp;%%%%????????
                t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                sig_noise=sig_noise.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                figure,plot(0:ts:(N-1)*ts,real(sig_noise(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
   
           case '????????????'
               load data/data_JDParameter;
               load data/data_target0Parameter
                ts=1/fs_JD;
                lamta=c/fz_JD;
                f_doppler=2*v/lamta;
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));
                Prs=(Pt_JD*(10^((Gt_JD*0.1)))*(10^((Gr_JD*0.1)))*(lamta^2)*sigma)/((4*pi)^3*R^4*10^(L_JD*0.1)); %????????????????
                A=sqrt(Prs);%????????????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
        %       [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
%                 Pn=(B_JD/(2.5*Kfm))^2;Bn=B_JD/2;
                [vRadarSig]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD);
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                sig_noise=vRadarSig+s_noise+noise_szp;%%%%????????
                t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD; 
                sig_noise=sig_noise.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                figure,plot(0:ts:(N-1)*ts,real(sig_noise(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                
        end
    case   '????????????'
        
        load data/data_smart
        Kfm=4e6;tau=1e-6; c=3e8;
        str=get(handles.radarSelect,'string');
        n=get(handles.radarSelect,'Value');
        
        switch str{n}   
            case '????????????'
     
                load data/data_LFMParameter;
                load data/data_target0Parameter;
                lamta=c/fz_LFM;%???? 
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;                                 %????????????????????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;  
                N=length(tm);
                f_doppler=2*v/lamta;%????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [vRadarSig]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                Pn=(B1_LFM/(2.5*Kfm))^2;Bn=B1_LFM/2;
                % [vRadarSig]=LFMhuiboxinhao(R,c,A,N,frame,fs,0,tm,f0,B1,tau,k);
                [vSmartNoiseSig]=jam_smartnoise( vRadarSig,Pn,Prj_smart,Bn,Kfm,fs_smart );
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);                
                sig_noise=vRadarSig+s_noise+vSmartNoiseSig;%%%%????????
                t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                sig_noise=sig_noise.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                figure,plot(0:ts:(N-1)*ts,real(sig_noise(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
            case '????????????'
                load data/data_BFParameter;
                load data/data_target0Parameter; 
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                number1=length(code);
                N=length(tm_B);%????????????????????????????
                ts=1/fs_B;
                lamta=c/fz_B;
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));
                 Pn=(B_B/(2.5*Kfm))^2;Bn=B_B/2;
                 [~,y1,~]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [vRadarSig]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B);
                [vSmartNoiseSig]=jam_smartnoise( vRadarSig,Pn,Prj_smart,Bn,Kfm,fs_smart );
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                sig_noise=vRadarSig+s_noise+vSmartNoiseSig;%%%%????????
                t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                sig_noise=sig_noise.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                figure,plot(0:ts:(N-1)*ts,real(sig_noise(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');  
            case '????????????'
                load data/data_JDParameter;
                load data/data_target0Parameter;
                ts=1/fs_JD;
                lamta=c/fz_JD;
                f_doppler=2*v/lamta;
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));
                Prs=(Pt_JD*(10^((Gt_JD*0.1)))*(10^((Gr_JD*0.1)))*(lamta^2)*sigma)/((4*pi)^3*R^4*10^(L_JD*0.1)); %????????????????
                A=sqrt(Prs);%????????????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                Pn=(B_JD/(2.5*Kfm))^2;Bn=B_JD/2;
                [vRadarSig]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD);
                [vSmartNoiseSig]=jam_smartnoise( vRadarSig,Pn,Prj_smart,Bn,Kfm,fs_smart );
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                sig_noise=vRadarSig+s_noise+vSmartNoiseSig;%%%%????????
                t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD; 
                sig_noise=sig_noise.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                figure,plot(0:ts:(N-1)*ts,real(sig_noise(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
        end
    case '????????'
               load data/data_saopin;
                % T_fr=2*Tr_saopin
                 [ noise_saopin,t_noise ] = jam_sweepfrequency( fs_saopin,Bj_saopin,fj_saopin,frame_saopin,Prj_saopin,Tr_saopin,T_fr_saopin,Time_begin_saopin,K_sweep_saopin );
                 Kfm=4e6;tau=1e-6; c=3e8;
            str=get(handles.radarSelect,'string');
            n=get(handles.radarSelect,'Value');
        switch str{n}   
            case '????????????'
        % Tr=40e-6;tau=1e-6;frame=64;fs=160e6;f0=10e6;B1=4e6;A=5;R=2e3;tm=0:1/fs:Tr-1/fs;k=B1/tau;N=length(tm);Pn=(B1/(2.5*Kfm))^2;Bn=B1/2;
                load data/data_LFMParameter;
                load data/data_target0Parameter;
                lamta=c/fz_LFM;%???? 
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;                                 %????????????????????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;  
                N=length(tm);
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));
                f_doppler=2*v/lamta;%????????????????
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [vRadarSig]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                sig_noise=vRadarSig+s_noise+noise_saopin;%%%%????????
                t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                sig_noise=sig_noise.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                figure,plot(0:ts:(N-1)*ts,real(sig_noise(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');

            case '????????????'
               load data/data_BFParameter;
               load data/data_target0Parameter; 
               tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                number1=length(code);
                N=length(tm_B);%????????????????????????????
                ts=1/fs_B;
                lamta=c/fz_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                Pn=(B_B/(2.5*Kfm))^2;Bn=B_B/2;
                [~,y1,~]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [vRadarSig]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                sig_noise=vRadarSig+s_noise+noise_saopin;%%%%????????
                t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                sig_noise=sig_noise.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                figure,plot(0:ts:(N-1)*ts,real(sig_noise(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');

           case '????????????'
               load data/data_JDParameter;
               load data/data_target0Parameter;
                ts=1/fs_JD;
                lamta=c/fz_JD;
                f_doppler=2*v/lamta;
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));
                Prs=(Pt_JD*(10^((Gt_JD*0.1)))*(10^((Gr_JD*0.1)))*(lamta^2)*sigma)/((4*pi)^3*R^4*10^(L_JD*0.1)); %????????????????
                A=sqrt(Prs);%????????????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
        %       [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                Pn=(B_JD/(2.5*Kfm))^2;Bn=B_JD/2;
                [vRadarSig]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD);
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                sig_noise=vRadarSig+s_noise+noise_saopin;%%%%????????
                t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD; 
                sig_noise=sig_noise.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                figure,plot(0:ts:(N-1)*ts,real(sig_noise(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');

        end

      end
       
    %end
elseif  n0==4

        switch strO{nO}
            case '??????????????'
                
            load data/data_dopplerblink;
            Kfm=4e6;tau=1e-6; c=3e8;
            c=3*10^8;
            R0=rand(100);
            str=get(handles.radarSelect,'string');
            n=get(handles.radarSelect,'Value');
            
            switch str{n}
                case '????????????'
                    load data/data_LFMParameter;
                    load data/data_target0Parameter;
                    lamta=c/fz_LFM;%???? 
                    ts=1/fs_LFM;
                    k=B1_LFM/tau_LFM;                                 %????????????????????
                    tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;  
                    N=length(tm);
                    An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));
                    f_doppler=2*v/lamta;%????????????????
                    Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                    A=sqrt(Prs);%????????????
                    [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k);
                    [ sig_jam,t_jam ] = jam_dopplerblink( fd_dopplerblink,Td_dopplerblink,R0,s_echo_2,fs_dopplerblink,Pj_dopplerblink,flagT_dopplerblink );
                    [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                    s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                    view_jam_dopplerblink( s_echo_1,sig_jam,t_jam,fs_dopplerblink );
                    figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                    
                case '????????????'
                    load data/data_BFParameter;
                    load data/data_target0Parameter ;
                    tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                    number1=length(code);
                    N=length(tm_B);%????????????????????????????
                    ts=1/fs_B;
                    lamta=c/fz_B;
                    An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));
                    Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                    A=sqrt(Prs);%????????????
                    f_doppler=2*v/lamta;%????????????????
                    Pn=(B_B/(2.5*Kfm))^2;Bn=B_B/2;
                    [~,y1,~]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                    [vRadarSig]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B);
                    [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                    [ sig_jam,t_jam ] = jam_dopplerblink( fd_dopplerblink,Td_dopplerblink,R0,vRadarSig,fs_dopplerblink,Pj_dopplerblink,flagT_dopplerblink );
                    sig_noise=vRadarSig+s_noise+sig_jam;%%%%????????
                    view_jam_dopplerblink( sig_noise,sig_jam,t_jam,fs_dopplerblink );
                    figure,plot(0:ts:(N-1)*ts,real(sig_noise(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                    
                case '????????????'
                    load data/data_JDParameter;
                    load data/data_target0Parameter;
                    ts=1/fs_JD;
                    lamta=c/fz_JD;
                    f_doppler=2*v/lamta;
                    An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));
                    Prs=(Pt_JD*(10^((Gt_JD*0.1)))*(10^((Gr_JD*0.1)))*(lamta^2)*sigma)/((4*pi)^3*R^4*10^(L_JD*0.1)); %????????????????
                    A=sqrt(Prs);%????????????
                    tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                    N=length(tm);%????????????????????????????
                    Pn=(B_JD/(2.5*Kfm))^2;Bn=B_JD/2;
                    [vRadarSig]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD);
                    [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                    [ sig_jam,t_jam ] = jam_dopplerblink( fd_dopplerblink,Td_dopplerblink,R0,vRadarSig,fs_dopplerblink,Pj_dopplerblink,flagT_dopplerblink );
                    sig_noise=vRadarSig+s_noise+sig_jam;%%%%????????
                    view_jam_dopplerblink( sig_noise,sig_jam,t_jam,fs_dopplerblink );
                    figure,plot(0:ts:(N-1)*ts,real(sig_noise(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                    
            end

            case '????????'
                load data/data_botiao
                Kfm=4e6;tau=1e-6; c=3e8;
                c=3*10^8;
                str=get(handles.radarSelect,'string');
                n=get(handles.radarSelect,'Value');
                switch str{n}
                    case '????????????'
                        load data/data_LFMParameter;
                        load data/data_target0Parameter;
                        lamta=c/fz_LFM;%????
                        tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                        N=length(tm);%????????????????????????????
                        An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                        R0=rand(100);
                        ts=1/fs_LFM;
                        k=B1_LFM/tau_LFM;   
                        Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                        A=sqrt(Prs);%????????????
                        [y,~]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                        [~,~]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                        f_doppler=2*v/lamta;%????????????????%????????????????????  
                        [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                        [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);

                         %????????????
                        px=1e3;py=1e3;pz=1e3;%????????
                        vx=10;vy=10;vz=0;%????????
                        ax=0;ay=0;az=0;%??????????
                        phi=pi/180;
                        [ TargetStatus ] = paraset_targetstatus( px,py,pz,vx,vy,vz,ax,ay,az,phi );
                        [ PassivePara ] = paraset_passivejaming( tf_botiao,sf_botiao,vl_botiao,vf_botiao,ts_botiao,bt_botiao,al_botiao,sref_botiao,smax_botiao );
                          vx=1;vy=1;
                        [ WindV ] = paraset_windvelocity( vx,vy );
                        CurrentT=2.2;
                        [ sig_jam,t_jam ] = jam_passive( s_echo_2,fs_LFM,f0_LFM,CurrentT,TargetStatus,WindV,PassivePara );
                        s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                        view_jam_passive( s_echo_1,t_jam,fs_LFM);
                        figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                        
                    case'????????????'
                        load data/data_BFParameter;
                        load data/data_target0Parameter ;
                        tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                        number1=length(code);
                        N=length(tm_B);%????????????????????????????
                        ts=1/fs_B;
                        lamta=c/fz_B;
                        An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));
                        Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                        A=sqrt(Prs);%????????????
                        f_doppler=2*v/lamta;%????????????????
                        Pn=(B_B/(2.5*Kfm))^2;Bn=B_B/2;
                        [~,y1,~]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                        [vRadarSig]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B);
                        [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                        
                         %????????????
                        px=1e3;py=1e3;pz=1e3;%????????
                        vx=10;vy=10;vz=0;%????????
                        ax=0;ay=0;az=0;%??????????
                        phi=pi/180;
                        [ TargetStatus ] = paraset_targetstatus( px,py,pz,vx,vy,vz,ax,ay,az,phi );
                        [ PassivePara ] = paraset_passivejaming( tf_botiao,sf_botiao,vl_botiao,vf_botiao,ts_botiao,bt_botiao,al_botiao,sref_botiao,smax_botiao );
                          vx=1;vy=1;
                        [ WindV ] = paraset_windvelocity( vx,vy );
                        CurrentT=2.2;
                        [ sig_jam,t_jam ] = jam_passive( vRadarSig,fs_B,f0_B,CurrentT,TargetStatus,WindV,PassivePara );
                        s_echo_1=vRadarSig+s_noise+sig_jam;%%%%????????
                        view_jam_passive( s_echo_1,t_jam,fs_B ); 
                        figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                        
                    case '????????????'
                        load data/data_JDParameter;
                        load data/data_target0Parameter;
                        ts=1/fs_JD;
                        lamta=c/fz_JD;
                        f_doppler=2*v/lamta;
                        An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));
                        Prs=(Pt_JD*(10^((Gt_JD*0.1)))*(10^((Gr_JD*0.1)))*(lamta^2)*sigma)/((4*pi)^3*R^4*10^(L_JD*0.1)); %????????????????
                        A=sqrt(Prs);%????????????
                        tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                        N=length(tm);%????????????????????????????
                        Pn=(B_JD/(2.5*Kfm))^2;Bn=B_JD/2;
                        [vRadarSig]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD);
                        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);

                         %????????????
                        px=1e3;py=1e3;pz=1e3;%????????
                        vx=10;vy=10;vz=0;%????????
                        ax=0;ay=0;az=0;%??????????
                        phi=pi/180;
                        [ TargetStatus ] = paraset_targetstatus( px,py,pz,vx,vy,vz,ax,ay,az,phi );
                        [ PassivePara ] = paraset_passivejaming( tf_botiao,sf_botiao,vl_botiao,vf_botiao,ts_botiao,bt_botiao,al_botiao,sref_botiao,smax_botiao );
                          vx=1;vy=1;
                        [ WindV ] = paraset_windvelocity( vx,vy );
                        CurrentT=2.2;
                        [ sig_jam,t_jam ] = jam_passive( vRadarSig,fs_JD,f0_JD,CurrentT,TargetStatus,WindV,PassivePara );
                        s_echo_1=vRadarSig+s_noise+sig_jam;%%%%????????
                        view_jam_passive( s_echo_1,t_jam,fs_JD );
                        figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                        
                end

            case 'AGC????'
                load data/data_AGC
                Kfm=4e6;tau=1e-6; c=3e8;
                c=3*10^8;
                str=get(handles.radarSelect,'string');
                n=get(handles.radarSelect,'Value');
                switch str{n}
                    case '????????????'
                        load data/data_LFMParameter;
                        load data/data_target0Parameter;
                        lamta=c/fz_LFM;%????
                        tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                        N=length(tm);%????????????????????????????
                        An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                        R0=rand(100);
                        ts=1/fs_LFM;
                        k=B1_LFM/tau_LFM;   
                        Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                        A=sqrt(Prs);%????????????
                        [y,~]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                        [~,~]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                        f_doppler=2*v/lamta;%????????????????%????????????????????  
                        [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                        [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                        [ sig_jam,t_jam ] = jam_AGC( CurrentT_AGC,Pj_AGC,Period_AGC,D_AGC,radio_AGC,s_echo_2,fs_AGC );
                        s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                        view_jam_AGC( s_echo_1,t_jam,fs_LFM );
                        figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                        
                    case'????????????'
                        load data/data_BFParameter;
                        load data/data_target0Parameter ;
                        tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                        number1=length(code);
                        N=length(tm_B);%????????????????????????????
                        ts=1/fs_B;
                        lamta=c/fz_B;
                        An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));
                        Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                        A=sqrt(Prs);%????????????
                        f_doppler=2*v/lamta;%????????????????
                        Pn=(B_B/(2.5*Kfm))^2;Bn=B_B/2;
                        [~,y1,~]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                        [vRadarSig]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B);
                        [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                        [ sig_jam,t_jam ] = jam_AGC( CurrentT_AGC,Pj_AGC,Period_AGC,D_AGC,radio_AGC,vRadarSig,fs_AGC );
                        s_echo_1=vRadarSig+s_noise+sig_jam;%%%%????????
                        view_jam_AGC( s_echo_1,t_jam,fs_B ); 
                        figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                        
                    case '????????????'
                        load data/data_JDParameter;
                        load data/data_target0Parameter;
                        ts=1/fs_JD;
                        lamta=c/fz_JD;
                        f_doppler=2*v/lamta;
                        An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));
                        Prs=(Pt_JD*(10^((Gt_JD*0.1)))*(10^((Gr_JD*0.1)))*(lamta^2)*sigma)/((4*pi)^3*R^4*10^(L_JD*0.1)); %????????????????
                        A=sqrt(Prs);%????????????
                        tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                        N=length(tm);%????????????????????????????
                        Pn=(B_JD/(2.5*Kfm))^2;Bn=B_JD/2;
                        [vRadarSig]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD);
                        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                        [ sig_jam,t_jam ] = jam_AGC( CurrentT_AGC,Pj_AGC,Period_AGC,D_AGC,radio_AGC,vRadarSig,fs_AGC );
                        s_echo_1=vRadarSig+s_noise+sig_jam;%%%%????????
                        view_jam_AGC( s_echo_1,t_jam,fs_JD );
                        figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                end
                
        end
      
end

% --- Executes on button press in Jmatchfilter.
function Jmatchfilter_Callback(~, ~, handles)%??????????????
% hObject    handle to Jmatchfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sigma;
global rcsk;
global sigma0;
sigma = rcs(rcsk,sigma0);
global temp1;
global string2;
global temp2;
global str;
global strOPSet;
global nOPSet;
global strO;
global nO;
str0=get(handles.jammingSelect,'string');
n0=get(handles.jammingSelect,'Value');
if n0==3
        str=get(handles.radarSelect,'string');
        n=get(handles.radarSelect,'Value');
        global c;
        c=3e8;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%????????
        global code ;
        
        load data/data_DeceptionJammingParameter 
        load data/data_target0Parameter
    switch str{n}
        case '????????????'
        load data/data_LFMParameter;
        switch string2{temp2}
            case {'??' ,'??','????' }
                    lamta=c/fz_LFM;%????
                    tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                    N=length(tm);%????????????????????????????
                    An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                    ts=1/fs_LFM;
                    k=B1_LFM/tau_LFM;   
                    Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                    A=sqrt(Prs);%????????????
                    f_doppler=2*v/lamta;%????????????????%????????????????????  
                    f_doppler1=2*v1/lamta;
                    Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
                    Aj=sqrt(Prj);
                    [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                    [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                    [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                    [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                    [s_ft,~]=LFMDeceptionJamming(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao,temp1);
                    s_echo_1=s_echo_2+s_noise+s_ft;%%%%????????
                    [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                    [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                    [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                    [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                    figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('LFM????????');
                     
            case '????'
                f1=10e6; %????????
                c=3e8;
                lamta=c/fz_LFM;
                Rj=2e3;
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L_LFM/10)));
                Aj=sqrt(Prj);
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                f_doppler=2*v/lamta;
                f_doppler1=f_doppler;
                k=B1_LFM/tau_LFM;%????????????????????
                ts=1/fs_LFM;%????????
                %????????????
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L/10))); %????????????????
                A=sqrt(Prs);
                Te=1143;%????
                F=4;%????????
                 tf=10;
                An=10*log10((1.382e-23)*Te*B_LFM*10^(F/10));%????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k); 
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                %????????????
               [s_echo_2,~]=LFMtuoyinhuiboxinhao(y,R,tr_LFM,ts,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
              %figure,plot(0:ts:(N-1)*ts,real(s_echo_2(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????');
                %figure,plot((0:fs/length(abs(fft(abs(fftshift(echo(1,:)))))):fs-fs/length(abs(fft(abs(fftshift(echo(1,:))))))),abs((fft(echo(1,:))))),xlabel('????f(??????Hz)'), ylabel('y(????????)'),title('??????????????');

                %????????????
                [s_ft,~]=LFMtuoyinganrao(R,frame_LFM,tf,Aj,tm,tau_LFM,f0_LFM,B1_LFM,N,f_doppler1,f_doppler,vf,lamta,k,ts,c,fs_LFM,temp1);

                %????????
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                %??????????????????????????????????????????????????
                s_echo_1=s_echo_2+s_noise+s_ft;
                 %????
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1,fs_LFM,f0_LFM);
                %????????????
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                %????????????????
                [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                figure,plot(0:ts*c/2:(M-1)*ts*c/2,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????');
                            
        end
        case '????????????'
            load data/data_BFParameter;

            switch string2{temp2}
            case {'??' ,'??','????' }
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
                Aj=sqrt(Prj);
                %????
                f_doppler1=2*v1/lamta; 
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                [s_ft,~]=BKDeceptionJamming(D,y,R1,tr_B,ts,c,Aj,N,frame_B,fs_B,f_doppler1,tm_B,f0_B,tau_B,congmubiao,y1,temp1);
                s_echo_1=s_echo_2+s_noise+s_ft;%%%%???????? 
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
                [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
                figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????????????');
                
            case '????'
                 Rj=2e3;
                 f1=10e6;
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                tm=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm);%????????????????????????????
                f_doppler=2*v/lamta;%????????????????
                f_doppler1= f_doppler;%????????????????
                Te=1143;%????
                An=10*log10((1.382e-23)*Te*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L/10)));
                Aj=sqrt(Prj);
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                %figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(????????)'), ylabel('y(????????)'),title('????????');
                %figure,plot((0:fs_B/length(abs(fft(abs(fftshift(y))))):fs_B-fs_B/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('????f(??????Hz)'), ylabel('y(????????)'),title('??????????????');
                %%%%%%%%%%%%1.1.2????????????????????
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                 %%%%%%%%%%%%%1.3????????????%%%%%%%%%%%%
                [s_echo_2,~]=BKtuoyinhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tm,tau_B); 
                %????????
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                %????????
                [s_ft,~]=BKtuoyinganrao(R,frame_B,tf,Aj,tau_B,N,f_doppler1,f_doppler,vf,lamta,ts,c,fs_B,y1,tr_B,temp1);
                %??????????????????????????????????????????????????
                s_echo_1=s_echo_2+s_noise+s_ft;
                %????
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1,fs_B,f0_B);
                %????????????
                [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
                %????????????????
                [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
                figure,plot(0:ts*c/2:(M-1)*ts*c/2,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????');

            end
       case '????????????'
           load data/data_JDParameter;
        f1=10e6; %????????
         c=3e8;
         Rj=2e3;
        switch string2{temp2}
        case {'??' ,'??','????' }
            fr=1/tr_JD;
            lamta=c/fz_JD;%????
            tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
            N=length(tm);%????????????????????????????
            An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
            [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
            ts=1/fs_JD; 
            Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
            A=sqrt(Prs);%????????????
            f_doppler=2*v/lamta;%????????????????
            Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
            Aj=sqrt(Prj);
            %????
             f_doppler1=2*v1/lamta; 
            [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
            [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
            [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
            [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
            [s_ft,~]=JDDeceptionJamming(D,R1,tr_JD,c,Aj,N,frame_JD,fs_JD,f_doppler1,tm,f0_JD,tau_JD,congmubiao,temp1);
            s_echo_1=s_echo_2+s_noise+s_ft;%%%%????????
            [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
            [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
            [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
            [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
            figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????????????');
           
        case '????'
            lamta=c/fz_JD;%????
            tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
            N=length(tm);%????????????????????????????
            An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
            [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
            ts=1/fs_JD; 
            Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
            A=sqrt(Prs);%????????????
            f_doppler=2*v/lamta;%????????????????

            %????
            f_doppler=2*v/lamta;
            f_doppler1=f_doppler;
            lamta=c/fz_JD;
            f_doppler=2*v/lamta;
            [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
            ts=1/fs_JD;%????????
            tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
            N=length(tm);%????????????????????????????

            Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L/10)));
            Aj=sqrt(Prj);
            %         figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(????????)'), ylabel('y(????????)'),title('????????');
            %         figure,plot((0:fs/length(abs(fft(abs(fftshift(y))))):fs-fs/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('????f(??????Hz)'), ylabel('y(????????)'),title('??????????????');
            %????????????????????
            [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);

            %????????????
            [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                 [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
            %????????????
            [s_ft,~]=JDtuoyinganrao(R,frame_JD,tf,Aj,tm,tau_JD,f0_JD,N,f_doppler1,f_doppler,vf,lamta,ts,c,fs_JD,temp1);

            %??????????????????????????????????????????????????
            s_echo_1=s_echo_2+s_noise+s_ft;
            [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1,fs_JD,f0_JD);
            %????????????
            [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
            %????????????????
            [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
            figure,plot(0:ts*c/2:(M-1)*ts*c/2,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????');

        end
            
    end
elseif n0==2

    str=get(handles.radarSelect,'string');
    n=get(handles.radarSelect,'Value');
    switch  str{n}
        case '????????????'
                 c=3e8;
                 load data/data_LFMParameter;
                 load data/data_target0Parameter
%                 global strOPSet;
%                 global nOPSet;
                switch strOPSet{nOPSet}
            case '????????????'
                load data/data_sp
                lamta=c/fz_LFM;%????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                f_doppler=2*v/lamta;%????????????????%????????????????????  
                [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                [noise_sp] = shepinzaosheng(fs_sp,Bj_sp,fj_sp,frame_sp,Prj_sp,Tr_sp);
                s_echo_1=s_echo_2+s_noise+noise_sp;%%%%????????
                t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('LFM????????');
       
           case '????????????'
                load data/data_tx
                lamta=c/fz_LFM;%????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                f_doppler=2*v/lamta;%????????????????%????????????????????  
                [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                [noise_tx] =zaoshengtiaoxiang(fs_tx,Bj_tx,fj_tx,Prj_tx,Tr_tx,frame_tx);
                s_echo_1=s_echo_2+s_noise+noise_tx;%%%%????????
                t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('LFM????????');
          

          case '????????????'
                            load data/data_tf
%                 load data/data_tx
                lamta=c/fz_LFM;%????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                f_doppler=2*v/lamta;%????????????????%????????????????????  
                [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
%                 [s_ft,echo3]=LFMganraoxinhao(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao);
                [noise_tf] = zaoshengtiaofu(fs_tf,Bj_tf,fj_tf,Prj_tf,Tr_tf,frame_tf);
                s_echo_1=s_echo_2+s_noise+noise_tf;%%%%????????
                t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('LFM????????');
          
            case '????????????'
                 load data/data_tp
                 Pn=(Bj_tp/2/(2.5*Kfm))^2;
                lamta=c/fz_LFM;%????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                f_doppler=2*v/lamta;%????????????????%????????????????????  
                [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
%                 [s_ft,echo3]=LFMganraoxinhao(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao);
                [noise_tp] =zaoshengtiaopin(fs_tp,Kfm,Prj_tp,Pn,Bn,fj_tp,frame_tp,Tr_tp);
                s_echo_1=s_echo_2+s_noise+noise_tp;%%%%????????
                t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('LFM????????');
                
             case '????????????'
                    load data/data_smart
                    Kfm=4e6;tau=1e-6; c=3e8;
                   Pn=(B_LFM/(2.5*Kfm))^2;Bn=B_LFM/2;
                lamta=c/fz_LFM;%????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                f_doppler=2*v/lamta;%????????????????%????????????????????  
                [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                [vSmartNoiseSig]=jam_smartnoise( s_echo_2,Pn,Prj_smart,Bn,Kfm,fs_smart );
                s_echo_1=s_echo_2+s_noise+vSmartNoiseSig;%%%%????????
                t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('LFM????????');
            case '??????????'
                   load data/data_shuzhuangpu
                    fj=[0.4e6,0.8e6,1.2e6];

                lamta=c/fz_LFM;%????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                f_doppler=2*v/lamta;%????????????????%????????????????????  
                [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
%                 [s_ft,echo3]=LFMganraoxinhao(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao);
               [noise_szp,~] = jam_combspectrum(fs_shuzhuangpu,Bj_shuzhuangpu,Ns_shuzhuangpu,fj,frame_shuzhuangpu,Prj_shuzhuangpu,Tr_shuzhuangpu);
                s_echo_1=s_echo_2+s_noise+noise_szp;%%%%????????
                t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('LFM????????');
           case '????????'
                load data/data_saopin
                T_fr=2*Tr_saopin;
                lamta=c/fz_LFM;%????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                f_doppler=2*v/lamta;%????????????????%????????????????????  
                [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
%                 [s_ft,echo3]=LFMganraoxinhao(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao);
                [ sig_noise,~ ] = jam_sweepfrequency( fs_saopin,Bj_saopin,fj_saopin,frame_saopin,Prj_saopin,Tr_saopin,T_fr,Time_begin_saopin,K_sweep_saopin );
                s_echo_1=s_echo_2+s_noise+sig_noise;%%%%????????
                t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('LFM????????');
                
                end 

            
    case '????????????'
                load data/data_BFParameter;
                load data/data_target0Parameter
                c=3e8;
                switch strOPSet{nOPSet}
            case '????????????'
                load data/data_sp
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
%                 Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
%                 Aj=sqrt(Prj);
                %????
                
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
%                 [s_ft,echo3]=BKganraoxinhao(D,y,R1,tr_B,ts,c,Aj,N,frame_B,fs_B,f_doppler1,tm_B,f0_B,tau_B,congmubiao,y1);
                [noise_sp] = shepinzaosheng(fs_sp,Bj_sp,fj_sp,frame_sp,Prj_sp,Tr_sp);
                s_echo_1=s_echo_2+s_noise+noise_sp;%%%%????????
                t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
               [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
               [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
               figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????????????');

            case '????????????'
                load data/data_tp
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
%                 Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
%                 Aj=sqrt(Prj);
                %????
                
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
%                 [s_ft,echo3]=BKganraoxinhao(D,y,R1,tr_B,ts,c,Aj,N,frame_B,fs_B,f_doppler1,tm_B,f0_B,tau_B,congmubiao,y1);
                Pn=(Bj_tp/2/(2.5*Kfm))^2;
                [noise_tp]=zaoshengtiaopin(fs_tp,Kfm,Prj_tp,Pn,Bn,fj_tp,frame_tp,Tr_tp);
                s_echo_1=s_echo_2+s_noise+noise_tp;%%%%????????
                t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
               [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
               [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
               figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????????????');

            case '????????????'
                    load data/data_tf
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                [noise_tf] =zaoshengtiaofu(fs_tf,Bj_tf,fj_tf,Prj_tf,Tr_tf,frame_tf);
                s_echo_1=s_echo_2+s_noise+noise_tf;%%%%????????
                t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
               [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
               [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
               figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????????????');

            case '????????????'
                     load data/data_tx
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                [noise_tx] =zaoshengtiaoxiang(fs_tx,Bj_tx,fj_tx,Prj_tx,Tr_tx,frame_tx);
                s_echo_1=s_echo_2+s_noise+noise_tx;%%%%????????
                t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
               [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
               [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
               figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????????????');
            case '????????????'
                load data/data_smart
                Kfm=4e6;tau=1e-6; c=3e8;
                 Pn=(Bj_smart/(2.5*Kfm))^2;Bn=Bj_smart/2;
                fr=1/tr_B;  
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                [vSmartNoiseSig]=jam_smartnoise( s_echo_2,Pn,Prj_smart,Bn,Kfm,fs_smart );
                s_echo_1=s_echo_2+s_noise+vSmartNoiseSig;%%%%????????
                t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
               [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
               [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
               figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????????????');
           case '??????????'
                load data/data_shuzhuangpu
                 fj=[0.4e6,0.8e6,1.2e6];
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                [noise_szp,~] = jam_combspectrum(fs_shuzhuangpu,Bj_shuzhuangpu,Ns_shuzhuangpu,fj,frame_shuzhuangpu,Prj_shuzhuangpu,Tr_shuzhuangpu);
                s_echo_1=s_echo_2+s_noise+noise_szp;%%%%????????
                t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
               [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
               [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
               figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????????????');
           case '????????'
                load data/data_saopin
                T_fr=2*Tr_saopin;
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
%                 Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
%                 Aj=sqrt(Prj);
                %????
                
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
%                 [s_ft,echo3]=BKganraoxinhao(D,y,R1,tr_B,ts,c,Aj,N,frame_B,fs_B,f_doppler1,tm_B,f0_B,tau_B,congmubiao,y1);
                [sig_noise,~] = jam_sweepfrequency( fs_saopin,Bj_saopin,fj_saopin,frame_saopin,Prj_saopin,Tr_saopin,T_fr,Time_begin_saopin,K_sweep_saopin );             
                s_echo_1=s_echo_2+s_noise+sig_noise;%%%%????????
                t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
               [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
               [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
               figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????????????');

                end
        
    case '????????????'
        load data/data_JDParameter;
        load data/data_target0Parameter
        c=3e8;
        switch strOPSet{nOPSet}
            case '????????????'
                load data/data_sp
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                [noise_sp] = shepinzaosheng(fs_sp,Bj_sp,fj_sp,frame_sp,Prj_sp,Tr_sp);
                s_echo_1=s_echo_2+s_noise+noise_sp;%%%%????????
                t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD;
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????????????');
           
            case '????????????'
                load data/data_tp
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                Pn=(Bj_tp/2/(2.5*Kfm))^2;
                [noise_tp]=zaoshengtiaopin(fs_tp,Kfm,Prj_tp,Pn,Bn,fj_tp,frame_tp,Tr_tp);
                s_echo_1=s_echo_2+s_noise+noise_tp;%%%%????????
                t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD;
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????????????');

            case '????????????'
                load data/data_tf
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                [noise_tf] =zaoshengtiaofu(fs_tf,Bj_tf,fj_tf,Prj_tf,Tr_tf,frame_tf);
                s_echo_1=s_echo_2+s_noise+noise_tf;%%%%????????
                t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD;
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????????????');

            case '????????????'
                load data/data_tx
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                [noise_tx] =zaoshengtiaoxiang(fs_tx,Bj_tx,fj_tx,Prj_tx,Tr_tx,frame_tx);
                s_echo_1=s_echo_2+s_noise+noise_tx;%%%%????????
                t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD;
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????????????');

            case '????????????'
                load data/data_smart
                Kfm=4e6;
                Pn=(Bj_smart/(2.5*Kfm))^2;Bn=Bj_smart/2;
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                [vSmartNoiseSig]=jam_smartnoise( s_echo_2,Pn,Prj_smart,Bn,Kfm,fs_smart );
                s_echo_1=s_echo_2+s_noise+vSmartNoiseSig;%%%%????????
                t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD;
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????????????');
               
            case '??????????'
                load data/data_shuzhuangpu
                fj=[0.4e6,0.8e6,1.2e6];
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                [noise_szp,~] = jam_combspectrum(fs_shuzhuangpu,Bj_shuzhuangpu,Ns_shuzhuangpu,fj,frame_shuzhuangpu,Prj_shuzhuangpu,Tr_shuzhuangpu);
                s_echo_1=s_echo_2+s_noise+noise_szp;%%%%????????
                t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD;
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????????????');
                
            case '????????'
                load data/data_saopin
                T_fr=2*Tr_saopin;
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                [sig_noise,~] = jam_sweepfrequency( fs_saopin,Bj_saopin,fj_saopin,frame_saopin,Prj_saopin,Tr_saopin,T_fr,Time_begin_saopin,K_sweep_saopin );               
                s_echo_1=s_echo_2+s_noise+sig_noise;%%%%????????
                t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD;
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????????????');

        end
    end
elseif n0==4
     str=get(handles.radarSelect,'string');
     n=get(handles.radarSelect,'Value');

     switch  str{n}
     case '????????????'
        c=3e8;
        load data/data_LFMParameter;
        load data/data_target0Parameter
        switch strO{nO}
            case '??????????????'
                load data/data_dopplerblink
                Kfm=4e6;tau=1e-6; c=3e8;
                lamta=c/fz_LFM;%????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                R0=rand(100);
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                f_doppler=2*v/lamta;%????????????????%????????????????????  
                [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
%                 [s_ft,echo3]=LFMganraoxinhao(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao);
                [ sig_jam,~ ] = jam_dopplerblink( fd_dopplerblink,Td_dopplerblink,R0,s_echo_2,fs_dopplerblink,Pj_dopplerblink,flagT_dopplerblink );
                s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('LFM????????');
            
            case '????????'
                load data/data_botiao
                Kfm=4e6;tau=1e-6; c=3e8;
                lamta=c/fz_LFM;%????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                R0=rand(100);
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                f_doppler=2*v/lamta;%????????????????%????????????????????  
                [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                
                 %????????????
                px=1e3;py=1e3;pz=1e3;%????????
                vx=10;vy=10;vz=0;%????????
                ax=0;ay=0;az=0;%??????????
                phi=pi/180;
                [ TargetStatus ] = paraset_targetstatus( px,py,pz,vx,vy,vz,ax,ay,az,phi );
                [ PassivePara ] = paraset_passivejaming( tf_botiao,sf_botiao,vl_botiao,vf_botiao,ts_botiao,bt_botiao,al_botiao,sref_botiao,smax_botiao );
                  vx=1;vy=1;
                [ WindV ] = paraset_windvelocity( vx,vy );
                CurrentT=2.2;
                [ sig_jam,~ ] = jam_passive( s_echo_2,fs_LFM,f0_LFM,CurrentT,TargetStatus,WindV,PassivePara );
                s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('LFM????????');
   
            case 'AGC????' 
                load data/data_AGC
                Kfm=4e6;tau=1e-6; c=3e8;
                lamta=c/fz_LFM;%????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                R0=rand(100);
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                f_doppler=2*v/lamta;%????????????????%????????????????????  
                [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                [ sig_jam,~ ] = jam_AGC( CurrentT_AGC,Pj_AGC,Period_AGC,D_AGC,radio_AGC,s_echo_2,fs_AGC );
                s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('LFM????????');
            
        end
        
    case '????????????'
        load data/data_BFParameter;
        load data/data_target0Parameter
        c=3e8;
        switch strO{nO}
            case '??????????????'
                load data/data_dopplerblink
                Kfm=4e6;tau=1e-6; c=3e8;
%                 global sigma;
%                 global str;
%                 global RadarS;
%                 global code;
                fr=1/tr_B;  
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                R0=rand(100);
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                [ sig_jam,~ ] = jam_dopplerblink( fd_dopplerblink,Td_dopplerblink,R0,s_echo_2,fs_dopplerblink,Pj_dopplerblink,flagT_dopplerblink );
                s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
               [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
               [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
               figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????????????');
            case '????????'
                load data/data_botiao
                Kfm=4e6;tau=1e-6; c=3e8;
%                 global sigma;
%                 global str;
%                 global RadarS;
%                 global code;
                fr=1/tr_B;  
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                R0=rand(100);
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                
                %????????????
                px=1e3;py=1e3;pz=1e3;%????????
                vx=10;vy=10;vz=0;%????????
                ax=0;ay=0;az=0;%??????????
                phi=pi/180;
                [ TargetStatus ] = paraset_targetstatus( px,py,pz,vx,vy,vz,ax,ay,az,phi );
                [ PassivePara ] = paraset_passivejaming( tf_botiao,sf_botiao,vl_botiao,vf_botiao,ts_botiao,bt_botiao,al_botiao,sref_botiao,smax_botiao );
                  vx=1;vy=1;
                [ WindV ] = paraset_windvelocity( vx,vy );
                CurrentT=2.2;
                [ sig_jam,~ ] = jam_passive( s_echo_2,fs_B,f0_B,CurrentT,TargetStatus,WindV,PassivePara );
                
                s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
               [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
               [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
               figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????????????');
            case 'AGC????'
                load data/data_AGC
                Kfm=4e6;tau=1e-6; c=3e8;
%                 global sigma;
%                 global str;
%                 global RadarS;
%                 global code;
                fr=1/tr_B;  
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                R0=rand(100);
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                [ sig_jam,~ ] = jam_AGC( CurrentT_AGC,Pj_AGC,Period_AGC,D_AGC,radio_AGC,s_echo_2,fs_AGC );
                s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
               [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
               [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
               figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????????????');
                
        end
        
        
                
     case '????????????'
        load data/data_JDParameter;
        load data/data_target0Parameter

        c=3e8;
        switch strO{nO}
            case '??????????????'
                load data/data_dopplerblink
                Kfm=4e6;
                R0=rand(100);
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                [ sig_jam,~ ] = jam_dopplerblink( fd_dopplerblink,Td_dopplerblink,R0,s_echo_2,fs_dopplerblink,Pj_dopplerblink,flagT_dopplerblink );
                s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????????????');
            
            case '????????'
                load data/data_botiao
                Kfm=4e6;
                R0=rand(100);
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                %????????????
                px=1e3;py=1e3;pz=1e3;%????????
                vx=10;vy=10;vz=0;%????????
                ax=0;ay=0;az=0;%??????????
                phi=pi/180;
                [ TargetStatus ] = paraset_targetstatus( px,py,pz,vx,vy,vz,ax,ay,az,phi );
                [ PassivePara ] = paraset_passivejaming( tf_botiao,sf_botiao,vl_botiao,vf_botiao,ts_botiao,bt_botiao,al_botiao,sref_botiao,smax_botiao );
                  vx=1;vy=1;
                [ WindV ] = paraset_windvelocity( vx,vy );
                CurrentT=2.2;
                [ sig_jam,~ ] = jam_passive( s_echo_2,fs_JD,f0_JD,CurrentT,TargetStatus,WindV,PassivePara );
                s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????????????');
            
                
            case 'AGC????'
                load data/data_AGC
                Kfm=4e6;
                R0=rand(100);
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                [ sig_jam,~ ] = jam_AGC( CurrentT_AGC,Pj_AGC,Period_AGC,D_AGC,radio_AGC,s_echo_2,fs_AGC );
                s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                [~,pc_result1,~]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                figure,plot(0:ts:(M-1)*ts,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????????????');
            
        end
        
     end
end

% --- Executes on button press in Jmtd.
function Jmtd_Callback(~, ~, handles)%??????MTD
% hObject    handle to Jmtd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sigma;
global rcsk;
global sigma0;
sigma = rcs(rcsk,sigma0);
global temp1;
global string2;
global temp2;
global strO;
global nO;
global strOPSet;
global nOPSet;
str0=get(handles.jammingSelect,'string');
n0=get(handles.jammingSelect,'Value');
global str;
if n0==3
str=get(handles.radarSelect,'string');
n=get(handles.radarSelect,'Value');

global code;
       c=3e8;
       load data/data_target0Parameter
       load data/data_DeceptionJammingParameter
    switch str{n}
        case '????????????'      
             load data/data_LFMParameter
%             global temp1;
%             global string2;
%             global temp2;
            switch string2{temp2}
            case {'??' ,'??','????' }
                fr=1/tr_LFM;%??????????
                lamta=c/fz_LFM;%????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????%????????????????????  
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
                Aj=sqrt(Prj);
                f_doppler1=2*v1/lamta;
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                [s_ft,~]=LFMDeceptionJamming(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao,temp1);
                s_echo_1=s_echo_2+s_noise+s_ft;%%%%????????
                [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
                figure,mesh(1:fr/num_tongdao_LFM:fr,0:ts*(tau_LFM/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_LFM/D/ts)-ts*(tau_LFM/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('LFM????????MTD????');
                
            case '????'
                fr=1/tr_LFM;%??????????
                f1=10e6; %????????
                c=3e8;
                lamta=c/fz_LFM;
                Rj=2e3;
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L_LFM/10)));
                Aj=sqrt(Prj);
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                f_doppler=2*v/lamta;
                f_doppler1=f_doppler;
                k=B1_LFM/tau_LFM;%????????????????????
                ts=1/fs_LFM;%????????
                %????????????
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L/10))); %????????????????
                A=sqrt(Prs);
                Te=1143;%????
                F=4;%????????
                 tf=10;
                An=10*log10((1.382e-23)*Te*B_LFM*10^(F/10));%????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k); 
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                %????????????
               [s_echo_2,~]=LFMtuoyinhuiboxinhao(y,R,tr_LFM,ts,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
              %figure,plot(0:ts:(N-1)*ts,real(s_echo_2(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????');
                %figure,plot((0:fs/length(abs(fft(abs(fftshift(echo(1,:)))))):fs-fs/length(abs(fft(abs(fftshift(echo(1,:))))))),abs((fft(echo(1,:))))),xlabel('????f(??????Hz)'), ylabel('y(????????)'),title('??????????????');

                %????????????
                [s_ft,~]=LFMtuoyinganrao(R,frame_LFM,tf,Aj,tm,tau_LFM,f0_LFM,B1_LFM,N,f_doppler1,f_doppler,vf,lamta,k,ts,c,fs_LFM,temp1);

                %????????
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                %??????????????????????????????????????????????????
                s_echo_1=s_echo_2+s_noise+s_ft;
                 %????
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1,fs_LFM,f0_LFM);
                %????????????
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                %????????????????
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                s_pc_result=reshape(pc_result,1,M1*frame_LFM);
                %??????????
                [s_mtd,s_mtd1]=MTD_tuoyin(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
                % figure,mesh(1:fr/num_tongdao*lamta/2:fr*lamta/2,0:ts*(tau/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau/D/ts)-ts*(tau/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('MTD????');
                figure,subplot(121),mesh(1:fr/num_tongdao_LFM*lamta/2:fr*lamta/2,0:ts*(tau_LFM/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_LFM/D/ts)-ts*(tau_LFM/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????/??'),ylabel('??????????????'),zlabel('y(????????)'),title('????????1??16??????????MTD????');
                subplot(122),mesh(1:fr/num_tongdao_LFM*lamta/2:fr*lamta/2,0:ts*(tau_LFM/D/ts)*c/2:(length(abs(s_mtd1(:,1)))*ts*(tau_LFM/D/ts)-ts*(tau_LFM/D/ts))*c/2,abs(s_mtd1)),xlabel('??????????????/??'),ylabel('??????????????'),zlabel('y(????????)'),title('????????113??128??????????MTD????');

            end
        case '????????????' 
            load data/data_BFParameter
            switch string2{temp2}
            case {'??' ,'??','????' }
                load data/data_BFParameter;
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
                Aj=sqrt(Prj);
                %????
                f_doppler1=2*v1/lamta; 
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                [s_ft,~]=BKDeceptionJamming(D,y,R1,tr_B,ts,c,Aj,N,frame_B,fs_B,f_doppler1,tm_B,f0_B,tau_B,congmubiao,y1,temp1);
                s_echo_1=s_echo_2+s_noise+s_ft;%%%%????????
                %figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                %figure,plot((0:fs_B/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_B-fs_B/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('????f(??????Hz)'), ylabel('y(????????)'),title('??????????????????');
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_B,num_tongdao_B);
                figure,mesh(1:fr/num_tongdao_B:fr,0:ts*(tau_B/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_B/D/ts)-ts*(tau_B/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('????????MTD????');
                
            case '????'
                 fr=1/tr_B;
                 Rj=2e3;
                 f1=10e6;
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                tm=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm);%????????????????????????????
                f_doppler=2*v/lamta;%????????????????
                f_doppler1= f_doppler;%????????????????
                Te=1143;%????
                An=10*log10((1.382e-23)*Te*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L/10)));
                Aj=sqrt(Prj);
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                %             figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(????????)'), ylabel('y(????????)'),title('????????');
                %             figure,plot((0:fs_B/length(abs(fft(abs(fftshift(y))))):fs_B-fs_B/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('????f(??????Hz)'), ylabel('y(????????)'),title('??????????????');
                %             %%%%%%%%%%%1.1.2????????????????????
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                 %%%%%%%%%%%%%1.3????????????%%%%%%%%%%%%
                [s_echo_2,~]=BKtuoyinhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tm,tau_B); 
                %????????
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                %????????
                [s_ft,~]=BKtuoyinganrao(R,frame_B,tf,Aj,tau_B,N,f_doppler1,f_doppler,vf,lamta,ts,c,fs_B,y1,tr_B,temp1);
                %??????????????????????????????????????????????????
                s_echo_1=s_echo_2+s_noise+s_ft;
                %????
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1,fs_B,f0_B);
                %????????????
                [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
                %????????????????
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
                s_pc_result=reshape(pc_result,1,M1*frame_B);
                %??????????
                [s_mtd,s_mtd1]=MTD_tuoyin(pc_result1.',M1,num_jilei_B,num_tongdao_B);
                figure,subplot(121),mesh(1:fr/num_tongdao_B*lamta/2:fr*lamta/2,0:ts*(tau_B/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_B/D/ts)-ts*(tau_B/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????/??'),ylabel('??????????????'),zlabel('y(????????)'),title('????????1??16??????????MTD????');
                subplot(122),mesh(1:fr/num_tongdao_B*lamta/2:fr*lamta/2,0:ts*(tau_B/D/ts)*c/2:(length(abs(s_mtd1(:,1)))*ts*(tau_B/D/ts)-ts*(tau_B/D/ts))*c/2,abs(s_mtd1)),xlabel('??????????????/??'),ylabel('??????????????'),zlabel('y(????????)'),title('????????113??128??????????MTD????');

            end
       case '????????????'
           
           load data/data_JDParameter;
            f1=10e6; %????????
             c=3e8;
             Rj=2e3;
        switch string2{temp2}
            case {'??' ,'??','????' }
                fr=1/tr_JD;
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
                Aj=sqrt(Prj);
                %????
                 f_doppler1=2*v1/lamta; 
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                [s_ft,~]=JDDeceptionJamming(D,R1,tr_JD,c,Aj,N,frame_JD,fs_JD,f_doppler1,tm,f0_JD,tau_JD,congmubiao,temp1);
                s_echo_1=s_echo_2+s_noise+s_ft;%%%%????????
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0_JD,fs_JD);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
                figure,mesh(1:fr/num_tongdao_JD:fr,0:ts*(tau_JD/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_JD/D/ts)-ts*(tau_JD/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('????????MTD????');
                
            case '????'
                fr=1/tr_JD;
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????

                %????
                f_doppler=2*v/lamta;
                f_doppler1=f_doppler;
                lamta=c/fz_JD;
                f_doppler=2*v/lamta;
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                ts=1/fs_JD;%????????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????

                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L/10)));
                Aj=sqrt(Prj);
                %figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(????????)'), ylabel('y(????????)'),title('????????');
                %figure,plot((0:fs/length(abs(fft(abs(fftshift(y))))):fs-fs/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('????f(??????Hz)'), ylabel('y(????????)'),title('??????????????');
                %????????????????????
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);

                %????????????
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                     [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                %????????????
                [s_ft,~]=JDtuoyinganrao(R,frame_JD,tf,Aj,tm,tau_JD,f0_JD,N,f_doppler1,f_doppler,vf,lamta,ts,c,fs_JD,temp1);

                %??????????????????????????????????????????????????
                s_echo_1=s_echo_2+s_noise+s_ft;
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1,fs_JD,f0_JD);
                %????????????
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                %????????????????
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                s_pc_result=reshape(pc_result,1,M1*frame_JD);
                %??????????
                [s_mtd,s_mtd1]=MTD_tuoyin(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
                % figure,mesh(1:fr/num_tongdao*lamta/2:fr*lamta/2,0:ts*(tau/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau/D/ts)-ts*(tau/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('MTD????');
                figure, subplot(121),mesh(1:fr/num_tongdao_JD*lamta/2:fr*lamta/2,0:ts*(tau_JD/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_JD/D/ts)-ts*(tau_JD/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????/??'),ylabel('??????????????'),zlabel('y(????????)'),title('????????1??16??????????MTD????');
                      subplot(122),mesh(1:fr/num_tongdao_JD*lamta/2:fr*lamta/2,0:ts*(tau_JD/D/ts)*c/2:(length(abs(s_mtd1(:,1)))*ts*(tau_JD/D/ts)-ts*(tau_JD/D/ts))*c/2,abs(s_mtd1)),xlabel('??????????????/??'),ylabel('??????????????'),zlabel('y(????????)'),title('????????113??128??????????MTD????');

        end
  
    end
elseif n0==2
    str=get(handles.radarSelect,'string');
    n=get(handles.radarSelect,'Value');

    switch  str{n}
        case '????????????'
                 c=3e8; 
                 load data/data_LFMParameter;
                 load data/data_target0Parameter
%                 global strOPSet;
%                 global nOPSet;
                switch strOPSet{nOPSet}
            case '????????????'  
                load data/data_sp
                fr=1/tr_LFM;
                lamta=c/fz_LFM;%????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                f_doppler=2*v/lamta;%????????????????%????????????????????  
                [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                [noise_sp] = shepinzaosheng(fs_sp,Bj_sp,fj_sp,frame_sp,Prj_sp,Tr_sp);
                s_echo_1=s_echo_2+s_noise+noise_sp;%%%%????????
                t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
                figure,mesh(1:fr/num_tongdao_LFM:fr,0:ts*(tau_LFM/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_LFM/D/ts)-ts*(tau_LFM/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('LFM????????MTD????');
  
            case '????????????'
                load data/data_tx
                fr=1/tr_LFM;
                lamta=c/fz_LFM;%????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                f_doppler=2*v/lamta;%????????????????%????????????????????  
                [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
%                 [s_ft,echo3]=LFMganraoxinhao(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao);
                [noise_tx] =zaoshengtiaoxiang(fs_tx,Bj_tx,fj_tx,Prj_tx,Tr_tx,frame_tx);
                s_echo_1=s_echo_2+s_noise+noise_tx;%%%%????????
                t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
                figure,mesh(1:fr/num_tongdao_LFM:fr,0:ts*(tau_LFM/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_LFM/D/ts)-ts*(tau_LFM/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('LFM????????MTD????');
          

            case '????????????'
                load data/data_tf
                fr=1/tr_LFM;
                lamta=c/fz_LFM;%????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                f_doppler=2*v/lamta;%????????????????%????????????????????  
                [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
%                 [s_ft,echo3]=LFMganraoxinhao(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao);
                [noise_tf] =zaoshengtiaofu(fs_tf,Bj_tf,fj_tf,Prj_tf,Tr_tf,frame_tf);
                s_echo_1=s_echo_2+s_noise+noise_tf;%%%%????????
                t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
                figure,mesh(1:fr/num_tongdao_LFM:fr,0:ts*(tau_LFM/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_LFM/D/ts)-ts*(tau_LFM/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('LFM????????MTD????');
          
             case '????????????'
                 load data/data_tp
                 fr=1/tr_LFM;
                 Pn=(Bj_tp/2/(2.5*Kfm))^2;
                lamta=c/fz_LFM;%????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                f_doppler=2*v/lamta;%????????????????%????????????????????  
                [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
%                 [s_ft,echo3]=LFMganraoxinhao(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao);
                [noise_tp] =zaoshengtiaopin(fs_tp,Kfm,Prj_tp,Pn,Bn,fj_tp,frame_tp,Tr_tp);
                s_echo_1=s_echo_2+s_noise+noise_tp;%%%%????????
                t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
                figure,mesh(1:fr/num_tongdao_LFM:fr,0:ts*(tau_LFM/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_LFM/D/ts)-ts*(tau_LFM/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('LFM????????MTD????');
          
            case '????????????'
                load data/data_smart
                    fr=1/tr_LFM;
                    Kfm=4e6;tau=1e-6; c=3e8;
                   Pn=(B_LFM/(2.5*Kfm))^2;Bn=B_LFM/2;
                lamta=c/fz_LFM;%????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                f_doppler=2*v/lamta;%????????????????%????????????????????  
                [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
%                 [s_ft,echo3]=LFMganraoxinhao(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao);
                [vSmartNoiseSig]=jam_smartnoise( s_echo_2,Pn,Prj_smart,Bn,Kfm,fs_smart );
                s_echo_1=s_echo_2+s_noise+vSmartNoiseSig;%%%%????????
                t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
                figure,mesh(1:fr/num_tongdao_LFM:fr,0:ts*(tau_LFM/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_LFM/D/ts)-ts*(tau_LFM/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('LFM????????MTD????');
           case '??????????'
                   load data/data_shuzhuangpu
                    fj=[0.4e6,0.8e6,1.2e6];
                fr=1/tr_LFM;
                lamta=c/fz_LFM;%????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                f_doppler=2*v/lamta;%????????????????%????????????????????  
                [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
%                 [s_ft,echo3]=LFMganraoxinhao(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao);
               [noise_szp,~] = jam_combspectrum(fs_shuzhuangpu,Bj_shuzhuangpu,Ns_shuzhuangpu,fj,frame_shuzhuangpu,Prj_shuzhuangpu,Tr_shuzhuangpu);
                s_echo_1=s_echo_2+s_noise+noise_szp;%%%%????????
                t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
               [s_mtd]=mtd(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
                figure,mesh(1:fr/num_tongdao_LFM:fr,0:ts*(tau_LFM/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_LFM/D/ts)-ts*(tau_LFM/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('LFM????????MTD????');
            case '????????'
                load data/data_saopin
                 fr=1/tr_LFM;
                T_fr=2*Tr_saopin;
                lamta=c/fz_LFM;%????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                f_doppler=2*v/lamta;%????????????????%????????????????????  
                [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
%                 [s_ft,echo3]=LFMganraoxinhao(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao);
                [ sig_noise,~ ] = jam_sweepfrequency( fs_saopin,Bj_saopin,fj_saopin,frame_saopin,Prj_saopin,Tr_saopin,T_fr,Time_begin_saopin,K_sweep_saopin );
                s_echo_1=s_echo_2+s_noise+sig_noise;%%%%????????
                t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                 [s_mtd]=mtd(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
                figure,mesh(1:fr/num_tongdao_LFM:fr,0:ts*(tau_LFM/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_LFM/D/ts)-ts*(tau_LFM/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('LFM????????MTD????');
            
                
                end 
  
          
    case '????????????'
        load data/data_BFParameter;
        load data/data_target0Parameter

        c=3e8;
        
        switch strOPSet{nOPSet}
            case '????????????'
                load data/data_sp
                fr=1/tr_B;  
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
%                 Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
%                 Aj=sqrt(Prj);
                %????
                
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
%                 [s_ft,echo3]=BKganraoxinhao(D,y,R1,tr_B,ts,c,Aj,N,frame_B,fs_B,f_doppler1,tm_B,f0_B,tau_B,congmubiao,y1);
                [noise_sp] = shepinzaosheng(fs_sp,Bj_sp,fj_sp,frame_sp,Prj_sp,Tr_sp);
                s_echo_1=s_echo_2+s_noise+noise_sp;%%%%????????
                t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_B,num_tongdao_B);
                figure,mesh(1:fr/num_tongdao_B:fr,0:ts*(tau_B/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_B/D/ts)-ts*(tau_B/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('????????MTD????');

            case '????????????'
                load data/data_tp
                fr=1/tr_B;  
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
%                 Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
%                 Aj=sqrt(Prj);
                %????
                
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
%                 [s_ft,echo3]=BKganraoxinhao(D,y,R1,tr_B,ts,c,Aj,N,frame_B,fs_B,f_doppler1,tm_B,f0_B,tau_B,congmubiao,y1);
                Pn=(Bj_tp/2/(2.5*Kfm))^2;
                [noise_tp]=zaoshengtiaopin(fs_tp,Kfm,Prj_tp,Pn,Bn,fj_tp,frame_tp,Tr_tp);
                s_echo_1=s_echo_2+s_noise+noise_tp;%%%%????????
                t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
               [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
               [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
               [s_mtd]=mtd(pc_result1.',M1,num_jilei_B,num_tongdao_B);
                figure,mesh(1:fr/num_tongdao_B:fr,0:ts*(tau_B/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_B/D/ts)-ts*(tau_B/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('????????MTD????');

            case '????????????'
                    load data/data_tf
                    fr=1/tr_B;  
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                [noise_tf] =zaoshengtiaofu(fs_tf,Bj_tf,fj_tf,Prj_tf,Tr_tf,frame_tf);
                s_echo_1=s_echo_2+s_noise+noise_tf;%%%%????????
                t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_B,num_tongdao_B);
                figure,mesh(1:fr/num_tongdao_B:fr,0:ts*(tau_B/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_B/D/ts)-ts*(tau_B/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('????????MTD????');

            case '????????????'
                load data/data_tx
                fr=1/tr_B;  
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                [noise_tx] =zaoshengtiaoxiang(fs_tx,Bj_tx,fj_tx,Prj_tx,Tr_tx,frame_tx);
                s_echo_1=s_echo_2+s_noise+noise_tx;%%%%????????
                t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_B,num_tongdao_B);
                figure,mesh(1:fr/num_tongdao_B:fr,0:ts*(tau_B/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_B/D/ts)-ts*(tau_B/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('????????MTD????');
           case '????????????'
                load data/data_smart
                Kfm=4e6;tau=1e-6; c=3e8;
%                 global sigma;
%                 global str;
%                 global RadarS;
%                 global code;
                 Pn=(Bj_smart/(2.5*Kfm))^2;Bn=Bj_smart/2;
                fr=1/tr_B;  
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                [vSmartNoiseSig]=jam_smartnoise( s_echo_2,Pn,Prj_smart,Bn,Kfm,fs_smart );
                s_echo_1=s_echo_2+s_noise+vSmartNoiseSig;%%%%????????
                t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
               [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
               [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_B,num_tongdao_B);
                figure,mesh(1:fr/num_tongdao_B:fr,0:ts*(tau_B/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_B/D/ts)-ts*(tau_B/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('????????MTD????');
            case '??????????'
                load data/data_shuzhuangpu
                 fj=[0.4e6,0.8e6,1.2e6];
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                [noise_szp,~] = jam_combspectrum(fs_shuzhuangpu,Bj_shuzhuangpu,Ns_shuzhuangpu,fj,frame_shuzhuangpu,Prj_shuzhuangpu,Tr_shuzhuangpu);
                s_echo_1=s_echo_2+s_noise+noise_szp;%%%%????????
                t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_B,num_tongdao_B);
                figure,mesh(1:fr/num_tongdao_B:fr,0:ts*(tau_B/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_B/D/ts)-ts*(tau_B/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('????????MTD????');
            case '????????'
                load data/data_saopin
                T_fr=2*Tr_saopin;
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
%                 Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
%                 Aj=sqrt(Prj);
                %????
                
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
%                 [s_ft,echo3]=BKganraoxinhao(D,y,R1,tr_B,ts,c,Aj,N,frame_B,fs_B,f_doppler1,tm_B,f0_B,tau_B,congmubiao,y1);
                [sig_noise,~] = jam_sweepfrequency( fs_saopin,Bj_saopin,fj_saopin,frame_saopin,Prj_saopin,Tr_saopin,T_fr,Time_begin_saopin,K_sweep_saopin );             
                s_echo_1=s_echo_2+s_noise+sig_noise;%%%%????????
                t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
               [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
               [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
              [s_mtd]=mtd(pc_result1.',M1,num_jilei_B,num_tongdao_B);
                figure,mesh(1:fr/num_tongdao_B:fr,0:ts*(tau_B/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_B/D/ts)-ts*(tau_B/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('????????MTD????');

         end
        
    case '????????????'
        load data/data_JDParameter;
        load data/data_target0Parameter
        c=3e8;

        switch strOPSet{nOPSet}
            case '????????????'
                load data/data_sp
                lamta=c/fz_JD;%????
                fr=1/tr_JD;
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                [noise_sp] = shepinzaosheng(fs_sp,Bj_sp,fj_sp,frame_sp,Prj_sp,Tr_sp);
                s_echo_1=s_echo_2+s_noise+noise_sp;%%%%????????
                t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD;
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
                figure,mesh(1:fr/num_tongdao_JD:fr,0:ts*(tau_JD/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_JD/D/ts)-ts*(tau_JD/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('????????MTD????');
               
            case '????????????'
                load data/data_tp
                fr=1/tr_JD;
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                Pn=(Bj_tp/2/(2.5*Kfm))^2;
                [noise_tp]=zaoshengtiaopin(fs_tp,Kfm,Prj_tp,Pn,Bn,fj_tp,frame_tp,Tr_tp);
                s_echo_1=s_echo_2+s_noise+noise_tp;%%%%????????
                t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD;
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
                figure,mesh(1:fr/num_tongdao_JD:fr,0:ts*(tau_JD/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_JD/D/ts)-ts*(tau_JD/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('????????MTD????');

            case '????????????'
                load data/data_tf
                fr=1/tr_JD;
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                [noise_tf] =zaoshengtiaofu(fs_tf,Bj_tf,fj_tf,Prj_tf,Tr_tf,frame_tf);
                s_echo_1=s_echo_2+s_noise+noise_tf;%%%%????????
                t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD;
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
                figure,mesh(1:fr/num_tongdao_JD:fr,0:ts*(tau_JD/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_JD/D/ts)-ts*(tau_JD/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('????????MTD????');

            case '????????????'
                load data/data_tx
                fr=1/tr_JD;
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                [noise_tx] =zaoshengtiaoxiang(fs_tx,Bj_tx,fj_tx,Prj_tx,Tr_tx,frame_tx);
                s_echo_1=s_echo_2+s_noise+noise_tx;%%%%????????
                t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD;
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
                figure,mesh(1:fr/num_tongdao_JD:fr,0:ts*(tau_JD/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_JD/D/ts)-ts*(tau_JD/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('????????MTD????');
               
            case '????????????'
                load data/data_smart
                fr=1/tr_JD;
                Kfm=4e6;
                Pn=(Bj_smart/(2.5*Kfm))^2;Bn=Bj_smart/2;
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                [vSmartNoiseSig]=jam_smartnoise( s_echo_2,Pn,Prj_smart,Bn,Kfm,fs_smart );
                s_echo_1=s_echo_2+s_noise+vSmartNoiseSig;%%%%????????
                t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD;
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
                figure,mesh(1:fr/num_tongdao_JD:fr,0:ts*(tau_JD/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_JD/D/ts)-ts*(tau_JD/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('????????MTD????');
                
            case '??????????'
                load data/data_shuzhuangpu
                fj=[0.4e6,0.8e6,1.2e6];
                fr=1/tr_JD;
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                [noise_szp,~] = jam_combspectrum(fs_shuzhuangpu,Bj_shuzhuangpu,Ns_shuzhuangpu,fj,frame_shuzhuangpu,Prj_shuzhuangpu,Tr_shuzhuangpu);
                s_echo_1=s_echo_2+s_noise+noise_szp;%%%%????????
                t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD;
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
                figure,mesh(1:fr/num_tongdao_JD:fr,0:ts*(tau_JD/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_JD/D/ts)-ts*(tau_JD/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('????????MTD????');
           
            case '????????'
                load data/data_saopin
                T_fr=2*Tr_saopin;
                fr=1/tr_JD;
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                [sig_noise,~] = jam_sweepfrequency( fs_saopin,Bj_saopin,fj_saopin,frame_saopin,Prj_saopin,Tr_saopin,T_fr,Time_begin_saopin,K_sweep_saopin );               
                s_echo_1=s_echo_2+s_noise+sig_noise;%%%%????????
                t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD;
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                  [s_mtd]=mtd(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
                figure,mesh(1:fr/num_tongdao_JD:fr,0:ts*(tau_JD/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_JD/D/ts)-ts*(tau_JD/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('????????MTD????');

        end
    end
elseif n0==4
    str=get(handles.radarSelect,'string');
    n=get(handles.radarSelect,'Value');
    switch  str{n}
        case '????????????'
            c=3e8; 
            load data/data_LFMParameter;
            load data/data_target0Parameter

            switch strO{nO}
                case '??????????????'
                    load data/data_dopplerblink
                    R0=rand(100);
                    fr=1/tr_LFM;
                    lamta=c/fz_LFM;%????
                    tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                    N=length(tm);%????????????????????????????
                    An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                    ts=1/fs_LFM;
                    k=B1_LFM/tau_LFM;   
                    Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                    A=sqrt(Prs);%????????????
                    [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                    [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                    f_doppler=2*v/lamta;%????????????????%????????????????????  
                    [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                    [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                    [ sig_jam,~ ] = jam_dopplerblink( fd_dopplerblink,Td_dopplerblink,R0,s_echo_2,fs_dopplerblink,Pj_dopplerblink,flagT_dopplerblink );
                    s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                    [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                    [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                    [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                    [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                     [s_mtd]=mtd(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
                    figure,mesh(1:fr/num_tongdao_LFM:fr,0:ts*(tau_LFM/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_LFM/D/ts)-ts*(tau_LFM/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('LFM????????MTD????');
            
                case '????????'
                    load data/data_botiao  
                     fr=1/tr_LFM;
                    lamta=c/fz_LFM;%????
                    tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                    N=length(tm);%????????????????????????????
                    An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                    ts=1/fs_LFM;
                    k=B1_LFM/tau_LFM;   
                    Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                    A=sqrt(Prs);%????????????
                    [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                    [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                    f_doppler=2*v/lamta;%????????????????%????????????????????  
                    [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                    [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                    px=1e3;py=1e3;pz=1e3;%????????
                    vx=10;vy=10;vz=0;%????????
                    ax=0;ay=0;az=0;%??????????
                    phi=pi/180;
                    [ TargetStatus ] = paraset_targetstatus( px,py,pz,vx,vy,vz,ax,ay,az,phi );
                    [ PassivePara ] = paraset_passivejaming( tf_botiao,sf_botiao,vl_botiao,vf_botiao,ts_botiao,bt_botiao,al_botiao,sref_botiao,smax_botiao );
                      vx=1;vy=1;
                    [ WindV ] = paraset_windvelocity( vx,vy );
                    CurrentT=2.2;
                    [ sig_jam,~ ] = jam_passive( s_echo_2,fs_LFM,f0_LFM,CurrentT,TargetStatus,WindV,PassivePara );
                    s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                    [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                    [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                    [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                    [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                     [s_mtd]=mtd(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
                    figure,mesh(1:fr/num_tongdao_LFM:fr,0:ts*(tau_LFM/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_LFM/D/ts)-ts*(tau_LFM/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('LFM????????MTD????');
            
                case 'AGC????' 
                    load data/data_AGC 
                     fr=1/tr_LFM;
                    lamta=c/fz_LFM;%????
                    tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                    N=length(tm);%????????????????????????????
                    An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                    ts=1/fs_LFM;
                    k=B1_LFM/tau_LFM;   
                    Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                    A=sqrt(Prs);%????????????
                    [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                    [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                    f_doppler=2*v/lamta;%????????????????%????????????????????  
                    [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                    [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                    [ sig_jam,~ ] = jam_AGC( CurrentT_AGC,Pj_AGC,Period_AGC,D_AGC,radio_AGC,s_echo_2,fs_AGC );
                    s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                    [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                    [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                    [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                    [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                     [s_mtd]=mtd(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
                    figure,mesh(1:fr/num_tongdao_LFM:fr,0:ts*(tau_LFM/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_LFM/D/ts)-ts*(tau_LFM/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('LFM????????MTD????');
            
            end
        case '????????????'
            load data/data_BFParameter;
            load data/data_target0Parameter

            c=3e8; 
            switch strO{nO}
                case '??????????????'
                    load data/data_dopplerblink
                    R0=rand(100);
                    fj=[0.4e6,0.8e6,1.2e6];
                    tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                    N=length(tm_B);%????????????????????????????
                    number1=length(code);
                    ts=1/fs_B;
                    lamta=c/fz_B;
                    fr=1/tr_B;
                    An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                    Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                    A=sqrt(Prs);%????????????
                    f_doppler=2*v/lamta;%????????????????
                    [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                    [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                    [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                    [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                    [ sig_jam,~ ] = jam_dopplerblink( fd_dopplerblink,Td_dopplerblink,R0,s_echo_2,fs_dopplerblink,Pj_dopplerblink,flagT_dopplerblink );
                    s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                    [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                    [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
                    [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
                    [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
                    [s_mtd]=mtd(pc_result1.',M1,num_jilei_B,num_tongdao_B);
                    figure,mesh(1:fr/num_tongdao_B:fr,0:ts*(tau_B/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_B/D/ts)-ts*(tau_B/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('????????MTD????');
            
                case '????????'
                    load data/data_botiao  
                    fj=[0.4e6,0.8e6,1.2e6];
                    tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                    N=length(tm_B);%????????????????????????????
                    number1=length(code);
                    ts=1/fs_B;
                    lamta=c/fz_B;
                    fr=1/tr_B;
                    An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                    Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                    A=sqrt(Prs);%????????????
                    f_doppler=2*v/lamta;%????????????????
                    [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                    [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                    [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                    [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                        %????????????
                    px=1e3;py=1e3;pz=1e3;%????????
                    vx=10;vy=10;vz=0;%????????
                    ax=0;ay=0;az=0;%??????????
                    phi=pi/180;
                    [ TargetStatus ] = paraset_targetstatus( px,py,pz,vx,vy,vz,ax,ay,az,phi );
                    [ PassivePara ] = paraset_passivejaming( tf_botiao,sf_botiao,vl_botiao,vf_botiao,ts_botiao,bt_botiao,al_botiao,sref_botiao,smax_botiao );
                      vx=1;vy=1;
                    [ WindV ] = paraset_windvelocity( vx,vy );
                    CurrentT=2.2;
                    [ sig_jam,~ ] = jam_passive( s_echo_2,fs_B,f0_B,CurrentT,TargetStatus,WindV,PassivePara );
                    s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                    [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                    [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
                    [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
                    [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
                    [s_mtd]=mtd(pc_result1.',M1,num_jilei_B,num_tongdao_B);
                    figure,mesh(1:fr/num_tongdao_B:fr,0:ts*(tau_B/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_B/D/ts)-ts*(tau_B/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('????????MTD????');
            
                case 'AGC????' 
                    load data/data_AGC 
                    fj=[0.4e6,0.8e6,1.2e6];
                    tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                    N=length(tm_B);%????????????????????????????
                    number1=length(code);
                    ts=1/fs_B;
                    lamta=c/fz_B;
                    fr=1/tr_B;
                    An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                    Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                    A=sqrt(Prs);%????????????
                    f_doppler=2*v/lamta;%????????????????
                    [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                    [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                    [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                    [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                    [ sig_jam,~ ] = jam_AGC( CurrentT_AGC,Pj_AGC,Period_AGC,D_AGC,radio_AGC,s_echo_2,fs_AGC );
                    s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                    [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                    [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
                    [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
                    [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
                    [s_mtd]=mtd(pc_result1.',M1,num_jilei_B,num_tongdao_B);
                    figure,mesh(1:fr/num_tongdao_B:fr,0:ts*(tau_B/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_B/D/ts)-ts*(tau_B/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('????????MTD????');
            
            end
            
        case '????????????'
            load data/data_JDParameter;
            load data/data_target0Parameter
            c=3e8;
            switch strO{nO}
                case '??????????????'
                    load data/data_dopplerblink
                    R0=rand(100);
                    fj=[0.4e6,0.8e6,1.2e6];
                    fr=1/tr_JD;
                    lamta=c/fz_JD;%????
                    tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                    N=length(tm);%????????????????????????????
                    An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                    [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                    ts=1/fs_JD; 
                    Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                    A=sqrt(Prs);%????????????
                    f_doppler=2*v/lamta;%????????????????
                    [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                    [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                    [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                    [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                    [ sig_jam,~ ] = jam_dopplerblink( fd_dopplerblink,Td_dopplerblink,R0,s_echo_2,fs_dopplerblink,Pj_dopplerblink,flagT_dopplerblink );
                    s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                    [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                    [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                    [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                    [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                    [s_mtd]=mtd(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
                    figure,mesh(1:fr/num_tongdao_JD:fr,0:ts*(tau_JD/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_JD/D/ts)-ts*(tau_JD/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('????????MTD????');

                case '????????'
                    load data/data_botiao  
                    fj=[0.4e6,0.8e6,1.2e6];
                    fr=1/tr_JD;
                    lamta=c/fz_JD;%????
                    tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                    N=length(tm);%????????????????????????????
                    An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                    [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                    ts=1/fs_JD; 
                    Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                    A=sqrt(Prs);%????????????
                    f_doppler=2*v/lamta;%????????????????
                    [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                    [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                    [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                    [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                        %????????????
                    px=1e3;py=1e3;pz=1e3;%????????
                    vx=10;vy=10;vz=0;%????????
                    ax=0;ay=0;az=0;%??????????
                    phi=pi/180;
                    [ TargetStatus ] = paraset_targetstatus( px,py,pz,vx,vy,vz,ax,ay,az,phi );
                    [ PassivePara ] = paraset_passivejaming( tf_botiao,sf_botiao,vl_botiao,vf_botiao,ts_botiao,bt_botiao,al_botiao,sref_botiao,smax_botiao );
                      vx=1;vy=1;
                    [ WindV ] = paraset_windvelocity( vx,vy );
                    CurrentT=2.2;
                    [ sig_jam,~ ] = jam_passive( s_echo_2,fs_JD,f0_JD,CurrentT,TargetStatus,WindV,PassivePara );
                    s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                    [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                    [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                    [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                    [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                    [s_mtd]=mtd(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
                    figure,mesh(1:fr/num_tongdao_JD:fr,0:ts*(tau_JD/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_JD/D/ts)-ts*(tau_JD/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('????????MTD????');

                case 'AGC????' 
                    load data/data_AGC 
                    fj=[0.4e6,0.8e6,1.2e6];
                    fr=1/tr_JD;
                    lamta=c/fz_JD;%????
                    tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                    N=length(tm);%????????????????????????????
                    An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                    [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                    ts=1/fs_JD; 
                    Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                    A=sqrt(Prs);%????????????
                    f_doppler=2*v/lamta;%????????????????
                    [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                    [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                    [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                    [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                    [ sig_jam,~ ] = jam_AGC( CurrentT_AGC,Pj_AGC,Period_AGC,D_AGC,radio_AGC,s_echo_2,fs_AGC );
                    s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                    [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                    [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                    [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                    [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                    [s_mtd]=mtd(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
                    figure,mesh(1:fr/num_tongdao_JD:fr,0:ts*(tau_JD/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_JD/D/ts)-ts*(tau_JD/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('????????MTD????');

            end     
    end
    
end

% --- Executes on button press in Jhengxujing.
function Jhengxujing_Callback(~, ~, handles)%????????????
% hObject    handle to Jhengxujing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sigma;
global rcsk;
global sigma0;
sigma = rcs(rcsk,sigma0);
global string2;
global temp2;
global temp1;
global strOPSet;
global nOPSet;
global strO;
global nO;
str0=get(handles.jammingSelect,'string');
n0=get(handles.jammingSelect,'Value');
if n0==3
global str;
str=get(handles.radarSelect,'string');
n=get(handles.radarSelect,'Value');

global code;
       c=3e8;
       load data/data_target0Parameter
       load data/data_DeceptionJammingParameter
    switch str{n}
        case '????????????'      
            load data/data_LFMParameter
%             global string2;
%             global temp2;
            switch string2{temp2}
            case {'??' ,'??','????' }
                fr=1/tr_LFM;%??????????
                lamta=c/fz_LFM;%????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????%????????????????????  
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
                Aj=sqrt(Prj);
                f_doppler1=2*v1/lamta;
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                [s_ft,~]=LFMDeceptionJamming(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao,temp1);
                s_echo_1=s_echo_2+s_noise+s_ft;%%%%????????
                [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
                hengxujing(M1,Pfa_LFM,s_mtd,num_cankao_LFM,num_tongdao_LFM,num_baohu_LFM,ts,c,tau_LFM,D);
          
            case '????'
                fr=1/tr_LFM;%??????????
                f1=10e6; %????????
                c=3e8;
                lamta=c/fz_LFM;
                Rj=2e3;
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L_LFM/10)));
                Aj=sqrt(Prj);
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                f_doppler=2*v/lamta;
                f_doppler1=f_doppler;
                k=B1_LFM/tau_LFM;%????????????????????
                ts=1/fs_LFM;%????????
                %????????????
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L/10))); %????????????????
                A=sqrt(Prs);
                Te=1143;%????
                F=4;%????????
                 tf=10;
                An=10*log10((1.382e-23)*Te*B_LFM*10^(F/10));%????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k); 
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                %????????????
               [s_echo_2,~]=LFMtuoyinhuiboxinhao(y,R,tr_LFM,ts,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
              %figure,plot(0:ts:(N-1)*ts,real(s_echo_2(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????');
                %figure,plot((0:fs/length(abs(fft(abs(fftshift(echo(1,:)))))):fs-fs/length(abs(fft(abs(fftshift(echo(1,:))))))),abs((fft(echo(1,:))))),xlabel('????f(??????Hz)'), ylabel('y(????????)'),title('??????????????');

                %????????????
                [s_ft,~]=LFMtuoyinganrao(R,frame_LFM,tf,Aj,tm,tau_LFM,f0_LFM,B1_LFM,N,f_doppler1,f_doppler,vf,lamta,k,ts,c,fs_LFM,temp1);

                %????????
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                %??????????????????????????????????????????????????
                s_echo_1=s_echo_2+s_noise+s_ft;
                 %????
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1,fs_LFM,f0_LFM);
                %????????????
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                %????????????????
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                s_pc_result=reshape(pc_result,1,M1*frame_LFM);
                %??????????
                [s_mtd,s_mtd1]=MTD_tuoyin(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
                % figure,mesh(1:fr/num_tongdao*lamta/2:fr*lamta/2,0:ts*(tau/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau/D/ts)-ts*(tau/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('MTD????');
                hengxujing_tuoyin(M1,Pfa_LFM,s_mtd,num_cankao_LFM,num_tongdao_LFM,num_baohu_LFM,ts,c,tau_LFM,D,s_mtd1);

            end
            case '????????????' 
                load data/data_BFParameter
%                 global string2;
%                 global temp2;
%                 global temp1;
                switch string2{temp2}
                case {'??' ,'??','????' }
                load data/data_BFParameter;
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
                Aj=sqrt(Prj);
                %????
                f_doppler1=2*v1/lamta; 
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                [s_ft,~]=BKDeceptionJamming(D,y,R1,tr_B,ts,c,Aj,N,frame_B,fs_B,f_doppler1,tm_B,f0_B,tau_B,congmubiao,y1,temp1);
                s_echo_1=s_echo_2+s_noise+s_ft;%%%%????????
                %figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                %figure,plot((0:fs_B/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_B-fs_B/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('????f(??????Hz)'), ylabel('y(????????)'),title('??????????????????');
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_B,num_tongdao_B);
                hengxujing(M1,Pfa_B,s_mtd,num_cankao_B,num_tongdao_B,num_baohu_B,ts,c,tau_B,D);
             case '????'
                fr=1/tr_B;
                Rj=2e3;
                f1=10e6;
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                tm=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm);%????????????????????????????
                f_doppler=2*v/lamta;%????????????????
                f_doppler1= f_doppler;%????????????????
                Te=1143;%????
                An=10*log10((1.382e-23)*Te*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L/10)));
                Aj=sqrt(Prj);
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                %             figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(????????)'), ylabel('y(????????)'),title('????????');
                %             figure,plot((0:fs_B/length(abs(fft(abs(fftshift(y))))):fs_B-fs_B/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('????f(??????Hz)'), ylabel('y(????????)'),title('??????????????');
                %             %%%%%%%%%%%1.1.2????????????????????
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                 %%%%%%%%%%%%%1.3????????????%%%%%%%%%%%%
                [s_echo_2,~]=BKtuoyinhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tm,tau_B); 
                %????????
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                %????????
                [s_ft,~]=BKtuoyinganrao(R,frame_B,tf,Aj,tau_B,N,f_doppler1,f_doppler,vf,lamta,ts,c,fs_B,y1,tr_B,temp1);
                %??????????????????????????????????????????????????
                s_echo_1=s_echo_2+s_noise+s_ft;
                %????
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1,fs_B,f0_B);
                %????????????
                [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
                %????????????????
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
                s_pc_result=reshape(pc_result,1,M1*frame_B);
                %??????????
                [s_mtd,s_mtd1]=MTD_tuoyin(pc_result1.',M1,num_jilei_B,num_tongdao_B);
                hengxujing_tuoyin(M1,Pfa_B,s_mtd,num_cankao_B,num_tongdao_B,num_baohu_B,ts,c,tau_B,D,s_mtd1);

                end
        case '????????????'
           
           load data/data_JDParameter;
           f1=10e6; %????????
           c=3e8;
           Rj=2e3;
          
          switch string2{temp2}
            case {'??' ,'??','????' }
                load data/data_JDParameter

                fr=1/tr_JD;
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
                Aj=sqrt(Prj);
                 f_doppler1=2*v1/lamta;
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????%???????????????????? 
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                [s_ft,~]=JDDeceptionJamming(D,R1,tr_JD,c,Aj,N,frame_JD,fs_JD,f_doppler1,tm,f0_JD,tau_JD,congmubiao,temp1);
                s_echo_1=s_echo_2+s_noise+s_ft;
                [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
                hengxujing(M1,Pfa_JD,s_mtd,num_cankao_JD,num_tongdao_JD,num_baohu_JD,ts,c,tau_JD,D);
            case '????'
                fr=1/tr_JD;
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????

                %????
                f_doppler=2*v/lamta;
                f_doppler1=f_doppler;
                lamta=c/fz_JD;
                f_doppler=2*v/lamta;
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                ts=1/fs_JD;%????????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????

                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L/10)));
                Aj=sqrt(Prj);
        %         figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(????????)'), ylabel('y(????????)'),title('????????');
        %         figure,plot((0:fs/length(abs(fft(abs(fftshift(y))))):fs-fs/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('????f(??????Hz)'), ylabel('y(????????)'),title('??????????????');
                %????????????????????
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);

                %????????????
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                     [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                %????????????
                [s_ft,~]=JDtuoyinganrao(R,frame_JD,tf,Aj,tm,tau_JD,f0_JD,N,f_doppler1,f_doppler,vf,lamta,ts,c,fs_JD,temp1);

                %??????????????????????????????????????????????????
                s_echo_1=s_echo_2+s_noise+s_ft;
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1,fs_JD,f0_JD);
                %????????????
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                %????????????????
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                s_pc_result=reshape(pc_result,1,M1*frame_JD);
                %??????????
                [s_mtd,s_mtd1]=MTD_tuoyin(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
                % figure,mesh(1:fr/num_tongdao*lamta/2:fr*lamta/2,0:ts*(tau/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau/D/ts)-ts*(tau/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('MTD????');
                hengxujing_tuoyin(M1,Pfa_JD,s_mtd,num_cankao_JD,num_tongdao_JD,num_baohu_JD,ts,c,tau_JD,D,s_mtd1);

          end
  
    end
elseif n0==2
  str=get(handles.radarSelect,'string');
    n=get(handles.radarSelect,'Value');

    switch  str{n}
        case '????????????'
                 c=3e8; 
                 load data/data_LFMParameter;
                 load data/data_target0Parameter
%                  global strOPSet;
%                  global nOPSet;
            switch strOPSet{nOPSet}
            case '????????????'  
                load data/data_sp
                fr=1/tr_LFM;
                lamta=c/fz_LFM;%????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                f_doppler=2*v/lamta;%????????????????%????????????????????  
                [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                [noise_sp] = shepinzaosheng(fs_sp,Bj_sp,fj_sp,frame_sp,Prj_sp,Tr_sp);
                s_echo_1=s_echo_2+s_noise+noise_sp;%%%%????????
                t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
                 hengxujing(M1,Pfa_LFM,s_mtd,num_cankao_LFM,num_tongdao_LFM,num_baohu_LFM,ts,c,tau_LFM,D);  
            case '????????????'
                load data/data_tx
                fr=1/tr_LFM;
                lamta=c/fz_LFM;%????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                f_doppler=2*v/lamta;%????????????????%????????????????????  
                [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
%                 [s_ft,echo3]=LFMganraoxinhao(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao);
                [noise_tx] =zaoshengtiaoxiang(fs_tx,Bj_tx,fj_tx,Prj_tx,Tr_tx,frame_tx);
                s_echo_1=s_echo_2+s_noise+noise_tx;%%%%????????
                t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
                 hengxujing(M1,Pfa_LFM,s_mtd,num_cankao_LFM,num_tongdao_LFM,num_baohu_LFM,ts,c,tau_LFM,D);  
          

            case '????????????'
                load data/data_tf
                fr=1/tr_LFM;
                lamta=c/fz_LFM;%????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                f_doppler=2*v/lamta;%????????????????%????????????????????  
                [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
%                 [s_ft,echo3]=LFMganraoxinhao(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao);
                [noise_tf] =zaoshengtiaofu(fs_tf,Bj_tf,fj_tf,Prj_tf,Tr_tf,frame_tf);
                s_echo_1=s_echo_2+s_noise+noise_tf;%%%%????????
                t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
                 hengxujing(M1,Pfa_LFM,s_mtd,num_cankao_LFM,num_tongdao_LFM,num_baohu_LFM,ts,c,tau_LFM,D);  
          
             case '????????????'
                 load data/data_tp
                 fr=1/tr_LFM;
                 Pn=(Bj_tp/2/(2.5*Kfm))^2;
                lamta=c/fz_LFM;%????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                f_doppler=2*v/lamta;%????????????????%????????????????????  
                [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
%                 [s_ft,echo3]=LFMganraoxinhao(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao);
                [noise_tp] =zaoshengtiaopin(fs_tp,Kfm,Prj_tp,Pn,Bn,fj_tp,frame_tp,Tr_tp);
                s_echo_1=s_echo_2+s_noise+noise_tp;%%%%????????
                t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
                 hengxujing(M1,Pfa_LFM,s_mtd,num_cankao_LFM,num_tongdao_LFM,num_baohu_LFM,ts,c,tau_LFM,D);  
                 
             case  '????????????'
                load data/data_smart
                fr=1/tr_LFM;
                Kfm=4e6;tau=1e-6; c=3e8;
                Pn=(B_LFM/(2.5*Kfm))^2;Bn=B_LFM/2;
                lamta=c/fz_LFM;%????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                f_doppler=2*v/lamta;%????????????????%????????????????????  
                [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
%                 [s_ft,echo3]=LFMganraoxinhao(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao);
                [vSmartNoiseSig]=jam_smartnoise( s_echo_2,Pn,Prj_smart,Bn,Kfm,fs_smart );
                s_echo_1=s_echo_2+s_noise+vSmartNoiseSig;%%%%????????
                t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
                 hengxujing(M1,Pfa_LFM,s_mtd,num_cankao_LFM,num_tongdao_LFM,num_baohu_LFM,ts,c,tau_LFM,D);
                 
            case '??????????'
                   load data/data_shuzhuangpu
                    fj=[0.4e6,0.8e6,1.2e6];

                lamta=c/fz_LFM;%????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                f_doppler=2*v/lamta;%????????????????%????????????????????  
                [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
%                 [s_ft,echo3]=LFMganraoxinhao(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao);
               [noise_szp,~] = jam_combspectrum(fs_shuzhuangpu,Bj_shuzhuangpu,Ns_shuzhuangpu,fj,frame_shuzhuangpu,Prj_shuzhuangpu,Tr_shuzhuangpu);
                s_echo_1=s_echo_2+s_noise+noise_szp;%%%%????????
                t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
               [s_mtd]=mtd(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
                 hengxujing(M1,Pfa_LFM,s_mtd,num_cankao_LFM,num_tongdao_LFM,num_baohu_LFM,ts,c,tau_LFM,D); 
                 
           case '????????'
                load data/data_saopin
                T_fr=2*Tr_saopin;
                lamta=c/fz_LFM;%????
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                A=sqrt(Prs);%????????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                f_doppler=2*v/lamta;%????????????????%????????????????????  
                [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
%                 [s_ft,echo3]=LFMganraoxinhao(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao);
                [ sig_noise,~ ] = jam_sweepfrequency( fs_saopin,Bj_saopin,fj_saopin,frame_saopin,Prj_saopin,Tr_saopin,T_fr,Time_begin_saopin,K_sweep_saopin );
                s_echo_1=s_echo_2+s_noise+sig_noise;%%%%????????
                t=0:1/fs_LFM:frame_LFM*tr_LFM-1/fs_LFM; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_LFM/2,tau_LFM);
                [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
               [s_mtd]=mtd(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
                hengxujing(M1,Pfa_LFM,s_mtd,num_cankao_LFM,num_tongdao_LFM,num_baohu_LFM,ts,c,tau_LFM,D);  
                
            end 

            
    case '????????????'
                load data/data_BFParameter;
                load data/data_target0Parameter
                
%                 global code;
%                 global sigma;
                c=3e8;
%                 global strOPSet;
%                 global nOPSet;
                switch strOPSet{nOPSet}
            case '????????????'
                load data/data_sp
                fr=1/tr_B;  
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
%                 Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
%                 Aj=sqrt(Prj);
                %????
                
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
%                 [s_ft,echo3]=BKganraoxinhao(D,y,R1,tr_B,ts,c,Aj,N,frame_B,fs_B,f_doppler1,tm_B,f0_B,tau_B,congmubiao,y1);
                [noise_sp] = shepinzaosheng(fs_sp,Bj_sp,fj_sp,frame_sp,Prj_sp,Tr_sp);
                s_echo_1=s_echo_2+s_noise+noise_sp;%%%%????????
                t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
               [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
               [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_B,num_tongdao_B);
                hengxujing(M1,Pfa_B,s_mtd,num_cankao_B,num_tongdao_B,num_baohu_B,ts,c,tau_B,D);
            case '????????????'
                load data/data_tp
                fr=1/tr_B;  
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
%                 Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
%                 Aj=sqrt(Prj);
                %????
                
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
%                 [s_ft,echo3]=BKganraoxinhao(D,y,R1,tr_B,ts,c,Aj,N,frame_B,fs_B,f_doppler1,tm_B,f0_B,tau_B,congmubiao,y1);
                Pn=(Bj_tp/2/(2.5*Kfm))^2;
                [noise_tp]=zaoshengtiaopin(fs_tp,Kfm,Prj_tp,Pn,Bn,fj_tp,frame_tp,Tr_tp);
                s_echo_1=s_echo_2+s_noise+noise_tp;%%%%????????
                t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
               [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
               [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
               [s_mtd]=mtd(pc_result1.',M1,num_jilei_B,num_tongdao_B);
                 hengxujing(M1,Pfa_B,s_mtd,num_cankao_B,num_tongdao_B,num_baohu_B,ts,c,tau_B,D);  

            case '????????????'
                    load data/data_tf
                    fr=1/tr_B;  
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                [noise_tf] =zaoshengtiaofu(fs_tf,Bj_tf,fj_tf,Prj_tf,Tr_tf,frame_tf);
                s_echo_1=s_echo_2+s_noise+noise_tf;%%%%????????
                t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
               [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
               [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
             [s_mtd]=mtd(pc_result1.',M1,num_jilei_B,num_tongdao_B);
                 hengxujing(M1,Pfa_B,s_mtd,num_cankao_B,num_tongdao_B,num_baohu_B,ts,c,tau_B,D);  

            case '????????????'
                load data/data_tx
                fr=1/tr_B;  
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                [noise_tx] =zaoshengtiaoxiang(fs_tx,Bj_tx,fj_tx,Prj_tx,Tr_tx,frame_tx);
                s_echo_1=s_echo_2+s_noise+noise_tx;%%%%????????
                t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
               [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
               [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
              [s_mtd]=mtd(pc_result1.',M1,num_jilei_B,num_tongdao_B);
                 hengxujing(M1,Pfa_B,s_mtd,num_cankao_B,num_tongdao_B,num_baohu_B,ts,c,tau_B,D);
                 
            case '????????????'
               load data/data_smart
                Kfm=4e6;tau=1e-6; c=3e8;
%                 global sigma;
%                 global str;
%                 global RadarS;
%                 global code;
                 Pn=(Bj_smart/(2.5*Kfm))^2;Bn=Bj_smart/2;
                fr=1/tr_B;  
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                [vSmartNoiseSig]=jam_smartnoise( s_echo_2,Pn,Prj_smart,Bn,Kfm,fs_smart );
                s_echo_1=s_echo_2+s_noise+vSmartNoiseSig;%%%%????????
                t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
               [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
               [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_B,num_tongdao_B);
%                 figure,mesh(1:fr/num_tongdao_B:fr,0:ts*(tau_B/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_B/D/ts)-ts*(tau_B/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('????????MTD????');
              hengxujing(M1,Pfa_B,s_mtd,num_cankao_B,num_tongdao_B,num_baohu_B,ts,c,tau_B,D);  
           case '??????????'
               load data/data_shuzhuangpu
                fj=[0.4e6,0.8e6,1.2e6];
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                [noise_szp,~] = jam_combspectrum(fs_shuzhuangpu,Bj_shuzhuangpu,Ns_shuzhuangpu,fj,frame_shuzhuangpu,Prj_shuzhuangpu,Tr_shuzhuangpu);
                s_echo_1=s_echo_2+s_noise+noise_szp;%%%%????????
                t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_B,num_tongdao_B);
                hengxujing(M1,Pfa_B,s_mtd,num_cankao_B,num_tongdao_B,num_baohu_B,ts,c,tau_B,D);  
            case '????????'
                load data/data_saopin
                T_fr=2*Tr_saopin;
                tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                N=length(tm_B);%????????????????????????????
                number1=length(code);
                ts=1/fs_B;
                lamta=c/fz_B;
                fr=1/tr_B;
                An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
%                 Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
%                 Aj=sqrt(Prj);
                %????
                
                [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
%                 [s_ft,echo3]=BKganraoxinhao(D,y,R1,tr_B,ts,c,Aj,N,frame_B,fs_B,f_doppler1,tm_B,f0_B,tau_B,congmubiao,y1);
                [sig_noise,~] = jam_sweepfrequency( fs_saopin,Bj_saopin,fj_saopin,frame_saopin,Prj_saopin,Tr_saopin,T_fr,Time_begin_saopin,K_sweep_saopin );             
                s_echo_1=s_echo_2+s_noise+sig_noise;%%%%????????
                t=0:1/fs_B:frame_B*tr_B-1/fs_B; 
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_B/2,tau_B);
                [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
               [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
               [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
               [s_mtd]=mtd(pc_result1.',M1,num_jilei_B,num_tongdao_B);
              hengxujing(M1,Pfa_B,s_mtd,num_cankao_B,num_tongdao_B,num_baohu_B,ts,c,tau_B,D);  

                end
        
    case '????????????'
        load data/data_JDParameter;
        load data/data_target0Parameter
%         global code;
%         global sigma;
        c=3e8;
%         global strOPSet;
%         global nOPSet;
        switch strOPSet{nOPSet}
           case '????????????'
                load data/data_sp
                lamta=c/fz_JD;%????
                fr=1/tr_JD;
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                [noise_sp] = shepinzaosheng(fs_sp,Bj_sp,fj_sp,frame_sp,Prj_sp,Tr_sp);
                s_echo_1=s_echo_2+s_noise+noise_sp;%%%%????????
                t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD;
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
                hengxujing(M1,Pfa_JD,s_mtd,num_cankao_JD,num_tongdao_JD,num_baohu_JD,ts,c,tau_JD,D);           
                
            case '????????????'
                load data/data_tp
                fr=1/tr_JD;
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                Pn=(Bj_tp/2/(2.5*Kfm))^2;
                [noise_tp]=zaoshengtiaopin(fs_tp,Kfm,Prj_tp,Pn,Bn,fj_tp,frame_tp,Tr_tp);
                s_echo_1=s_echo_2+s_noise+noise_tp;%%%%????????
                t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD;
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
                hengxujing(M1,Pfa_JD,s_mtd,num_cankao_JD,num_tongdao_JD,num_baohu_JD,ts,c,tau_JD,D);

            case '????????????'
                load data/data_tf
                fr=1/tr_JD;
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                [noise_tf] =zaoshengtiaofu(fs_tf,Bj_tf,fj_tf,Prj_tf,Tr_tf,frame_tf);
                s_echo_1=s_echo_2+s_noise+noise_tf;%%%%????????
                t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD;
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
                hengxujing(M1,Pfa_JD,s_mtd,num_cankao_JD,num_tongdao_JD,num_baohu_JD,ts,c,tau_JD,D);

            case '????????????'
                load data/data_tx
                fr=1/tr_JD;
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                [noise_tx] =zaoshengtiaoxiang(fs_tx,Bj_tx,fj_tx,Prj_tx,Tr_tx,frame_tx);
                s_echo_1=s_echo_2+s_noise+noise_tx;%%%%????????
                t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD;
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
                hengxujing(M1,Pfa_JD,s_mtd,num_cankao_JD,num_tongdao_JD,num_baohu_JD,ts,c,tau_JD,D);
               
            case '????????????'
                load data/data_smart
                fr=1/tr_JD;
                Kfm=4e6;
                Pn=(Bj_smart/(2.5*Kfm))^2;Bn=Bj_smart/2;
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                [vSmartNoiseSig]=jam_smartnoise( s_echo_2,Pn,Prj_smart,Bn,Kfm,fs_smart );
                s_echo_1=s_echo_2+s_noise+vSmartNoiseSig;%%%%????????
                t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD;
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
                hengxujing(M1,Pfa_JD,s_mtd,num_cankao_JD,num_tongdao_JD,num_baohu_JD,ts,c,tau_JD,D);
                case '??????????'
                 load data/data_shuzhuangpu
                fj=[0.4e6,0.8e6,1.2e6];
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                [noise_szp,~] = jam_combspectrum(fs_shuzhuangpu,Bj_shuzhuangpu,Ns_shuzhuangpu,fj,frame_shuzhuangpu,Prj_shuzhuangpu,Tr_shuzhuangpu);
                s_echo_1=s_echo_2+s_noise+noise_szp;%%%%????????
                t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD;
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
               [s_mtd]=mtd(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
                hengxujing(M1,Pfa_JD,s_mtd,num_cankao_JD,num_tongdao_JD,num_baohu_JD,ts,c,tau_JD,D);
              
            case '????????'
                load data/data_saopin
                T_fr=2*Tr_saopin;
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                [sig_noise,~] = jam_sweepfrequency( fs_saopin,Bj_saopin,fj_saopin,frame_saopin,Prj_saopin,Tr_saopin,T_fr,Time_begin_saopin,K_sweep_saopin );               
                s_echo_1=s_echo_2+s_noise+sig_noise;%%%%????????
                t=0:1/fs_JD:frame_JD*tr_JD-1/fs_JD;
                s_echo_1=s_echo_1.*rectpuls(t-2*R/c-tau_JD/2,tau_JD);
                [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                [s_mtd]=mtd(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
                hengxujing(M1,Pfa_JD,s_mtd,num_cankao_JD,num_tongdao_JD,num_baohu_JD,ts,c,tau_JD,D);

                
        end
        
    end 
    
elseif n0==4
    str=get(handles.radarSelect,'string');
    n=get(handles.radarSelect,'Value');
    switch  str{n}
        case '????????????'
            c=3e8; 
%             global sigma;
            load data/data_LFMParameter;
            load data/data_target0Parameter

            switch strO{nO}
                case '??????????????'
                    load data/data_dopplerblink
                    R0=rand(100);
                    fr=1/tr_LFM;
                    lamta=c/fz_LFM;%????
                    tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                    N=length(tm);%????????????????????????????
                    An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                    ts=1/fs_LFM;
                    k=B1_LFM/tau_LFM;   
                    Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                    A=sqrt(Prs);%????????????
                    [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                    [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                    f_doppler=2*v/lamta;%????????????????%????????????????????  
                    [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                    [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
    %                 [s_ft,echo3]=LFMganraoxinhao(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao);
                    [ sig_jam,~ ] = jam_dopplerblink( fd_dopplerblink,Td_dopplerblink,R0,s_echo_2,fs_dopplerblink,Pj_dopplerblink,flagT_dopplerblink );
                    s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                    [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                    [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                    [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                    [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                    [s_mtd]=mtd(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
                    hengxujing(M1,Pfa_LFM,s_mtd,num_cankao_LFM,num_tongdao_LFM,num_baohu_LFM,ts,c,tau_LFM,D); 
                case '????????'
                    load data/data_botiao  
                     fr=1/tr_LFM;
                    lamta=c/fz_LFM;%????
                    tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                    N=length(tm);%????????????????????????????
                    An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                    ts=1/fs_LFM;
                    k=B1_LFM/tau_LFM;   
                    Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                    A=sqrt(Prs);%????????????
                    [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                    [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                    f_doppler=2*v/lamta;%????????????????%????????????????????  
                    [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                    [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
    %                 [s_ft,echo3]=LFMganraoxinhao(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao);
                        %????????????
                    px=1e3;py=1e3;pz=1e3;%????????
                    vx=10;vy=10;vz=0;%????????
                    ax=0;ay=0;az=0;%??????????
                    phi=pi/180;
                    [ TargetStatus ] = paraset_targetstatus( px,py,pz,vx,vy,vz,ax,ay,az,phi );
                    [ PassivePara ] = paraset_passivejaming( tf_botiao,sf_botiao,vl_botiao,vf_botiao,ts_botiao,bt_botiao,al_botiao,sref_botiao,smax_botiao );
                      vx=1;vy=1;
                    [ WindV ] = paraset_windvelocity( vx,vy );
                    CurrentT=2.2;
                    [ sig_jam,~ ] = jam_passive( s_echo_2,fs_LFM,f0_LFM,CurrentT,TargetStatus,WindV,PassivePara );
                    s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                    [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                    [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                    [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                    [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                     [s_mtd]=mtd(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
                    hengxujing(M1,Pfa_LFM,s_mtd,num_cankao_LFM,num_tongdao_LFM,num_baohu_LFM,ts,c,tau_LFM,D);
                    
                case 'AGC????' 
                    load data/data_AGC 
                     fr=1/tr_LFM;
                    lamta=c/fz_LFM;%????
                    tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                    N=length(tm);%????????????????????????????
                    An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
                    ts=1/fs_LFM;
                    k=B1_LFM/tau_LFM;   
                    Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
                    A=sqrt(Prs);%????????????
                    [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                    [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                    f_doppler=2*v/lamta;%????????????????%????????????????????  
                    [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                    [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
    %                 [s_ft,echo3]=LFMganraoxinhao(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao);
                    [ sig_jam,~ ] = jam_AGC( CurrentT_AGC,Pj_AGC,Period_AGC,D_AGC,radio_AGC,s_echo_2,fs_AGC );
                    s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                    [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
                    [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
                    [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                    [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                     [s_mtd]=mtd(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
                    hengxujing(M1,Pfa_LFM,s_mtd,num_cankao_LFM,num_tongdao_LFM,num_baohu_LFM,ts,c,tau_LFM,D);
                    
            end
        case '????????????'
            load data/data_BFParameter;
            load data/data_target0Parameter

%             global code;
%             global sigma;
            c=3e8; 
%             global strO;
%             global nO;
            switch strO{nO}
                case '??????????????'
                    load data/data_dopplerblink
                    R0=rand(100);
                    fj=[0.4e6,0.8e6,1.2e6];
                    tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                    N=length(tm_B);%????????????????????????????
                    number1=length(code);
                    ts=1/fs_B;
                    lamta=c/fz_B;
                    fr=1/tr_B;
                    An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                    Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                    A=sqrt(Prs);%????????????
                    f_doppler=2*v/lamta;%????????????????
                    [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                    [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                    [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                    [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                    [ sig_jam,~ ] = jam_dopplerblink( fd_dopplerblink,Td_dopplerblink,R0,s_echo_2,fs_dopplerblink,Pj_dopplerblink,flagT_dopplerblink );
                    s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                    [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                    [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
                    [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
                    [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
                    [s_mtd]=mtd(pc_result1.',M1,num_jilei_B,num_tongdao_B);
                    hengxujing(M1,Pfa_B,s_mtd,num_cankao_B,num_tongdao_B,num_baohu_B,ts,c,tau_B,D);
                    
                case '????????'
                    load data/data_botiao  
                    fj=[0.4e6,0.8e6,1.2e6];
                    tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                    N=length(tm_B);%????????????????????????????
                    number1=length(code);
                    ts=1/fs_B;
                    lamta=c/fz_B;
                    fr=1/tr_B;
                    An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                    Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                    A=sqrt(Prs);%????????????
                    f_doppler=2*v/lamta;%????????????????
                    [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                    [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                    [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                    [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                        %????????????
                    px=1e3;py=1e3;pz=1e3;%????????
                    vx=10;vy=10;vz=0;%????????
                    ax=0;ay=0;az=0;%??????????
                    phi=pi/180;
                    [ TargetStatus ] = paraset_targetstatus( px,py,pz,vx,vy,vz,ax,ay,az,phi );
                    [ PassivePara ] = paraset_passivejaming( tf_botiao,sf_botiao,vl_botiao,vf_botiao,ts_botiao,bt_botiao,al_botiao,sref_botiao,smax_botiao );
                      vx=1;vy=1;
                    [ WindV ] = paraset_windvelocity( vx,vy );
                    CurrentT=2.2;
                    [ sig_jam,~ ] = jam_passive( s_echo_2,fs_B,f0_B,CurrentT,TargetStatus,WindV,PassivePara );
                    s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                    [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                    [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
                    [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
                    [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
                    [s_mtd]=mtd(pc_result1.',M1,num_jilei_B,num_tongdao_B);
                    hengxujing(M1,Pfa_B,s_mtd,num_cankao_B,num_tongdao_B,num_baohu_B,ts,c,tau_B,D);
                    
                case 'AGC????' 
                    load data/data_AGC 
                    fj=[0.4e6,0.8e6,1.2e6];
                    tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
                    N=length(tm_B);%????????????????????????????
                    number1=length(code);
                    ts=1/fs_B;
                    lamta=c/fz_B;
                    fr=1/tr_B;
                    An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
                    Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
                    A=sqrt(Prs);%????????????
                    f_doppler=2*v/lamta;%????????????????
                    [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
                    [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
                    [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
                    [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
                    [ sig_jam,~ ] = jam_AGC( CurrentT_AGC,Pj_AGC,Period_AGC,D_AGC,radio_AGC,s_echo_2,fs_AGC );
                    s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                    [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
                    [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
                    [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
                    [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
                    [s_mtd]=mtd(pc_result1.',M1,num_jilei_B,num_tongdao_B);
                    hengxujing(M1,Pfa_B,s_mtd,num_cankao_B,num_tongdao_B,num_baohu_B,ts,c,tau_B,D);
                    
            end
            
        case '????????????'
            load data/data_JDParameter;
            load data/data_target0Parameter
%             global code;
%             global sigma;
            c=3e8;
%             global strO;
%             global nO;
            switch strO{nO}
                case '??????????????'
                    load data/data_dopplerblink
                    R0=rand(100);
                    fj=[0.4e6,0.8e6,1.2e6];
                    fr=1/tr_JD;
                    lamta=c/fz_JD;%????
                    tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                    N=length(tm);%????????????????????????????
                    An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                    [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                    ts=1/fs_JD; 
                    Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                    A=sqrt(Prs);%????????????
                    f_doppler=2*v/lamta;%????????????????
                    [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                    [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                    [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                    [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                    [ sig_jam,~ ] = jam_dopplerblink( fd_dopplerblink,Td_dopplerblink,R0,s_echo_2,fs_dopplerblink,Pj_dopplerblink,flagT_dopplerblink );
                    s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                    [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                    [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                    [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                    [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                    [s_mtd]=mtd(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
                    hengxujing(M1,Pfa_JD,s_mtd,num_cankao_JD,num_tongdao_JD,num_baohu_JD,ts,c,tau_JD,D);

                case '????????'
                    load data/data_botiao  
                    fj=[0.4e6,0.8e6,1.2e6];
                    fr=1/tr_JD;
                    lamta=c/fz_JD;%????
                    tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                    N=length(tm);%????????????????????????????
                    An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                    [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                    ts=1/fs_JD; 
                    Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                    A=sqrt(Prs);%????????????
                    f_doppler=2*v/lamta;%????????????????
                    [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                    [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                    [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                    [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                        %????????????
                    px=1e3;py=1e3;pz=1e3;%????????
                    vx=10;vy=10;vz=0;%????????
                    ax=0;ay=0;az=0;%??????????
                    phi=pi/180;
                    [ TargetStatus ] = paraset_targetstatus( px,py,pz,vx,vy,vz,ax,ay,az,phi );
                    [ PassivePara ] = paraset_passivejaming( tf_botiao,sf_botiao,vl_botiao,vf_botiao,ts_botiao,bt_botiao,al_botiao,sref_botiao,smax_botiao );
                      vx=1;vy=1;
                    [ WindV ] = paraset_windvelocity( vx,vy );
                    CurrentT=2.2;
                    [ sig_jam,~ ] = jam_passive( s_echo_2,fs_JD,f0_JD,CurrentT,TargetStatus,WindV,PassivePara );
                    s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                    [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                    [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                    [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                    [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                    [s_mtd]=mtd(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
                    hengxujing(M1,Pfa_JD,s_mtd,num_cankao_JD,num_tongdao_JD,num_baohu_JD,ts,c,tau_JD,D);

                case 'AGC????' 
                    load data/data_AGC 
                    fj=[0.4e6,0.8e6,1.2e6];
                    fr=1/tr_JD;
                    lamta=c/fz_JD;%????
                    tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                    N=length(tm);%????????????????????????????
                    An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                    [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                    ts=1/fs_JD; 
                    Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                    A=sqrt(Prs);%????????????
                    f_doppler=2*v/lamta;%????????????????
                    [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                    [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
                    [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                    [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                    [ sig_jam,~ ] = jam_AGC( CurrentT_AGC,Pj_AGC,Period_AGC,D_AGC,radio_AGC,s_echo_2,fs_AGC );
                    s_echo_1=s_echo_2+s_noise+sig_jam;%%%%????????
                    [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
                    [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
                    [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                    [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                    [s_mtd]=mtd(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
                    hengxujing(M1,Pfa_JD,s_mtd,num_cankao_JD,num_tongdao_JD,num_baohu_JD,ts,c,tau_JD,D);

    
            end
    end
end
% --- Executes during object creation, after setting all properties.

% --- Executes during object creation, after setting all properties.
function uipanel4_CreateFcn(~, ~, ~)
% hObject    handle to uipanel4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in Jzongtifangzhen.
function Jzongtifangzhen_Callback(~, ~, handles)%??????????????
% hObject    handle to Jzongtifangzhen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global RadarS;

c=3e8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%????????

global rcsk;
% load data/data_LFMParameter
 load data/data_target0Parameter
% load data/data_DeceptionJammingParameter
%    LFMmainJ(fz_LFM,B_LFM,B1_LFM,Gt_LFM,Gr_LFM,F_LFM,Pt_LFM,L_LFM,Te_LFM,tau_LFM,tr_LFM,f0_LFM,f1_LFM,fs_LFM,frame_LFM,num_jilei_LFM,num_tongdao_LFM,num_cankao_LFM,num_baohu_LFM,Pfa_LFM,R,v,rcsk,sigma0) ;
% load data/data_target0Parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global temp1;
    global string2;
    global temp2;
    global strOPSet;
    global nOPSet;
    global strO;
    global nO;
if RadarS==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LFM????????       
    str0=get(handles.jammingSelect,'string');
    n0=get(handles.jammingSelect,'Value');
    if n0==3
         load data/data_LFMParameter
           load data/data_DeceptionJammingParameter;
     
           load data/data_target0Parameter;
            switch string2{temp2}
            case {'??' ,'??','????' }
               LFMmainJ(fz_LFM,B_LFM,B1_LFM,Gt_LFM,Gr_LFM,F_LFM,Pt_LFM,L_LFM,Te_LFM,tau_LFM,tr_LFM,f0_LFM,f1_LFM,fs_LFM,frame_LFM,num_jilei_LFM,num_tongdao_LFM,num_cankao_LFM,num_baohu_LFM,Pfa_LFM,R,v,rcsk,sigma0,temp1) ;
            case {'????'}
                  global sigma;             
               fr=1/tr_LFM;%??????????
                f1=10e6; %????????
                c=3e8;
                lamta=c/fz_LFM;
                Rj=2e3;
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L_LFM/10)));
                Aj=sqrt(Prj);
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
                N=length(tm);%????????????????????????????
                f_doppler=2*v/lamta;
                f_doppler1=f_doppler;
                k=B1_LFM/tau_LFM;%????????????????????
                ts=1/fs_LFM;%????????
                %????????????
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L/10))); %????????????????
                A=sqrt(Prs);
                Te=1143;%????
                F=4;%????????
                 tf=10;
                An=10*log10((1.382e-23)*Te*B_LFM*10^(F/10));%????????
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k); 
                figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(????????)'), ylabel('y(????????)'),title('????????');
                figure,plot((0:fs_LFM/length(abs(fft(abs(fftshift(y))))):fs_LFM-fs_LFM/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('????f(??????Hz)'), ylabel('y(????????/Hz)'),title('??????????????');
                [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
                %????????????
               [s_echo_2,~]=LFMtuoyinhuiboxinhao(y,R,tr_LFM,ts,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
              %figure,plot(0:ts:(N-1)*ts,real(s_echo_2(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????');
                %figure,plot((0:fs/length(abs(fft(abs(fftshift(echo(1,:)))))):fs-fs/length(abs(fft(abs(fftshift(echo(1,:))))))),abs((fft(echo(1,:))))),xlabel('????f(??????Hz)'), ylabel('y(????????)'),title('??????????????');

                %????????????
                [s_ft,echo3]=LFMtuoyinganrao(R,frame_LFM,tf,Aj,tm,tau_LFM,f0_LFM,B1_LFM,N,f_doppler1,f_doppler,vf,lamta,k,ts,c,fs_LFM,temp1);

                %????????
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                %??????????????????????????????????????????????????
                s_echo_1=s_echo_2+s_noise+s_ft;
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????????');
                figure,plot((0:fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_LFM-fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('????f(??????Hz)'), ylabel('y(????????/Hz)'),title('??????????????????');
                
                 %????
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1,fs_LFM,f0_LFM);
                %????????????
                [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
                %????????????????
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);
                figure,plot(0:ts*c/2:(M-1)*ts*c/2,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????');

%                 s_pc_result=reshape(pc_result,1,M1*frame_LFM);
                %??????????
                [s_mtd,s_mtd1]=MTD_tuoyin(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
                figure,subplot(121),mesh(1:fr/num_tongdao_LFM*lamta/2:fr*lamta/2,0:ts*(tau_LFM/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_LFM/D/ts)-ts*(tau_LFM/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????/??'),ylabel('??????????????'),zlabel('y(????????)'),title('????????1??16??????????MTD????');
                subplot(122),mesh(1:fr/num_tongdao_LFM*lamta/2:fr*lamta/2,0:ts*(tau_LFM/D/ts)*c/2:(length(abs(s_mtd1(:,1)))*ts*(tau_LFM/D/ts)-ts*(tau_LFM/D/ts))*c/2,abs(s_mtd1)),xlabel('??????????????/??'),ylabel('??????????????'),zlabel('y(????????)'),title('????????113??128??????????MTD????');

                % figure,mesh(1:fr/num_tongdao*lamta/2:fr*lamta/2,0:ts*(tau/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau/D/ts)-ts*(tau/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('MTD????');
                hengxujing_tuoyin(M1,Pfa_LFM,s_mtd,num_cankao_LFM,num_tongdao_LFM,num_baohu_LFM,ts,c,tau_LFM,D,s_mtd1);

            end
    elseif n0==2

                switch strOPSet{nOPSet}
                case '????????????'
                   load data/data_LFMParameter
                   load data/data_target0Parameter;
                   LFMmainJO_sp(fz_LFM,B_LFM,B1_LFM,Gt_LFM,Gr_LFM,F_LFM,Pt_LFM,L_LFM,Te_LFM,tau_LFM,tr_LFM,f0_LFM,f1_LFM,fs_LFM,frame_LFM,num_jilei_LFM,num_tongdao_LFM,num_cankao_LFM,num_baohu_LFM,Pfa_LFM,R,v,rcsk,sigma0) ;
                 case '????????????'
                   load data/data_LFMParameter
                   load data/data_target0Parameter;
                   LFMmainJO_tx(fz_LFM,B_LFM,B1_LFM,Gt_LFM,Gr_LFM,F_LFM,Pt_LFM,L_LFM,Te_LFM,tau_LFM,tr_LFM,f0_LFM,f1_LFM,fs_LFM,frame_LFM,num_jilei_LFM,num_tongdao_LFM,num_cankao_LFM,num_baohu_LFM,Pfa_LFM,R,v,rcsk,sigma0) ;
                case '????????????'
                   load data/data_LFMParameter
                   load data/data_target0Parameter;
                   LFMmainJO_tf(fz_LFM,B_LFM,B1_LFM,Gt_LFM,Gr_LFM,F_LFM,Pt_LFM,L_LFM,Te_LFM,tau_LFM,tr_LFM,f0_LFM,f1_LFM,fs_LFM,frame_LFM,num_jilei_LFM,num_tongdao_LFM,num_cankao_LFM,num_baohu_LFM,Pfa_LFM,R,v,rcsk,sigma0) ;                       
                case '????????????'
                   load data/data_LFMParameter
                   load data/data_target0Parameter;
                   LFMmainJO_tp(fz_LFM,B_LFM,B1_LFM,Gt_LFM,Gr_LFM,F_LFM,Pt_LFM,L_LFM,Te_LFM,tau_LFM,tr_LFM,f0_LFM,f1_LFM,fs_LFM,frame_LFM,num_jilei_LFM,num_tongdao_LFM,num_cankao_LFM,num_baohu_LFM,Pfa_LFM,R,v,rcsk,sigma0) ;
                case '????????????'
                   load data/data_LFMParameter
                   load data/data_target0Parameter;
                   LFMmainJO_smart(fz_LFM,B_LFM,B1_LFM,Gt_LFM,Gr_LFM,F_LFM,Pt_LFM,L_LFM,Te_LFM,tau_LFM,tr_LFM,f0_LFM,f1_LFM,fs_LFM,frame_LFM,num_jilei_LFM,num_tongdao_LFM,num_cankao_LFM,num_baohu_LFM,Pfa_LFM,R,v,rcsk,sigma0) ;
                case '??????????'
                   load data/data_LFMParameter
                   load data/data_target0Parameter;
                   LFMmainJO_shuzhuangpu(fz_LFM,B_LFM,B1_LFM,Gt_LFM,Gr_LFM,F_LFM,Pt_LFM,L_LFM,Te_LFM,tau_LFM,tr_LFM,f0_LFM,f1_LFM,fs_LFM,frame_LFM,num_jilei_LFM,num_tongdao_LFM,num_cankao_LFM,num_baohu_LFM,Pfa_LFM,R,v,rcsk,sigma0) ;
                case '????????'
                   load data/data_LFMParameter
                   load data/data_target0Parameter;
                   LFMmainJO_saopin(fz_LFM,B_LFM,B1_LFM,Gt_LFM,Gr_LFM,F_LFM,Pt_LFM,L_LFM,Te_LFM,tau_LFM,tr_LFM,f0_LFM,f1_LFM,fs_LFM,frame_LFM,num_jilei_LFM,num_tongdao_LFM,num_cankao_LFM,num_baohu_LFM,Pfa_LFM,R,v,rcsk,sigma0) ;

                end
        
        elseif n0==4
            switch strO{nO}
            case '??????????????'
                load data/data_LFMParameter;
                       load data/data_target0Parameter;
                       LFMmainJ_dopplerblink(fz_LFM,B_LFM,B1_LFM,Gt_LFM,Gr_LFM,F_LFM,Pt_LFM,L_LFM,Te_LFM,tau_LFM,tr_LFM,f0_LFM,f1_LFM,fs_LFM,frame_LFM,num_jilei_LFM,num_tongdao_LFM,num_cankao_LFM,num_baohu_LFM,Pfa_LFM,R,v,rcsk,sigma0);
            
            case '????????'
                load data/data_LFMParameter;
                       load data/data_target0Parameter;
                       LFMmainJ_passive(fz_LFM,B_LFM,B1_LFM,Gt_LFM,Gr_LFM,F_LFM,Pt_LFM,L_LFM,Te_LFM,tau_LFM,tr_LFM,f0_LFM,f1_LFM,fs_LFM,frame_LFM,num_jilei_LFM,num_tongdao_LFM,num_cankao_LFM,num_baohu_LFM,Pfa_LFM,R,v,rcsk,sigma0);
            case'AGC????'
              load data/data_LFMParameter;
                       load data/data_target0Parameter;
                       LFMmainJ_AGC(fz_LFM,B_LFM,B1_LFM,Gt_LFM,Gr_LFM,F_LFM,Pt_LFM,L_LFM,Te_LFM,tau_LFM,tr_LFM,f0_LFM,f1_LFM,fs_LFM,frame_LFM,num_jilei_LFM,num_tongdao_LFM,num_cankao_LFM,num_baohu_LFM,Pfa_LFM,R,v,rcsk,sigma0);   
           end
    end

elseif RadarS==2
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BF????????
    global code; 
%     global temp1;
%     global string2;
%     global temp2;
    str0=get(handles.jammingSelect,'string');
    n0=get(handles.jammingSelect,'Value');
    if n0==3
%         global temp1;
%         global string2;
%         global temp2;
        
         load data/data_BFParameter
         load data/data_DeceptionJammingParameter;
         load data/data_target0Parameter;
         global number1;
        number1=length(code);
           switch string2{temp2}
            case {'??' ,'??','????' }
                BKmainJ(fz_B,B_B,Gt_B,Gr_B,F_B,Pt_B,L_B,Pfa_B,R,v,rcsk,sigma0,Te_B,tau_B,tr_B,f0_B,f1_B,fs_B,frame_B,num_jilei_B,num_tongdao_B,num_cankao_B,num_baohu_B,flag,number1,code,temp1)  ;
            case {'????'}
               fr=1/tr_B;
               Rj=2e3;
               f1=10e6;    
            number1=length(code);
            ts=1/fs_B;
            lamta=c/fz_B;
            tm=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
            N=length(tm);%????????????????????????????
            f_doppler=2*v/lamta;%????????????????
            f_doppler1= f_doppler;%????????????????
            Te=1143;%????
            An=10*log10((1.382e-23)*Te*B_B*10^(F_B/10));%????????
            Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*0.5)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
            A=sqrt(Prs);%????????????
            Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L/10)));
            Aj=sqrt(Prj);
            [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
             figure,plot((0:fs_B/length(abs(fft(abs(fftshift(y))))):fs_B-fs_B/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('????f(??????Hz)'), ylabel('y(????????/Hz)'),title('??????????????');
%             %%%%%%%%%%%1.1.2????????????????????
            [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
             %%%%%%%%%%%%%1.3????????????%%%%%%%%%%%%
            [s_echo_2,echo]=BKtuoyinhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tm,tau_B); 
            figure,plot(0:ts:(N-1)*ts,real(s_echo_2(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????');
            figure,plot((0:fs_B/length(abs(fft(abs(fftshift(echo(1,:)))))):fs_B-fs_B/length(abs(fft(abs(fftshift(echo(1,:))))))),abs((fft(echo(1,:))))),xlabel('????f(??????Hz)'), ylabel('y(????????/Hz)'),title('??????????????');

            %????????
            [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
            %????????
            [s_ft,~]=BKtuoyinganrao(R,frame_B,tf,Aj,tau_B,N,f_doppler1,f_doppler,vf,lamta,ts,c,fs_B,y1,tr_B,temp1);
            %??????????????????????????????????????????????????
            s_echo_1=s_echo_2+s_noise+s_ft;
            %????
            [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1,fs_B,f0_B);
            %????????????
            [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
            %????????????????
            [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
            figure,plot(0:ts*c/2:(M-1)*ts*c/2,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????');
            s_pc_result=reshape(pc_result,1,M1*frame_B);
            %??????????
            [s_mtd,s_mtd1]=MTD_tuoyin(pc_result1.',M1,num_jilei_B,num_tongdao_B);
            subplot(121),mesh(1:fr/num_tongdao_B*lamta/2:fr*lamta/2,0:ts*(tau_B/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_B/D/ts)-ts*(tau_B/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????/??'),ylabel('??????????????'),zlabel('y(????????)'),title('????????1??16??????????MTD????');
            subplot(122),mesh(1:fr/num_tongdao_B*lamta/2:fr*lamta/2,0:ts*(tau_B/D/ts)*c/2:(length(abs(s_mtd1(:,1)))*ts*(tau_B/D/ts)-ts*(tau_B/D/ts))*c/2,abs(s_mtd1)),xlabel('??????????????/??'),ylabel('??????????????'),zlabel('y(????????)'),title('????????113??128??????????MTD????');

            hengxujing_tuoyin(M1,Pfa_B,s_mtd,num_cankao_B,num_tongdao_B,num_baohu_B,ts,c,tau_B,D,s_mtd1);

            end
    elseif n0==2
%                     global strOPSet;
%                     global nOPSet;
            number1=length(code);

            switch strOPSet{nOPSet}
            case '????????????'
               load data/data_BFParameter
               load data/data_target0Parameter;
               BKmainJO_sp(fz_B,B_B,Gt_B,Gr_B,F_B,Pt_B,L_B,Pfa_B,R,v,rcsk,sigma0,Te_B,tau_B,tr_B,f0_B,f1_B,fs_B,frame_B,num_jilei_B,num_tongdao_B,num_cankao_B,num_baohu_B,flag,number1,code)  ;
             case '????????????'
               load data/data_BFParameter
               load data/data_target0Parameter;
               BKmainJO_tx(fz_B,B_B,Gt_B,Gr_B,F_B,Pt_B,L_B,Pfa_B,R,v,rcsk,sigma0,Te_B,tau_B,tr_B,f0_B,f1_B,fs_B,frame_B,num_jilei_B,num_tongdao_B,num_cankao_B,num_baohu_B,flag,number1,code)  ;
            case '????????????'
               load data/data_BFParameter
               load data/data_target0Parameter;
               BKmainJO_tf(fz_B,B_B,Gt_B,Gr_B,F_B,Pt_B,L_B,Pfa_B,R,v,rcsk,sigma0,Te_B,tau_B,tr_B,f0_B,f1_B,fs_B,frame_B,num_jilei_B,num_tongdao_B,num_cankao_B,num_baohu_B,flag,number1,code)  ;
            case '????????????'
               load data/data_BFParameter
               load data/data_target0Parameter;
               BKmainJO_tp(fz_B,B_B,Gt_B,Gr_B,F_B,Pt_B,L_B,Pfa_B,R,v,rcsk,sigma0,Te_B,tau_B,tr_B,f0_B,f1_B,fs_B,frame_B,num_jilei_B,num_tongdao_B,num_cankao_B,num_baohu_B,flag,number1,code)  ;
           case '????????????'
               load data/data_BFParameter
               load data/data_target0Parameter;
               BKmainJO_smart(fz_B,B_B,Gt_B,Gr_B,F_B,Pt_B,L_B,Pfa_B,R,v,rcsk,sigma0,Te_B,tau_B,tr_B,f0_B,f1_B,fs_B,frame_B,num_jilei_B,num_tongdao_B,num_cankao_B,num_baohu_B,flag,number1,code)  ;
           case '??????????'
               load data/data_BFParameter
               load data/data_target0Parameter;
               BKmainJO_shuzhuangpu(fz_B,B_B,Gt_B,Gr_B,F_B,Pt_B,L_B,Pfa_B,R,v,rcsk,sigma0,Te_B,tau_B,tr_B,f0_B,f1_B,fs_B,frame_B,num_jilei_B,num_tongdao_B,num_cankao_B,num_baohu_B,flag,number1,code)  ;
           case '????????'
               load data/data_BFParameter
               load data/data_target0Parameter;
               BKmainJO_saopin(fz_B,B_B,Gt_B,Gr_B,F_B,Pt_B,L_B,Pfa_B,R,v,rcsk,sigma0,Te_B,tau_B,tr_B,f0_B,f1_B,fs_B,frame_B,num_jilei_B,num_tongdao_B,num_cankao_B,num_baohu_B,flag,number1,code)  ;


            end
                               
    elseif n0==4
        
        switch strO{nO}
            case '??????????????'
                load data/data_BFParameter
         load data/data_dopplerblink;
         load data/data_target0Parameter;
%           global number1;
         number1=length(code);
                BKmainJ_dopplerblink(fz_B,B_B,Gt_B,Gr_B,F_B,Pt_B,L_B,Pfa_B,R,v,rcsk,sigma0,Te_B,tau_B,tr_B,f0_B,f1_B,fs_B,frame_B,num_jilei_B,num_tongdao_B,num_cankao_B,num_baohu_B,flag,number1,code)  ;
            case '????????'
                
                load data/data_BFParameter
         load data/data_botiao;
         load data/data_target0Parameter;
%           global number1;
         number1=length(code);
                BKmainJ_passive(fz_B,B_B,Gt_B,Gr_B,F_B,Pt_B,L_B,Pfa_B,R,v,rcsk,sigma0,Te_B,tau_B,tr_B,f0_B,f1_B,fs_B,frame_B,num_jilei_B,num_tongdao_B,num_cankao_B,num_baohu_B,flag,number1,code)  ;
            case'AGC????'
                 load data/data_BFParameter
         load data/data_AGC;
         load data/data_target0Parameter;
%           global number1;
          number1=length(code);
                BKmainJ_AGC(fz_B,B_B,Gt_B,Gr_B,F_B,Pt_B,L_B,Pfa_B,R,v,rcsk,sigma0,Te_B,tau_B,tr_B,f0_B,f1_B,fs_B,frame_B,num_jilei_B,num_tongdao_B,num_cankao_B,num_baohu_B,flag,number1,code)  ;

        end
    end  
elseif RadarS==3
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%JD???????? 
            
          str0=get(handles.jammingSelect,'string');
          n0=get(handles.jammingSelect,'Value');
          if n0==3
               load data/data_JDParameter
             load data/data_DeceptionJammingParameter;
             load data/data_target0Parameter;
             switch string2{temp2}
                 case {'??' ,'??','????' }
                 JDmainJ(fz_JD,B_JD,Gt_JD,Gr_JD,F_JD,Pt_JD,L_JD,Pfa_JD,R,v,rcsk,sigma0,Te_JD,tau_JD,tr_JD,f0_JD,f1_JD,fs_JD,frame_JD,num_jilei_JD,num_tongdao_JD,num_cankao_JD,num_baohu_JD,temp1) ;
                 case {'????'}
                Rj=2e3;
                c=3e8;
                f1=10e6;
                fr=1/tr_JD;
                lamta=c/fz_JD;%????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????
                An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                ts=1/fs_JD; 
                Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
                A=sqrt(Prs);%????????????
                f_doppler=2*v/lamta;%????????????????

                %????
                f_doppler=2*v/lamta;
                f_doppler1=f_doppler;
                lamta=c/fz_JD;
                f_doppler=2*v/lamta;
                [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
                figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(????????)'), ylabel('y(????????)'),title('????????');
                figure,plot((0:fs_JD/length(abs(fft(abs(fftshift(y))))):fs_JD-fs_JD/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('????f(??????Hz)'), ylabel('y(????????/Hz)'),title('??????????????');

                ts=1/fs_JD;%????????
                tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
                N=length(tm);%????????????????????????????

                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L/10)));
                Aj=sqrt(Prj);
        %         figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(????????)'), ylabel('y(????????)'),title('????????');
        %         figure,plot((0:fs/length(abs(fft(abs(fftshift(y))))):fs-fs/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('????f(??????Hz)'), ylabel('y(????????)'),title('??????????????');
                %????????????????????
                [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);

                %????????????
                [s_echo_2,echo]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
                figure,plot(0:ts:(N-1)*ts,real(s_echo_2(1:N))),xlabel('t(??????s)'), ylabel('y(????????)'),title('????????');
                figure,plot((0:fs_JD/length(abs(fft(abs(fftshift(echo(1,:)))))):fs_JD-fs_JD/length(abs(fft(abs(fftshift(echo(1,:))))))),abs((fft(echo(1,:))))),xlabel('????f(??????Hz)'), ylabel('y(????????/Hz)'),title('??????????????');

                [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
                %????????????
                [s_ft,~]=JDtuoyinganrao(R,frame_JD,tf,Aj,tm,tau_JD,f0_JD,N,f_doppler1,f_doppler,vf,lamta,ts,c,fs_JD,temp1);

                %??????????????????????????????????????????????????
                s_echo_1=s_echo_2+s_noise+s_ft;
                [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1,fs_JD,f0_JD);
                %????????????
                [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
                %????????????????
                [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
                figure,plot(0:ts*c/2:(M-1)*ts*c/2,20*log10(abs(pc_result1(1,:)))),xlabel('t(??????s)'), ylabel('y(??????dB)'),title('????????');
                s_pc_result=reshape(pc_result,1,M1*frame_JD);
                %??????????
                [s_mtd,s_mtd1]=MTD_tuoyin(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
                figure,subplot(121),mesh(1:fr/num_tongdao_JD*lamta/2:fr*lamta/2,0:ts*(tau_JD/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau_JD/D/ts)-ts*(tau_JD/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????/??'),ylabel('??????????????'),zlabel('y(????????)'),title('????????1??16??????????MTD????');
                subplot(122),mesh(1:fr/num_tongdao_JD*lamta/2:fr*lamta/2,0:ts*(tau_JD/D/ts)*c/2:(length(abs(s_mtd1(:,1)))*ts*(tau_JD/D/ts)-ts*(tau_JD/D/ts))*c/2,abs(s_mtd1)),xlabel('??????????????/??'),ylabel('??????????????'),zlabel('y(????????)'),title('????????113??128??????????MTD????');

                % figure,mesh(1:fr/num_tongdao*lamta/2:fr*lamta/2,0:ts*(tau/D/ts)*c/2:(length(abs(s_mtd(:,1)))*ts*(tau/D/ts)-ts*(tau/D/ts))*c/2,abs(s_mtd)),xlabel('??????????????????Hz'),ylabel('??????????????'),zlabel('y(????????)'),title('MTD????');
                hengxujing_tuoyin(M1,Pfa_JD,s_mtd,num_cankao_JD,num_tongdao_JD,num_baohu_JD,ts,c,tau_JD,D,s_mtd1);

            end
          elseif n0==2
                    switch strOPSet{nOPSet}
                    case '????????????'
                       load data/data_JDParameter
                       load data/data_target0Parameter;
                 JDmainJ_sp(fz_JD,B_JD,Gt_JD,Gr_JD,F_JD,Pt_JD,L_JD,Pfa_JD,R,v,rcsk,sigma0,Te_JD,tau_JD,tr_JD,f0_JD,f1_JD,fs_JD,frame_JD,num_jilei_JD,num_tongdao_JD,num_cankao_JD,num_baohu_JD) ;
                     case '????????????'
                       load data/data_JDParameter
                       load data/data_target0Parameter;
                 JDmainJ_tx(fz_JD,B_JD,Gt_JD,Gr_JD,F_JD,Pt_JD,L_JD,Pfa_JD,R,v,rcsk,sigma0,Te_JD,tau_JD,tr_JD,f0_JD,f1_JD,fs_JD,frame_JD,num_jilei_JD,num_tongdao_JD,num_cankao_JD,num_baohu_JD) ;
                    case '????????????'
                       load data/data_JDParameter
                       load data/data_target0Parameter;
                 JDmainJ_tf(fz_JD,B_JD,Gt_JD,Gr_JD,F_JD,Pt_JD,L_JD,Pfa_JD,R,v,rcsk,sigma0,Te_JD,tau_JD,tr_JD,f0_JD,f1_JD,fs_JD,frame_JD,num_jilei_JD,num_tongdao_JD,num_cankao_JD,num_baohu_JD) ;
                    case '????????????'
                       load data/data_JDParameter
                       load data/data_target0Parameter;
                 JDmainJ_tp(fz_JD,B_JD,Gt_JD,Gr_JD,F_JD,Pt_JD,L_JD,Pfa_JD,R,v,rcsk,sigma0,Te_JD,tau_JD,tr_JD,f0_JD,f1_JD,fs_JD,frame_JD,num_jilei_JD,num_tongdao_JD,num_cankao_JD,num_baohu_JD) ;
                    case '????????????'
                       load data/data_JDParameter
                       load data/data_target0Parameter;
                JDmainJ_smart(fz_JD,B_JD,Gt_JD,Gr_JD,F_JD,Pt_JD,L_JD,Pfa_JD,R,v,rcsk,sigma0,Te_JD,tau_JD,tr_JD,f0_JD,f1_JD,fs_JD,frame_JD,num_jilei_JD,num_tongdao_JD,num_cankao_JD,num_baohu_JD) ;
                   case '??????????'
                       load data/data_JDParameter
                       load data/data_target0Parameter;
                JDmainJ_shuzhuangpu(fz_JD,B_JD,Gt_JD,Gr_JD,F_JD,Pt_JD,L_JD,Pfa_JD,R,v,rcsk,sigma0,Te_JD,tau_JD,tr_JD,f0_JD,f1_JD,fs_JD,frame_JD,num_jilei_JD,num_tongdao_JD,num_cankao_JD,num_baohu_JD) ;
                    case '????????'
                       load data/data_JDParameter
                       load data/data_target0Parameter;
                JDmainJ_saopin(fz_JD,B_JD,Gt_JD,Gr_JD,F_JD,Pt_JD,L_JD,Pfa_JD,R,v,rcsk,sigma0,Te_JD,tau_JD,tr_JD,f0_JD,f1_JD,fs_JD,frame_JD,num_jilei_JD,num_tongdao_JD,num_cankao_JD,num_baohu_JD) ;

                    end
          
          elseif n0==4
        
                switch strO{nO}
                    case '??????????????'
                 load data/data_JDParameter;
                 load data/data_dopplerblink;
                 load data/data_target0Parameter;

                          JDmainJ_dopplerblink(fz_JD,B_JD,Gt_JD,Gr_JD,F_JD,Pt_JD,L_JD,Pfa_JD,R,v,rcsk,sigma0,Te_JD,tau_JD,tr_JD,f0_JD,f1_JD,fs_JD,frame_JD,num_jilei_JD,num_tongdao_JD,num_cankao_JD,num_baohu_JD) ;

                     case '????????'
                 load data/data_JDParameter;
                 load data/data_botiao;
                 load data/data_target0Parameter;
                          JDmainJ_passive(fz_JD,B_JD,Gt_JD,Gr_JD,F_JD,Pt_JD,L_JD,Pfa_JD,R,v,rcsk,sigma0,Te_JD,tau_JD,tr_JD,f0_JD,f1_JD,fs_JD,frame_JD,num_jilei_JD,num_tongdao_JD,num_cankao_JD,num_baohu_JD) ;  

                    case 'AGC????'
                      load data/data_JDParameter;
                 load data/data_AGC;
                 load data/data_target0Parameter;
                          JDmainJ_AGC(fz_JD,B_JD,Gt_JD,Gr_JD,F_JD,Pt_JD,L_JD,Pfa_JD,R,v,rcsk,sigma0,Te_JD,tau_JD,tr_JD,f0_JD,f1_JD,fs_JD,frame_JD,num_jilei_JD,num_tongdao_JD,num_cankao_JD,num_baohu_JD) ;    

                end
          end
end


% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(~, ~, ~)%????????????????
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
antannaParameterSetting;


% --- Executes during object creation, after setting all properties.
function radar_CreateFcn(~, ~, ~)
% hObject    handle to radar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton27.
function pushbutton27_Callback(~, ~, ~)%????????
% hObject    handle to pushbutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%??????????????????


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global RadarS;
global c;
c=3e8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%????????
global sigma;
global rcsk;
global sigma0;
sigma = rcs(rcsk,sigma0);

load data/data_target0Parameter

load data/antannaParameter;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if RadarS==1
        %%%%%%%%%%%%%%????????????
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LFM????????
 load data/data_LFMParameter
        fr=1/tr_LFM;%????????????
      
        lamta=c/fz_LFM;%????
        tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%????????????????????????
        N=length(tm);%????????????????????????????
        An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%????????
        [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
        ts=1/fs_LFM;
        k=B1_LFM/tau_LFM;   
        Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %????????????????
        A=sqrt(Prs);%????????????
        f_doppler=2*v/lamta;%????????????????%???????????????????? 
        [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
        [M,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y/sqrt(Pt_LFM),tr_LFM,ts,N);
        [s_echo_2,~]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
        
        
         [~,RTAz,RTEl] = xyz2radar(radar.x,radar.y,radar.z,...
                    target.x,target.y,target.z);                        %????????????
%         [FTR,FTAz,FTEl] = xyz2radar(radar.x1,radar.y1,radar.z1,...
%                     jammer.x3,jammer.y3,jammer.z3);                        %????????????
        Gt = 10^(radar.Gt/10);                                                           %????????????????
%         Gj = jammer.Gain;                                                          %??????????????????
        Gr = gain(radar,RTAz,RTEl);                                                %??????????????????
%         Grj = gain(radar,FTAz,FTEl)/Gj;/Gt                                               %????????????????????
        echo1.sum=s_echo_2*Gr.Gainsum;
        echo1.Az=s_echo_2*Gr.GainAz;
        echo1.Ez=s_echo_2*Gr.GainEl;
        
        
        
        s_echo_2=echo1.sum;
        [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
        s_echo_1=s_echo_2+s_noise;
        [s_echo_1]=gaofang(f0_LFM,B_LFM,fs_LFM,s_echo_1);    
        [s_echo_1,f0]=hunpin(s_echo_1,N,frame_LFM,f1_LFM,fs_LFM,f0_LFM);
        [s_echo_mf]=jianbo(s_echo_1,N,frame_LFM,f0,fs_LFM);
        [~,~,~]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts);                  
        [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_LFM,match_filter_fft,tau_LFM,D,ts); 
        [s_mtd]=mtd(pc_result1.',M1,num_jilei_LFM,num_tongdao_LFM);
%        CFARsingle=hengxujing(M1,Pfa_LFM,s_mtd,num_cankao_LFM,num_tongdao_LFM,num_baohu_LFM,ts,c,tau_LFM,D);
%        
%        
%        
% [indexV,indexR] = find(abs(CFARsingle)>0);

%% ????????????
AzAngle = -radar.antennabeamwidth/2:0.01:radar.antennabeamwidth/2;
ElAngle = -radar.antennabeamwidth/2:0.01:radar.antennabeamwidth/2;
G = gain(radar,AzAngle,ElAngle);
echo1.El=echo1.Ez;
receivesingle=echo1;
radar.Fs=fs_LFM;
radar.TimeWidth=tau_LFM;
TempNum=fix(radar.Fs*radar.TimeWidth);
[~,NYY0]=find(receivesingle.sum);
NYY=NYY0(1);
su=abs(receivesingle.sum(1,NYY));
az=abs(receivesingle.Az(1,NYY));
ratio=az./su;
AzRatio = min(ratio);

ElRatio = min(abs(receivesingle.El(1,NYY))./abs(receivesingle.sum(1,NYY)));
[~, Azloc] = min((AzRatio - G.NDSAz(:,round(length(AzAngle)/2))).^2);
[~, Elloc] = min((ElRatio -	G.NDSEl(round(length(ElAngle)/2),:)).^2);
AOA.Az = AzAngle(Azloc(:));
AOA.El = AzAngle(Elloc(:));
figure,plot(AOA.Az,AOA.El,'*');grid on;hold on;
x=[0 AOA.Az];y=ones(1,2)*AOA.El;
plot(x,y,'r-.');hold on;
x=[AOA.Az AOA.Az];y=[0 AOA.El];
plot(x,y,'r-.');
nu=[AOA.Az AOA.El];
str2=num2str(nu);
str1='(';str3=')';
str6=cat(2,str1,str2);
str=cat(2,str6,str3);
text(AOA.Az,AOA.El,str);
title('????????????');
xlabel('??????/??');
ylabel('??????/??');
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      title('????????');xlabel('??????/??????');ylabel('??????/??????');





        %%%%%%%%%%%%%%%%%%%%%%%%%
elseif RadarS==2  %%%%%%%%%%%%%????????????
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BF????????
         global code 
 load data/data_BFParameter
         tm_B=0:1/fs_B:tr_B-1/fs_B;%????????????????????????
        N=length(tm_B);%????????????????????????????
        number1=length(code);
        ts=1/fs_B;
        lamta=c/fz_B;
        fr=1/tr_B;
        An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%????????
        Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %????????????????
        A=sqrt(Prs);%????????????
        f_doppler=2*v/lamta;%????????????????
        [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
        [s_echo_2,~]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
        
        
         [~,RTAz,RTEl] = xyz2radar(radar.x,radar.y,radar.z,...
                    target.x,target.y,target.z);                        %????????????
%         [FTR,FTAz,FTEl] = xyz2radar(radar.x1,radar.y1,radar.z1,...
%                     jammer.x3,jammer.y3,jammer.z3);                        %????????????
        Gt = 10^(radar.Gt/10);                                                           %????????????????
%         Gj = jammer.Gain;                                                          %??????????????????
        Gr = gain(radar,RTAz,RTEl);                                                %??????????????????
%         Grj = gain(radar,FTAz,FTEl)/Gj;/Gt                                               %????????????????????
        echo1.sum=s_echo_2*Gr.Gainsum;
        echo1.Az=s_echo_2*Gr.GainAz;
        echo1.Ez=s_echo_2*Gr.GainEl;
        
        
        
        s_echo_2=echo1.sum;
        [M,match_filter_fft]=maiyaxishu(f0_B,fs_B,y/sqrt(Pt_B),tr_B,ts,N);
        [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
        s_echo_1=s_echo_2+s_noise;
        [s_echo_1]=gaofang(f0_B,B_B,fs_B,s_echo_1);
        [s_echo_1,f0]=hunpin(s_echo_1,N,frame_B,f1_B,fs_B,f0_B);
        [s_echo_mf]=jianbo(s_echo_1,N,frame_B,f0,fs_B);
        [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_B,match_filter_fft,tau_B,D,ts);
        [s_mtd]=mtd(pc_result1.',M1,num_jilei_B,num_tongdao_B);
%         hengxujing(M1,Pfa_B,s_mtd,num_cankao_B,num_tongdao_B,num_baohu_B,ts,c,tau_B,D);
        
        
        AzAngle = -radar.antennabeamwidth/2:0.01:radar.antennabeamwidth/2;
        ElAngle = -radar.antennabeamwidth/2:0.01:radar.antennabeamwidth/2;
        G = gain(radar,AzAngle,ElAngle);
        echo1.El=echo1.Ez;
        receivesingle=echo1;
        radar.Fs=fs_B;
        radar.TimeWidth=tau_B;
        TempNum=fix(radar.Fs*radar.TimeWidth);
        [~,NYY0]=find(receivesingle.sum);
        NYY=NYY0(1);
        su=abs(receivesingle.sum(1,NYY));
        az=abs(receivesingle.Az(1,NYY));
        ratio=az./su;
        AzRatio = min(ratio);

        ElRatio = min(abs(receivesingle.El(1,NYY))./abs(receivesingle.sum(1,NYY)));
        [~, Azloc] = min((AzRatio - G.NDSAz(:,round(length(AzAngle)/2))).^2);
        [~, Elloc] = min((ElRatio -	G.NDSEl(round(length(ElAngle)/2),:)).^2);
        AOA.Az = AzAngle(Azloc(:));
        AOA.El = AzAngle(Elloc(:));
        figure,plot(AOA.Az,AOA.El,'*');grid on;hold on;
        x=[0 AOA.Az];
        y=ones(1,2)*AOA.El;
        plot(x,y,'r-.');hold on;
        x=[AOA.Az AOA.Az];
        y=[0 AOA.El];
        plot(x,y,'r-.');
        nu=[AOA.Az AOA.El];
        str2=num2str(nu);
        str1='(';str3=')';
        str6=cat(2,str1,str2);
        str=cat(2,str6,str3);
        text(AOA.Az,AOA.El,str);
        title('????????????');
        xlabel('??????/??');
        ylabel('??????/??');
    else
        %%%%%%%%%%%%%%%%%%%%%%%????????
         load data/data_JDParameter

        fr=1/tr_JD;
        lamta=c/fz_JD;%????
        tm=0:1/fs_JD:tr_JD-1/fs_JD;%????????????????????????
        N=length(tm);%????????????????????????????
        An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%????????
        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        ts=1/fs_JD; 
        Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %????????????????
        A=sqrt(Prs);%????????????
        f_doppler=2*v/lamta;%????????????????%???????????????????? 
        [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
        [M,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y/sqrt(Pt_JD),tr_JD,ts,N);
        [s_echo_2,~]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
        
        
         [~,RTAz,RTEl] = xyz2radar(radar.x,radar.y,radar.z,...
                    target.x,target.y,target.z);                        %????????????
%         [FTR,FTAz,FTEl] = xyz2radar(radar.x1,radar.y1,radar.z1,...
%                     jammer.x3,jammer.y3,jammer.z3);                        %????????????
        Gt = 10^(radar.Gt/10);                                                           %????????????????
%         Gj = jammer.Gain;                                                          %??????????????????
        Gr = gain(radar,RTAz,RTEl);                                                %??????????????????
%         Grj = gain(radar,FTAz,FTEl)/Gj;/Gt                                               %????????????????????
        echo1.sum=s_echo_2*Gr.Gainsum;
        echo1.Az=s_echo_2*Gr.GainAz;
        echo1.Ez=s_echo_2*Gr.GainEl;
        
        
        
        s_echo_2=echo1.sum;
        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        s_echo_1=s_echo_2+s_noise;
        [s_echo_1]=gaofang(f0_JD,B_JD,fs_JD,s_echo_1);
        [s_echo_1,f0]=hunpin(s_echo_1,N,frame_JD,f1_JD,fs_JD,f0_JD);
        [s_echo_mf]=jianbo(s_echo_1,N,frame_JD,f0,fs_JD);
        [pc_result,pc_result1,M1]=maichongyasuo(s_echo_mf,M,frame_JD,match_filter_fft,tau_JD,D,ts);
        [s_mtd]=mtd(pc_result1.',M1,num_jilei_JD,num_tongdao_JD);
%         hengxujing(M1,Pfa_JD,s_mtd,num_cankao_JD,num_tongdao_JD,num_baohu_JD,ts,c,tau_JD,D);
            AzAngle = -radar.antennabeamwidth/2:0.01:radar.antennabeamwidth/2;
            ElAngle = -radar.antennabeamwidth/2:0.01:radar.antennabeamwidth/2;
            G = gain(radar,AzAngle,ElAngle);
            echo1.El=echo1.Ez;
            receivesingle=echo1;
            radar.Fs=fs_JD;
            radar.TimeWidth=tau_JD;
            TempNum=fix(radar.Fs*radar.TimeWidth);
            [~,NYY0]=find(receivesingle.sum);
            NYY=NYY0(1);
            su=abs(receivesingle.sum(1,NYY));
            az=abs(receivesingle.Az(1,NYY));
            ratio=az./su;
            AzRatio = min(ratio);

            ElRatio = min(abs(receivesingle.El(1,NYY))./abs(receivesingle.sum(1,NYY)));
            [~, Azloc] = min((AzRatio - G.NDSAz(:,round(length(AzAngle)/2))).^2);
            [~, Elloc] = min((ElRatio -	G.NDSEl(round(length(ElAngle)/2),:)).^2);
            AOA.Az = AzAngle(Azloc(:));
            AOA.El = AzAngle(Elloc(:));
            figure,plot(AOA.Az,AOA.El,'*');grid on;hold on;
            x=[0 AOA.Az];
            y=ones(1,2)*AOA.El;
            plot(x,y,'r-.');hold on;
            x=[AOA.Az AOA.Az];
            y=[0 AOA.El];
            plot(x,y,'r-.');
            nu=[AOA.Az AOA.El];
            str2=num2str(nu);
            str1='(';str3=')';
            str6=cat(2,str1,str2);
            str=cat(2,str6,str3);
            text(AOA.Az,AOA.El,str);
            title('????????????');
            xlabel('??????/??');
            ylabel('??????/??');

end
