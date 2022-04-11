function varargout = oppression_smart(varargin)
% OPPRESSION_SMART MATLAB code for oppression_smart.fig
%      OPPRESSION_SMART, by itself, creates a new OPPRESSION_SMART or raises the existing
%      singleton*.
%
%      H = OPPRESSION_SMART returns the handle to a new OPPRESSION_SMART or the handle to
%      the existing singleton*.
%
%      OPPRESSION_SMART('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OPPRESSION_SMART.M with the given input arguments.
%
%      OPPRESSION_SMART('Property','Value',...) creates a new OPPRESSION_SMART or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before oppression_smart_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to oppression_smart_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help oppression_smart

% Last Modified by GUIDE v2.5 14-Oct-2016 20:31:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @oppression_smart_OpeningFcn, ...
                   'gui_OutputFcn',  @oppression_smart_OutputFcn, ...
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


% --- Executes just before oppression_smart is made visible.
function oppression_smart_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to oppression_smart (see VARARGIN)

% Choose default command line output for oppression_smart
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global sigma0;
sigma0 = 2;

% UIWAIT makes oppression_smart wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = oppression_smart_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function fs_smart_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to fs_smart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fs_smart as text
%        str2double(get(hObject,'String')) returns contents of fs_smart as a double


% --- Executes during object creation, after setting all properties.
function fs_smart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs_smart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Bj_smart_Callback(hObject, eventdata, handles)
% hObject    handle to Bj_smart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Bj_smart as text
%        str2double(get(hObject,'String')) returns contents of Bj_smart as a double


% --- Executes during object creation, after setting all properties.
function Bj_smart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bj_smart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fj_smart_Callback(hObject, eventdata, handles)
% hObject    handle to fj_smart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fj_smart as text
%        str2double(get(hObject,'String')) returns contents of fj_smart as a double


% --- Executes during object creation, after setting all properties.
function fj_smart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fj_smart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frame_smart_Callback(hObject, eventdata, handles)
% hObject    handle to frame_smart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_smart as text
%        str2double(get(hObject,'String')) returns contents of frame_smart as a double


% --- Executes during object creation, after setting all properties.
function frame_smart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_smart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Prj_smart_Callback(hObject, eventdata, handles)
% hObject    handle to Prj_smart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Prj_smart as text
%        str2double(get(hObject,'String')) returns contents of Prj_smart as a double


% --- Executes during object creation, after setting all properties.
function Prj_smart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Prj_smart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tr_smart_Callback(hObject, eventdata, handles)
% hObject    handle to Tr_smart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tr_smart as text
%        str2double(get(hObject,'String')) returns contents of Tr_smart as a double


% --- Executes during object creation, after setting all properties.
function Tr_smart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tr_smart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fs_smart=str2double(get(handles.fs_smart,'string'));
Bj_smart=str2double(get(handles.Bj_smart,'string'));
fj_smart=str2double(get(handles.fj_smart,'string'));
frame_smart=str2double(get(handles.frame_smart,'string'));
Prj_smart=str2double(get(handles.Prj_smart,'string'));
Tr_smart=str2double(get(handles.Tr_smart,'string'));
save data/data_smart fs_smart Bj_smart fj_smart frame_smart Prj_smart Tr_smart
h1=msgbox('载入目标参数成功');
pause(1);
close(h1);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.fs_smart,'string','');
set(handles.Bj_smart,'string','');
set(handles.fj_smart,'string','');
set(handles.frame_smart,'string','');
set(handles.Prj_smart,'string','');
set(handles.Tr_smart,'string','');
h2=msgbox('重置目标参数成功');
pause(1);
close(h2);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.fs_smart,'string','160e6');
set(handles.Bj_smart,'string','20e6');
set(handles.fj_smart,'string','10e6');
set(handles.frame_smart,'string','64');
set(handles.Prj_smart,'string','40');
set(handles.Tr_smart,'string','40e-6');
h3=msgbox('初始化目标参数成功');
pause(1);
close (h3);

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load data/data_smart %#ok<*LOAD>
Kfm=4e6;tau=1e-6; c=3e8;
global sigma;
global sigma0;
global rcsk;
sigma = rcs(rcsk,sigma0);
global str;
global RadarS;
global code;
switch str{RadarS}
    case '线性调频信号'
     
        
% Tr=40e-6;tau=1e-6;frame=64;fs=160e6;f0=10e6;B1=4e6;A=5;R=2e3;tm=0:1/fs:Tr-1/fs;k=B1/tau;N=length(tm);Pn=(B1/(2.5*Kfm))^2;Bn=B1/2;
        load data/data_LFMParameter;
        load data/data_target0Parameter
        lamta=c/fz_LFM;%波长 
%         ts=1/fs_LFM;
        k=B1_LFM/tau_LFM;                                 %线性调频信号调制系数
        tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;  
        N=length(tm);
        f_doppler=2*v/lamta;%真目标多普勒频率
        Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %目标回波信号功率
        A=sqrt(Prs);%回波信号幅度
        [vRadarSig]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
        Pn=(B1_LFM/(2.5*Kfm))^2;Bn=B1_LFM/2;
% [vRadarSig]=LFMhuiboxinhao(R,c,A,N,frame,fs,0,tm,f0,B1,tau,k);
        [vSmartNoiseSig]=jam_smartnoise( vRadarSig,Pn,Prj_smart,Bn,Kfm,fs_smart );
        view_SmartNoise( vSmartNoiseSig,fs_smart,tau,R); 
   case '相位编码信号'
       load data/data_BFParameter;
       load data/data_target0Parameter 
       tm_B=0:1/fs_B:tr_B-1/fs_B;%一个脉冲重复周期采样序列
        number1=length(code);
        N=length(tm_B);%一个脉冲重复周期采样点数长度
        ts=1/fs_B;
        lamta=c/fz_B;
        Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %目标回波信号功率
        A=sqrt(Prs);%回波信号幅度
        f_doppler=2*v/lamta;%真目标多普勒频率
         Pn=(B_B/(2.5*Kfm))^2;Bn=B_B/2;
         [~,y1]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
        [vRadarSig]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B);
        [vSmartNoiseSig]=jam_smartnoise( vRadarSig,Pn,Prj_smart,Bn,Kfm,fs_smart );
        view_SmartNoise( vSmartNoiseSig,fs_smart,tau,R); 
   case '简单脉冲信号'
       load data/data_JDParameter;
       load data/data_target0Parameter
%         ts=1/fs_JD;
        lamta=c/fz_JD;
        f_doppler=2*v/lamta;
        Prs=(Pt_JD*(10^((Gt_JD*0.1)))*(10^((Gr_JD*0.1)))*(lamta^2)*sigma)/((4*pi)^3*R^4*10^(L_JD*0.1)); %目标回波信号功率
        A=sqrt(Prs);%回波信号幅度
        tm=0:1/fs_JD:tr_JD-1/fs_JD;%一个脉冲重复周期采样序列
        N=length(tm);%一个脉冲重复周期采样点数长度
%       [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
        Pn=(B_JD/(2.5*Kfm))^2;Bn=B_JD/2;
        [vRadarSig]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD);
        [vSmartNoiseSig]=jam_smartnoise( vRadarSig,Pn,Prj_smart,Bn,Kfm,fs_smart );
        view_SmartNoise( vSmartNoiseSig,fs_smart,tau,R); 
end
