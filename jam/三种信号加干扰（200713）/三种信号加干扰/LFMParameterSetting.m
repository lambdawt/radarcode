function varargout = LFMParameterSetting(varargin)
% LFMPARAMETERSETTING MATLAB code for LFMParameterSetting.fig
%      LFMPARAMETERSETTING, by itself, creates a new LFMPARAMETERSETTING or raises the existing
%      singleton*.
%
%      H = LFMPARAMETERSETTING returns the handle to a new LFMPARAMETERSETTING or the handle to
%      the existing singleton*.
%
%      LFMPARAMETERSETTING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LFMPARAMETERSETTING.M with the given input arguments.
%
%      LFMPARAMETERSETTING('Property','Value',...) creates a new LFMPARAMETERSETTING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LFMParameterSetting_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LFMParameterSetting_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LFMParameterSetting

% Last Modified by GUIDE v2.5 07-Aug-2016 15:12:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LFMParameterSetting_OpeningFcn, ...
                   'gui_OutputFcn',  @LFMParameterSetting_OutputFcn, ...
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


% --- Executes just before LFMParameterSetting is made visible.
function LFMParameterSetting_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LFMParameterSetting (see VARARGIN)

% Choose default command line output for LFMParameterSetting
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes LFMParameterSetting wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LFMParameterSetting_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function B_LFM_Callback(hObject, eventdata, handles) %#ok<*INUSD,*DEFNU>
% hObject    handle to B_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of B_LFM as text
%        str2double(get(hObject,'String')) returns contents of B_LFM as a double


% --- Executes during object creation, after setting all properties.
function B_LFM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Gt_LFM_Callback(hObject, eventdata, handles)
% hObject    handle to Gt_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Gt_LFM as text
%        str2double(get(hObject,'String')) returns contents of Gt_LFM as a double


% --- Executes during object creation, after setting all properties.
function Gt_LFM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gt_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Gr_LFM_Callback(hObject, eventdata, handles)
% hObject    handle to Gr_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Gr_LFM as text
%        str2double(get(hObject,'String')) returns contents of Gr_LFM as a double


% --- Executes during object creation, after setting all properties.
function Gr_LFM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gr_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fz_LFM_Callback(hObject, eventdata, handles)
% hObject    handle to fz_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fz_LFM as text
%        str2double(get(hObject,'String')) returns contents of fz_LFM as a double


% --- Executes during object creation, after setting all properties.
function fz_LFM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fz_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function F_LFM_Callback(hObject, eventdata, handles)
% hObject    handle to F_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of F_LFM as text
%        str2double(get(hObject,'String')) returns contents of F_LFM as a double


% --- Executes during object creation, after setting all properties.
function F_LFM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to F_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pt_LFM_Callback(hObject, eventdata, handles)
% hObject    handle to Pt_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pt_LFM as text
%        str2double(get(hObject,'String')) returns contents of Pt_LFM as a double


% --- Executes during object creation, after setting all properties.
function Pt_LFM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pt_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function L_LFM_Callback(hObject, eventdata, handles)
% hObject    handle to L_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of L_LFM as text
%        str2double(get(hObject,'String')) returns contents of L_LFM as a double


% --- Executes during object creation, after setting all properties.
function L_LFM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to L_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function R_Callback(hObject, eventdata, handles)
% hObject    handle to R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R as text
%        str2double(get(hObject,'String')) returns contents of R as a double


% --- Executes during object creation, after setting all properties.
function R_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function v_Callback(hObject, eventdata, handles)
% hObject    handle to v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of v as text
%        str2double(get(hObject,'String')) returns contents of v as a double


% --- Executes during object creation, after setting all properties.
function v_CreateFcn(hObject, eventdata, handles)
% hObject    handle to v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigma0_Callback(hObject, eventdata, handles)
% hObject    handle to sigma0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigma0 as text
%        str2double(get(hObject,'String')) returns contents of sigma0 as a double


% --- Executes during object creation, after setting all properties.
function sigma0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigma0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Te_LFM_Callback(hObject, eventdata, handles)
% hObject    handle to Te_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Te_LFM as text
%        str2double(get(hObject,'String')) returns contents of Te_LFM as a double


% --- Executes during object creation, after setting all properties.
function Te_LFM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Te_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tau_LFM_Callback(hObject, eventdata, handles)
% hObject    handle to tau_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tau_LFM as text
%        str2double(get(hObject,'String')) returns contents of tau_LFM as a double


% --- Executes during object creation, after setting all properties.
function tau_LFM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tau_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tr_LFM_Callback(hObject, eventdata, handles)
% hObject    handle to tr_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tr_LFM as text
%        str2double(get(hObject,'String')) returns contents of tr_LFM as a double


% --- Executes during object creation, after setting all properties.
function tr_LFM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tr_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f0_LFM_Callback(hObject, eventdata, handles)
% hObject    handle to f0_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f0_LFM as text
%        str2double(get(hObject,'String')) returns contents of f0_LFM as a double


% --- Executes during object creation, after setting all properties.
function f0_LFM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f0_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fs_LFM_Callback(hObject, eventdata, handles)
% hObject    handle to fs_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fs_LFM as text
%        str2double(get(hObject,'String')) returns contents of fs_LFM as a double


% --- Executes during object creation, after setting all properties.
function fs_LFM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs_LFM (see GCBO)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%雷达参数

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%目标参数
% global R;
% global v;
% global sigma0;
% global rcsk;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B_LFM=str2double(get(handles.B_LFM,'string'));  
Gt_LFM=str2double(get(handles.Gt_LFM,'string'));   
Gr_LFM=str2double(get(handles.Gr_LFM,'string'));   
fz_LFM=str2double(get(handles.fz_LFM,'string'));   
F_LFM=str2double(get(handles.F_LFM,'string'));   
Pt_LFM=str2double(get(handles.Pt_LFM,'string')); 
L_LFM=str2double(get(handles.L_LFM,'string')); 
Te_LFM=str2double(get(handles.Te_LFM,'string')); 
tau_LFM=str2double(get(handles.tau_LFM,'string')); 
tr_LFM=str2double(get(handles.tr_LFM,'string')); 
f0_LFM=str2double(get(handles.f0_LFM,'string')); 
fs_LFM=str2double(get(handles.fs_LFM,'string')); 
B1_LFM=str2double(get(handles.B1_LFM,'string'));  
frame_LFM=str2double(get(handles.frame_LFM,'string')); 
f1_LFM=str2double(get(handles.f1_LFM,'string')); 
num_jilei_LFM=str2double(get(handles.num_jilei_LFM,'string'));
num_tongdao_LFM=str2double(get(handles.num_tongdao_LFM,'string'));
num_cankao_LFM=str2double(get(handles.num_cankao_LFM,'string'));
num_baohu_LFM=str2double(get(handles.num_baohu_LFM,'string'));
Pfa_LFM=str2double(get(handles.Pfa_LFM,'string'));
save data/data_LFMParameter B_LFM Gt_LFM Gr_LFM fz_LFM   F_LFM  Pt_LFM L_LFM Te_LFM tau_LFM tr_LFM f0_LFM fs_LFM  B1_LFM frame_LFM f1_LFM num_jilei_LFM num_tongdao_LFM num_cankao_LFM num_baohu_LFM Pfa_LFM
% R=str2double(get(handles.R,'string')); 
% v=str2double(get(handles.v,'string')); 
% sigma0=str2double(get(handles.sigma0,'string')); 
% rcsk=str2double(get(handles.rcsk,'string')); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%生成发射信号
ts=1/fs_LFM;
k=B1_LFM/tau_LFM;                                 %线性调频信号调制系数
tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;  
N=length(tm);
[y]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
figure,plot(0:00001*ts:(N-1)*ts,real(y)),xlabel('t(单位：秒)'), ylabel('y(单位：伏)'),title('LFM发射信号');
figure,plot((0:fs_LFM/length(abs(fft(abs(fftshift(y))))):fs_LFM-fs_LFM/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('LFM发射信号的频谱');




% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function B1_LFM_Callback(hObject, eventdata, handles)
% hObject    handle to B1_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of B1_LFM as text
%        str2double(get(hObject,'String')) returns contents of B1_LFM as a double


% --- Executes during object creation, after setting all properties.
function B1_LFM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B1_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function frame_LFM_Callback(hObject, eventdata, handles)
% hObject    handle to frame_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_LFM as text
%        str2double(get(hObject,'String')) returns contents of frame_LFM as a double


% --- Executes during object creation, after setting all properties.
function frame_LFM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f1_LFM_Callback(hObject, eventdata, handles)
% hObject    handle to f1_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f1_LFM as text
%        str2double(get(hObject,'String')) returns contents of f1_LFM as a double


% --- Executes during object creation, after setting all properties.
function f1_LFM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f1_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_jilei_LFM_Callback(hObject, eventdata, handles)
% hObject    handle to num_jilei_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_jilei_LFM as text
%        str2double(get(hObject,'String')) returns contents of num_jilei_LFM as a double


% --- Executes during object creation, after setting all properties.
function num_jilei_LFM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_jilei_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_tongdao_LFM_Callback(hObject, eventdata, handles)
% hObject    handle to num_tongdao_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_tongdao_LFM as text
%        str2double(get(hObject,'String')) returns contents of num_tongdao_LFM as a double


% --- Executes during object creation, after setting all properties.
function num_tongdao_LFM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_tongdao_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_cankao_LFM_Callback(hObject, eventdata, handles)
% hObject    handle to num_cankao_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_cankao_LFM as text
%        str2double(get(hObject,'String')) returns contents of num_cankao_LFM as a double


% --- Executes during object creation, after setting all properties.
function num_cankao_LFM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_cankao_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_baohu_LFM_Callback(hObject, eventdata, handles)
% hObject    handle to num_baohu_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_baohu_LFM as text
%        str2double(get(hObject,'String')) returns contents of num_baohu_LFM as a double


% --- Executes during object creation, after setting all properties.
function num_baohu_LFM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_baohu_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pfa_LFM_Callback(hObject, eventdata, handles)
% hObject    handle to Pfa_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pfa_LFM as text
%        str2double(get(hObject,'String')) returns contents of Pfa_LFM as a double


% --- Executes during object creation, after setting all properties.
function Pfa_LFM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pfa_LFM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load data/data_LFMParameter
set(handles.B_LFM,'string',B_LFM);
set(handles.Gt_LFM,'string',Gt_LFM);
set(handles.Gr_LFM,'string',Gr_LFM);
set(handles.fz_LFM,'string',fz_LFM);
set(handles.F_LFM,'string',F_LFM);
set(handles.Pt_LFM,'string',Pt_LFM);
set(handles.L_LFM,'string',L_LFM);
set(handles.B1_LFM,'string',B1_LFM);
set(handles.Te_LFM,'string',Te_LFM);
set(handles.tau_LFM,'string',tau_LFM);
set(handles.tr_LFM,'string',tr_LFM);
set(handles.f0_LFM,'string',f0_LFM);
set(handles.fs_LFM,'string',fs_LFM);
set(handles.Pfa_LFM,'string',Pfa_LFM);
set(handles.frame_LFM,'string',frame_LFM);
set(handles.f1_LFM,'string',f1_LFM);
set(handles.num_jilei_LFM,'string',num_jilei_LFM);
set(handles.num_tongdao_LFM,'string',num_tongdao_LFM);
set(handles.num_cankao_LFM,'string',num_cankao_LFM);
set(handles.num_baohu_LFM,'string',num_baohu_LFM);

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.B_LFM,'string','');
set(handles.Gt_LFM,'string','');
set(handles.Gr_LFM,'string','');
set(handles.fz_LFM,'string','');
set(handles.F_LFM,'string','');
set(handles.Pt_LFM,'string','');
set(handles.L_LFM,'string','');
set(handles.B1_LFM,'string','');
set(handles.Te_LFM,'string','');
set(handles.tau_LFM,'string','');
set(handles.tr_LFM,'string','');
set(handles.f0_LFM,'string','');
set(handles.fs_LFM,'string','');
set(handles.Pfa_LFM,'string','');
set(handles.frame_LFM,'string','');
set(handles.f1_LFM,'string','');
set(handles.num_jilei_LFM,'string','');
set(handles.num_tongdao_LFM,'string','');
set(handles.num_cankao_LFM,'string','');
set(handles.num_baohu_LFM,'string','');
