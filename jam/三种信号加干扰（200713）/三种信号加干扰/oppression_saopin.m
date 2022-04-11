function varargout = oppression_saopin(varargin)
% OPPRESSION_SAOPIN MATLAB code for oppression_saopin.fig
%      OPPRESSION_SAOPIN, by itself, creates a new OPPRESSION_SAOPIN or raises the existing
%      singleton*.
%
%      H = OPPRESSION_SAOPIN returns the handle to a new OPPRESSION_SAOPIN or the handle to
%      the existing singleton*.
%
%      OPPRESSION_SAOPIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OPPRESSION_SAOPIN.M with the given input arguments.
%
%      OPPRESSION_SAOPIN('Property','Value',...) creates a new OPPRESSION_SAOPIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before oppression_saopin_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to oppression_saopin_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help oppression_saopin

% Last Modified by GUIDE v2.5 11-Aug-2017 17:00:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @oppression_saopin_OpeningFcn, ...
                   'gui_OutputFcn',  @oppression_saopin_OutputFcn, ...
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


% --- Executes just before oppression_saopin is made visible.
function oppression_saopin_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to oppression_saopin (see VARARGIN)

% Choose default command line output for oppression_saopin
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes oppression_saopin wait for user response (see UIRESUME)
% uiwait(handles.figure1);
set(handles.fs_saopin,'string','160e6');
set(handles.Bj_saopin,'string','20e6');
set(handles.fj_saopin,'string','10e6');
set(handles.frame_saopin,'string','64');
set(handles.Prj_saopin,'string','10');
set(handles.Tr_saopin,'string','40e-6');
set(handles.K_sweep_saopin,'string','1e9');
set(handles.Time_begin_saopin,'string','0');
set(handles.T_fr_saopin,'string','1200e-6');


% --- Outputs from this function are returned to the command line.
function varargout = oppression_saopin_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fs_saopin=str2double(get(handles.fs_saopin,'string'));
Bj_saopin=str2double(get(handles.Bj_saopin,'string'));
fj_saopin=str2double(get(handles.fj_saopin,'string'));
frame_saopin=str2double(get(handles.frame_saopin,'string'));
Prj_saopin=str2double(get(handles.Prj_saopin,'string'));
Tr_saopin=str2double(get(handles.Tr_saopin,'string'));
K_sweep_saopin=str2double(get(handles.K_sweep_saopin,'string'));
Time_begin_saopin=str2double(get(handles.Time_begin_saopin,'string'));
T_fr_saopin=str2double(get(handles.T_fr_saopin,'string'));
save data/data_saopin fs_saopin Bj_saopin fj_saopin frame_saopin Prj_saopin Tr_saopin K_sweep_saopin Time_begin_saopin T_fr_saopin
h1=msgbox('载入目标参数成功');
pause(1);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.fs_saopin,'string',' ');
set(handles.Bj_saopin,'string',' ');
set(handles.fj_saopin,'string',' ');
set(handles.frame_saopin,'string',' ');
set(handles.Prj_saopin,'string',' ');
set(handles.Tr_saopin,'string',' ');
set(handles.K_sweep_saopin,'string',' ');
set(handles.Time_begin_saopin,'string',' ');
set(handles.T_fr_saopin,'string',' ');
h2=msgbox('重置目标参数成功');
pause(1);
close (h2);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 load data/data_saopin
 %T_fr_saopin=2*Tr_saopin
 [ sig_noise,t_noise ] = jam_sweepfrequency( fs_saopin,Bj_saopin,fj_saopin,frame_saopin,Prj_saopin,Tr_saopin,T_fr_saopin,Time_begin_saopin,K_sweep_saopin );
view_jam_sweepfrequency( sig_noise,t_noise,fs_saopin );

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles) %#ok<*INUSL>
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.fs_saopin,'string','160e6');
set(handles.Bj_saopin,'string','20e6');
set(handles.fj_saopin,'string','10e6');
set(handles.frame_saopin,'string','64');
set(handles.Prj_saopin,'string','10');
set(handles.Tr_saopin,'string','40e-6');
set(handles.K_sweep_saopin,'string','1e9');
set(handles.Time_begin_saopin,'string','0');
set(handles.T_fr_saopin,'string','1200e-6');
h3=msgbox('初始化目标参数成功');
pause(1);
close (h3);


function fs_saopin_Callback(hObject, eventdata, handles)
% hObject    handle to fs_saopin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fs_saopin as text
%        str2double(get(hObject,'String')) returns contents of fs_saopin as a double


% --- Executes during object creation, after setting all properties.
function fs_saopin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs_saopin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Bj_saopin_Callback(hObject, eventdata, handles)
% hObject    handle to Bj_saopin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Bj_saopin as text
%        str2double(get(hObject,'String')) returns contents of Bj_saopin as a double


% --- Executes during object creation, after setting all properties.
function Bj_saopin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bj_saopin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fj_saopin_Callback(hObject, eventdata, handles)
% hObject    handle to fj_saopin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fj_saopin as text
%        str2double(get(hObject,'String')) returns contents of fj_saopin as a double


% --- Executes during object creation, after setting all properties.
function fj_saopin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fj_saopin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frame_saopin_Callback(hObject, eventdata, handles)
% hObject    handle to frame_saopin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_saopin as text
%        str2double(get(hObject,'String')) returns contents of frame_saopin as a double


% --- Executes during object creation, after setting all properties.
function frame_saopin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_saopin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Prj_saopin_Callback(hObject, eventdata, handles)
% hObject    handle to Prj_saopin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Prj_saopin as text
%        str2double(get(hObject,'String')) returns contents of Prj_saopin as a double


% --- Executes during object creation, after setting all properties.
function Prj_saopin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Prj_saopin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tr_saopin_Callback(hObject, eventdata, handles)
% hObject    handle to Tr_saopin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tr_saopin as text
%        str2double(get(hObject,'String')) returns contents of Tr_saopin as a double


% --- Executes during object creation, after setting all properties.
function Tr_saopin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tr_saopin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function K_sweep_saopin_Callback(hObject, eventdata, handles)
% hObject    handle to K_sweep_saopin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of K_sweep_saopin as text
%        str2double(get(hObject,'String')) returns contents of K_sweep_saopin as a double


% --- Executes during object creation, after setting all properties.
function K_sweep_saopin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to K_sweep_saopin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Time_begin_saopin_Callback(hObject, eventdata, handles)
% hObject    handle to Time_begin_saopin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Time_begin_saopin as text
%        str2double(get(hObject,'String')) returns contents of Time_begin_saopin as a double


% --- Executes during object creation, after setting all properties.
function Time_begin_saopin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Time_begin_saopin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T_fr_saopin_Callback(hObject, eventdata, handles) %#ok<*DEFNU,*INUSD>
% hObject    handle to T_fr_saopin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_fr_saopin as text
%        str2double(get(hObject,'String')) returns contents of T_fr_saopin as a double


% --- Executes during object creation, after setting all properties.
function T_fr_saopin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_fr_saopin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
