function varargout = oppression_shuzhuangpu(varargin)
% OPPRESSION_SHUZHUANGPU MATLAB code for oppression_shuzhuangpu.fig
%      OPPRESSION_SHUZHUANGPU, by itself, creates a new OPPRESSION_SHUZHUANGPU or raises the existing
%      singleton*.
%
%      H = OPPRESSION_SHUZHUANGPU returns the handle to a new OPPRESSION_SHUZHUANGPU or the handle to
%      the existing singleton*.
%
%      OPPRESSION_SHUZHUANGPU('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OPPRESSION_SHUZHUANGPU.M with the given input arguments.
%
%      OPPRESSION_SHUZHUANGPU('Property','Value',...) creates a new OPPRESSION_SHUZHUANGPU or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before oppression_shuzhuangpu_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to oppression_shuzhuangpu_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help oppression_shuzhuangpu

% Last Modified by GUIDE v2.5 15-Oct-2016 10:06:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @oppression_shuzhuangpu_OpeningFcn, ...
                   'gui_OutputFcn',  @oppression_shuzhuangpu_OutputFcn, ...
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


% --- Executes just before oppression_shuzhuangpu is made visible.
function oppression_shuzhuangpu_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to oppression_shuzhuangpu (see VARARGIN)

% Choose default command line output for oppression_shuzhuangpu
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes oppression_shuzhuangpu wait for user response (see UIRESUME)
% uiwait(handles.figure1);
set(handles.fs_shuzhuangpu,'string','160e6');
set(handles.Bj_shuzhuangpu,'string','20e6');
set(handles.Ns_shuzhuangpu,'string','3');
set(handles.frame_shuzhuangpu,'string','64');
set(handles.Prj_shuzhuangpu,'string','10');
set(handles.Tr_shuzhuangpu,'string','40e-6');

% --- Outputs from this function are returned to the command line.
function varargout = oppression_shuzhuangpu_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fs_shuzhuangpu=str2double(get(handles.fs_shuzhuangpu,'string'));
Bj_shuzhuangpu=str2double(get(handles.Bj_shuzhuangpu,'string'));
Ns_shuzhuangpu=str2double(get(handles.Ns_shuzhuangpu,'string'));
frame_shuzhuangpu=str2double(get(handles.frame_shuzhuangpu,'string'));
Prj_shuzhuangpu=str2double(get(handles.Prj_shuzhuangpu,'string'));
Tr_shuzhuangpu=str2double(get(handles.Tr_shuzhuangpu,'string'));
save data/data_shuzhuangpu fs_shuzhuangpu Bj_shuzhuangpu Ns_shuzhuangpu frame_shuzhuangpu Prj_shuzhuangpu Tr_shuzhuangpu
h3=msgbox('初始化目标参数成功');
pause(1);
close (h3);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.fs_shuzhuangpu,'string','160e6');
set(handles.Bj_shuzhuangpu,'string','20e6');
set(handles.Ns_shuzhuangpu,'string','3');
set(handles.frame_shuzhuangpu,'string','64');
set(handles.Prj_shuzhuangpu,'string','10');
set(handles.Tr_shuzhuangpu,'string','40e-6');
h3=msgbox('初始化目标参数成功');
pause(1);
close (h3);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.fs_shuzhuangpu,'string',' ');
set(handles.Bj_shuzhuangpu,'string',' ');
set(handles.Ns_shuzhuangpu,'string',' ');
set(handles.frame_shuzhuangpu,'string',' ');
set(handles.Prj_shuzhuangpu,'string',' ');
set(handles.Tr_shuzhuangpu,'string',' ');
h2=msgbox('重置目标参数成功');
pause(1);
close(h2);
% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load data/data_shuzhuangpu
fj=[1e7,4e7,7e7];
[sig_noise,t_noise] = jam_combspectrum(fs_shuzhuangpu,Bj_shuzhuangpu,Ns_shuzhuangpu,fj,frame_shuzhuangpu,Prj_shuzhuangpu,Tr_shuzhuangpu);
view_jam_combspectrum( sig_noise,t_noise,fs_shuzhuangpu );

function fs_shuzhuangpu_Callback(hObject, eventdata, handles)
% hObject    handle to fs_shuzhuangpu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fs_shuzhuangpu as text
%        str2double(get(hObject,'String')) returns contents of fs_shuzhuangpu as a double


% --- Executes during object creation, after setting all properties.
function fs_shuzhuangpu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs_shuzhuangpu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Bj_shuzhuangpu_Callback(hObject, eventdata, handles)
% hObject    handle to Bj_shuzhuangpu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Bj_shuzhuangpu as text
%        str2double(get(hObject,'String')) returns contents of Bj_shuzhuangpu as a double


% --- Executes during object creation, after setting all properties.
function Bj_shuzhuangpu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bj_shuzhuangpu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frame_shuzhuangpu_Callback(hObject, eventdata, handles)
% hObject    handle to frame_shuzhuangpu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_shuzhuangpu as text
%        str2double(get(hObject,'String')) returns contents of frame_shuzhuangpu as a double


% --- Executes during object creation, after setting all properties.
function frame_shuzhuangpu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_shuzhuangpu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Prj_shuzhuangpu_Callback(hObject, eventdata, handles)
% hObject    handle to Prj_shuzhuangpu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Prj_shuzhuangpu as text
%        str2double(get(hObject,'String')) returns contents of Prj_shuzhuangpu as a double


% --- Executes during object creation, after setting all properties.
function Prj_shuzhuangpu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Prj_shuzhuangpu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tr_shuzhuangpu_Callback(hObject, eventdata, handles)
% hObject    handle to Tr_shuzhuangpu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tr_shuzhuangpu as text
%        str2double(get(hObject,'String')) returns contents of Tr_shuzhuangpu as a double


% --- Executes during object creation, after setting all properties.
function Tr_shuzhuangpu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tr_shuzhuangpu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ns_shuzhuangpu_Callback(hObject, eventdata, handles)
% hObject    handle to Ns_shuzhuangpu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ns_shuzhuangpu as text
%        str2double(get(hObject,'String')) returns contents of Ns_shuzhuangpu as a double


% --- Executes during object creation, after setting all properties.
function Ns_shuzhuangpu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ns_shuzhuangpu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
