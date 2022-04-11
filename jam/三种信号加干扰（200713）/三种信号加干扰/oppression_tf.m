function varargout = oppression_tf(varargin)
% OPPRESSION_TF MATLAB code for oppression_tf.fig
%      OPPRESSION_TF, by itself, creates a new OPPRESSION_TF or raises the existing
%      singleton*.
%
%      H = OPPRESSION_TF returns the handle to a new OPPRESSION_TF or the handle to
%      the existing singleton*.
%
%      OPPRESSION_TF('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OPPRESSION_TF.M with the given input arguments.
%
%      OPPRESSION_TF('Property','Value',...) creates a new OPPRESSION_TF or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before oppression_tf_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to oppression_tf_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help oppression_tf

% Last Modified by GUIDE v2.5 13-Oct-2016 20:42:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @oppression_tf_OpeningFcn, ...
                   'gui_OutputFcn',  @oppression_tf_OutputFcn, ...
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


% --- Executes just before oppression_tf is made visible.
function oppression_tf_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to oppression_tf (see VARARGIN)

% Choose default command line output for oppression_tf
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes oppression_tf wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = oppression_tf_OutputFcn(hObject, eventdata, handles) 
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

fs_tf=str2double(get(handles.fs_tf,'string'));
Bj_tf=str2double(get(handles.Bj_tf,'string'));
fj_tf=str2double(get(handles.fj_tf,'string'));
frame_tf=str2double(get(handles.frame_tf,'string'));
Prj_tf=str2double(get(handles.Prj_tf,'string'));
Tr_tf=str2double(get(handles.Tr_tf,'string'));

save data/data_tf fs_tf Bj_tf fj_tf frame_tf Prj_tf Tr_tf
h1=msgbox('载入目标参数成功');
pause(1);
close(h1);
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.fs_tf,'string','');
set(handles.Bj_tf,'string','');
set(handles.fj_tf,'string','');
set(handles.frame_tf,'string','');
set(handles.Prj_tf,'string','');
set(handles.Tr_tf,'string','');
h2=msgbox('重置目标参数成功');
pause(1);
close(h2);
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.fs_tf,'string','160e6');
set(handles.Bj_tf,'string','20e6');
set(handles.fj_tf,'string','10e6');
set(handles.frame_tf,'string','64');
set(handles.Prj_tf,'string','40');
set(handles.Tr_tf,'string','40e-6');
h3=msgbox('初始化目标参数成功');
pause(1);
close (h3);
% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

  load data/data_tf
 view_zaoshengtiaofu(fs_tf,Bj_tf,fj_tf,Prj_tf,Tr_tf,frame_tf);

function fs_tf_Callback(hObject, eventdata, handles)
% hObject    handle to fs_tf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fs_tf as text
%        str2double(get(hObject,'String')) returns contents of fs_tf as a double


% --- Executes during object creation, after setting all properties.
function fs_tf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs_tf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Bj_tf_Callback(hObject, eventdata, handles)
% hObject    handle to Bj_tf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Bj_tf as text
%        str2double(get(hObject,'String')) returns contents of Bj_tf as a double


% --- Executes during object creation, after setting all properties.
function Bj_tf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bj_tf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fj_tf_Callback(hObject, eventdata, handles)
% hObject    handle to fj_tf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fj_tf as text
%        str2double(get(hObject,'String')) returns contents of fj_tf as a double


% --- Executes during object creation, after setting all properties.
function fj_tf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fj_tf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frame_tf_Callback(hObject, eventdata, handles)
% hObject    handle to frame_tf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_tf as text
%        str2double(get(hObject,'String')) returns contents of frame_tf as a double


% --- Executes during object creation, after setting all properties.
function frame_tf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_tf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Prj_tf_Callback(hObject, eventdata, handles)
% hObject    handle to Prj_tf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Prj_tf as text
%        str2double(get(hObject,'String')) returns contents of Prj_tf as a double


% --- Executes during object creation, after setting all properties.
function Prj_tf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Prj_tf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tr_tf_Callback(hObject, eventdata, handles)
% hObject    handle to Tr_tf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tr_tf as text
%        str2double(get(hObject,'String')) returns contents of Tr_tf as a double


% --- Executes during object creation, after setting all properties.
function Tr_tf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tr_tf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
