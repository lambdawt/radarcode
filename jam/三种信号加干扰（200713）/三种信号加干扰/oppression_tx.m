function varargout = oppression_tx(varargin)
% OPPRESSION_TX MATLAB code for oppression_tx.fig
%      OPPRESSION_TX, by itself, creates a new OPPRESSION_TX or raises the existing
%      singleton*.
%
%      H = OPPRESSION_TX returns the handle to a new OPPRESSION_TX or the handle to
%      the existing singleton*.
%
%      OPPRESSION_TX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OPPRESSION_TX.M with the given input arguments.
%
%      OPPRESSION_TX('Property','Value',...) creates a new OPPRESSION_TX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before oppression_tx_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to oppression_tx_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help oppression_tx

% Last Modified by GUIDE v2.5 13-Oct-2016 16:59:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @oppression_tx_OpeningFcn, ...
                   'gui_OutputFcn',  @oppression_tx_OutputFcn, ...
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


% --- Executes just before oppression_tx is made visible.
function oppression_tx_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to oppression_tx (see VARARGIN)

% Choose default command line output for oppression_tx
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes oppression_tx wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = oppression_tx_OutputFcn(hObject, eventdata, handles) 
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

fs_tx=str2double(get(handles.fs_tx,'string'));
Bj_tx=str2double(get(handles.Bj_tx,'string'));
fj_tx=str2double(get(handles.fj_tx,'string'));
frame_tx=str2double(get(handles.frame_tx,'string'));
Prj_tx=str2double(get(handles.Prj_tx,'string'));
Tr_tx=str2double(get(handles.Tr_tx,'string'));
save data/data_tx fs_tx Bj_tx fj_tx frame_tx Prj_tx Tr_tx
h1=msgbox('载入目标参数成功');
pause(1);
close(h1);

function fs_tx_Callback(hObject, eventdata, handles)
% hObject    handle to fs_tx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fs_tx as text
%        str2double(get(hObject,'String')) returns contents of fs_tx as a double


% --- Executes during object creation, after setting all properties.
function fs_tx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs_tx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Bj_tx_Callback(hObject, eventdata, handles)
% hObject    handle to Bj_tx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Bj_tx as text
%        str2double(get(hObject,'String')) returns contents of Bj_tx as a double


% --- Executes during object creation, after setting all properties.
function Bj_tx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bj_tx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fj_tx_Callback(hObject, eventdata, handles)
% hObject    handle to fj_tx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fj_tx as text
%        str2double(get(hObject,'String')) returns contents of fj_tx as a double


% --- Executes during object creation, after setting all properties.
function fj_tx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fj_tx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frame_tx_Callback(hObject, eventdata, handles)
% hObject    handle to frame_tx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_tx as text
%        str2double(get(hObject,'String')) returns contents of frame_tx as a double


% --- Executes during object creation, after setting all properties.
function frame_tx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_tx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Prj_tx_Callback(hObject, eventdata, handles)
% hObject    handle to Prj_tx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Prj_tx as text
%        str2double(get(hObject,'String')) returns contents of Prj_tx as a double


% --- Executes during object creation, after setting all properties.
function Prj_tx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Prj_tx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tr_tx_Callback(hObject, eventdata, handles)
% hObject    handle to Tr_tx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tr_tx as text
%        str2double(get(hObject,'String')) returns contents of Tr_tx as a double


% --- Executes during object creation, after setting all properties.
function Tr_tx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tr_tx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.fs_tx,'string','');
set(handles.Bj_tx,'string','');
set(handles.fj_tx,'string','');
set(handles.frame_tx,'string','');
set(handles.Prj_tx,'string','');
set(handles.Tr_tx,'string','');
h2=msgbox('重置目标参数成功');
pause(1);
close(h2);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.fs_tx,'string','160e6');
set(handles.Bj_tx,'string','20e6');
set(handles.fj_tx,'string','10e6');
set(handles.frame_tx,'string','64');
set(handles.Prj_tx,'string','40');
set(handles.Tr_tx,'string','40e-6');
h3=msgbox('初始化目标参数成功');
pause(1);
close (h3);
% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  load data/data_tx
 view_zaoshengtiaoxiang(fs_tx,Bj_tx,fj_tx,Prj_tx,Tr_tx,frame_tx);
