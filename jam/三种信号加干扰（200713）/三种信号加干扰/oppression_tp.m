function varargout = oppression_tp(varargin)
%OPPRESSION_TP M-file for oppression_tp.fig
%      OPPRESSION_TP, by itself, creates a new OPPRESSION_TP or raises the existing
%      singleton*.
%
%      H = OPPRESSION_TP returns the handle to a new OPPRESSION_TP or the handle to
%      the existing singleton*.
%
%      OPPRESSION_TP('Property','Value',...) creates a new OPPRESSION_TP using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to oppression_tp_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      OPPRESSION_TP('CALLBACK') and OPPRESSION_TP('CALLBACK',hObject,...) call the
%      local function named CALLBACK in OPPRESSION_TP.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help oppression_tp

% Last Modified by GUIDE v2.5 28-May-2020 14:50:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @oppression_tp_OpeningFcn, ...
                   'gui_OutputFcn',  @oppression_tp_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before oppression_tp is made visible.
function oppression_tp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for oppression_tp
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes oppression_tp wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = oppression_tp_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fs_tp=str2double(get(handles.fs_tp,'string'));
Bj_tp=str2double(get(handles.Bj_tp,'string'));
fj_tp=str2double(get(handles.fj_tp,'string'));
frame_tp=str2double(get(handles.frame_tp,'string'));
Prj_tp=str2double(get(handles.Prj_tp,'string'));
Tr_tp=str2double(get(handles.Tr_tp,'string'));
Kfm=str2double(get(handles.Kfm,'string'));
Bn=str2double(get(handles.Bn,'string'));
save data/data_tp fs_tp Bj_tp fj_tp frame_tp Prj_tp Tr_tp Kfm Bn
h1=msgbox('载入目标参数成功');
pause(1);
close(h1);
% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.fs_tp,'string','');
set(handles.Bj_tp,'string','');
set(handles.fj_tp,'string','');
set(handles.frame_tp,'string','');
set(handles.Prj_tp,'string','');
set(handles.Tr_tp,'string','');
set(handles.Kfm,'string','');
set(handles.Bn,'string','');
h2=msgbox('重置目标参数成功');
pause(1);
close(h2);
% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.fs_tp,'string','160e6');
set(handles.Bj_tp,'string','20e6');
set(handles.fj_tp,'string','10e6');
set(handles.frame_tp,'string','64');
set(handles.Prj_tp,'string','40');
set(handles.Tr_tp,'string','40e-6');
set(handles.Kfm,'string','4e6');
set(handles.Bn,'string','10e6');
h3=msgbox('初始化目标参数成功');
pause(1);
close (h3);
% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 load data/data_tp
% zaoshengtiaopin(fs_tp,Kfm,Prj_tp,Bj_tp,Bn,fj_tp,frame_tp,Tr_tp)
Pn=(Bj_tp/2/(2.5*Kfm))^2;
view_zaoshengtiaopin(fs_tp,Kfm,Prj_tp,Pn,Bn,fj_tp,frame_tp,Tr_tp);

function fs_tp_Callback(hObject, eventdata, handles)
% hObject    handle to fs_tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fs_tp as text
%        str2double(get(hObject,'String')) returns contents of fs_tp as a double


% --- Executes during object creation, after setting all properties.
function fs_tp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs_tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Bj_tp_Callback(hObject, eventdata, handles)
% hObject    handle to Bj_tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Bj_tp as text
%        str2double(get(hObject,'String')) returns contents of Bj_tp as a double


% --- Executes during object creation, after setting all properties.
function Bj_tp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bj_tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fj_tp_Callback(hObject, eventdata, handles)
% hObject    handle to fj_tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fj_tp as text
%        str2double(get(hObject,'String')) returns contents of fj_tp as a double


% --- Executes during object creation, after setting all properties.
function fj_tp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fj_tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frame_tp_Callback(hObject, eventdata, handles)
% hObject    handle to frame_tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_tp as text
%        str2double(get(hObject,'String')) returns contents of frame_tp as a double


% --- Executes during object creation, after setting all properties.
function frame_tp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Prj_tp_Callback(hObject, eventdata, handles)
% hObject    handle to Prj_tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Prj_tp as text
%        str2double(get(hObject,'String')) returns contents of Prj_tp as a double


% --- Executes during object creation, after setting all properties.
function Prj_tp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Prj_tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tr_tp_Callback(hObject, eventdata, handles)
% hObject    handle to Tr_tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tr_tp as text
%        str2double(get(hObject,'String')) returns contents of Tr_tp as a double


% --- Executes during object creation, after setting all properties.
function Tr_tp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tr_tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Kfm_Callback(hObject, eventdata, handles)
% hObject    handle to Kfm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Kfm as text
%        str2double(get(hObject,'String')) returns contents of Kfm as a double


% --- Executes during object creation, after setting all properties.
function Kfm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Kfm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Bn_Callback(hObject, eventdata, handles)
% hObject    handle to Bn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Bn as text
%        str2double(get(hObject,'String')) returns contents of Bn as a double


% --- Executes during object creation, after setting all properties.
function Bn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on pushbutton8 and none of its controls.
function pushbutton8_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
