function varargout = oppression_sp(varargin)
% OPPRESSION_SP MATLAB code for oppression_sp.fig
%      OPPRESSION_SP, by itself, creates a new OPPRESSION_SP or raises the existing
%      singleton*.
%
%      H = OPPRESSION_SP returns the handle to a new OPPRESSION_SP or the handle to
%      the existing singleton*.
%
%      OPPRESSION_SP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OPPRESSION_SP.M with the given input arguments.
%
%      OPPRESSION_SP('Property','Value',...) creates a new OPPRESSION_SP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before oppression_sp_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to oppression_sp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help oppression_sp

% Last Modified by GUIDE v2.5 13-Oct-2016 20:10:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @oppression_sp_OpeningFcn, ...
                   'gui_OutputFcn',  @oppression_sp_OutputFcn, ...
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


% --- Executes just before oppression_sp is made visible.
function oppression_sp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to oppression_sp (see VARARGIN)

% Choose default command line output for oppression_sp
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes oppression_sp wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = oppression_sp_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function fs_sp_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to fs_sp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fs_sp as text
%        str2double(get(hObject,'String')) returns contents of fs_sp as a double


% --- Executes during object creation, after setting all properties.
function fs_sp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs_sp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Bj_sp_Callback(hObject, eventdata, handles)
% hObject    handle to Bj_sp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Bj_sp as text
%        str2double(get(hObject,'String')) returns contents of Bj_sp as a double


% --- Executes during object creation, after setting all properties.
function Bj_sp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bj_sp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fj_sp_Callback(hObject, eventdata, handles)
% hObject    handle to fj_sp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fj_sp as text
%        str2double(get(hObject,'String')) returns contents of fj_sp as a double


% --- Executes during object creation, after setting all properties.
function fj_sp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fj_sp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frame_sp_Callback(hObject, eventdata, handles)
% hObject    handle to frame_sp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_sp as text
%        str2double(get(hObject,'String')) returns contents of frame_sp as a double


% --- Executes during object creation, after setting all properties.
function frame_sp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_sp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Prj_sp_Callback(hObject, eventdata, handles)
% hObject    handle to Prj_sp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Prj_sp as text
%        str2double(get(hObject,'String')) returns contents of Prj_sp as a double


% --- Executes during object creation, after setting all properties.
function Prj_sp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Prj_sp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tr_sp_Callback(hObject, eventdata, handles)
% hObject    handle to Tr_sp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tr_sp as text
%        str2double(get(hObject,'String')) returns contents of Tr_sp as a double


% --- Executes during object creation, after setting all properties.
function Tr_sp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tr_sp (see GCBO)
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

fs_sp=str2double(get(handles.fs_sp,'string'));
Bj_sp=str2double(get(handles.Bj_sp,'string'));
fj_sp=str2double(get(handles.fj_sp,'string'));
frame_sp=str2double(get(handles.frame_sp,'string'));
Prj_sp=str2double(get(handles.Prj_sp,'string'));
Tr_sp=str2double(get(handles.Tr_sp,'string'));
save data/data_sp fs_sp Bj_sp fj_sp frame_sp Prj_sp Tr_sp

h1=msgbox('载入目标参数成功');
pause(1);
close(h1);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.fs_sp,'string','');
set(handles.Bj_sp,'string','');
set(handles.fj_sp,'string','');
set(handles.frame_sp,'string','');
set(handles.Prj_sp,'string','');
set(handles.Tr_sp,'string','');
h2=msgbox('重置目标参数成功');
pause(1);
close(h2);
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.fs_sp,'string','160e6');
set(handles.Bj_sp,'string','20e6');
set(handles.fj_sp,'string','10e6');
set(handles.frame_sp,'string','64');
set(handles.Prj_sp,'string','40');
set(handles.Tr_sp,'string','40e-6');
h3=msgbox('初始化目标参数成功');
pause(1);
close (h3);

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 load data/data_sp
view_shepinzaosheng(fs_sp,Bj_sp,fj_sp,frame_sp,Prj_sp,Tr_sp);
