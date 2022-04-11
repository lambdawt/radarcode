function varargout = antannaParameterSetting(varargin)
% antannaPARAMETERSETTING MATLAB code for antannaParameterSetting.fig
%      antannaPARAMETERSETTING, by itself, creates a new antannaPARAMETERSETTING or raises the existing
%      singleton*.
%
%      H = antannaPARAMETERSETTING returns the handle to a new antannaPARAMETERSETTING or the handle to
%      the existing singleton*.
%
%      antannaPARAMETERSETTING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in antannaPARAMETERSETTING.M with the given input arguments.
%
%      antannaPARAMETERSETTING('Property','Value',...) creates a new antannaPARAMETERSETTING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before antannaParameterSetting_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to antannaParameterSetting_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help antannaParameterSetting

% Last Modified by GUIDE v2.5 06-Jul-2020 10:03:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @antannaParameterSetting_OpeningFcn, ...
                   'gui_OutputFcn',  @antannaParameterSetting_OutputFcn, ...
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


% --- Executes just before antannaParameterSetting is made visible.
function antannaParameterSetting_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to antannaParameterSetting (see VARARGIN)

% Choose default command line output for antannaParameterSetting
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes antannaParameterSetting wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = antannaParameterSetting_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in saveParameter.
function saveParameter_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to saveParameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
 radar.x=str2num(get(handles.x_antanna,'string'));
 radar.y=str2num(get(handles.y_antanna,'string'));
 radar.z=str2num(get(handles.z_antanna,'string'));
 radar.Gt=str2num(get(handles.Gt_antanna,'string'));
 radar.Gr=str2num(get(handles.Gr_antanna,'string'));
 radar.antennaAzAngle=str2num(get(handles.AzAngle_antanna,'string'));
 radar.antennaELAngle=str2num(get(handles.ELAngle_antanna,'string'));
 radar.antennabeamwidth=str2num(get(handles.beamwidth_antanna,'string'));
 radar.AdjacentAngle=str2num(get(handles.AdjacentAngle_antanna,'string'));
 radar.strfunc=get(handles.func_antenna,'string');
 radar.nfunc=get(handles.func_antenna,'value');
 %radar.nanglemeasureMode=get(handles.anglemeasureMode_antanna,'value');
 radar.strfunc{radar.nfunc}; %#ok<*VUNUS>
 %radar.stranglemeasureMode{radar.nanglemeasureMode};

%  anglemeasureMode_antanna
 %radar.func=str2num(get(handles.func_antanna,'string'));
 
save  data/antannaParameter  radar;
 h1=msgbox('保存参数成功');
 pause(1);
close (h1);


% --- Executes on button press in reset.

function Gt_antanna_Callback(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to Gt_antanna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Gt_antanna as text
%        str2double(get(hObject,'String')) returns contents of Gt_antanna as a double


% --- Executes during object creation, after setting all properties.
function Gt_antanna_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gt_antanna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Gr_antanna_Callback(hObject, eventdata, handles)
% hObject    handle to Gr_antanna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Gr_antanna as text
%        str2double(get(hObject,'String')) returns contents of Gr_antanna as a double


% --- Executes during object creation, after setting all properties.
function Gr_antanna_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gr_antanna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in func_antanna.
function func_antanna_Callback(hObject, eventdata, handles)
% hObject    handle to func_antanna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns func_antanna contents as cell array
%        contents{get(hObject,'Value')} returns selected item from func_antanna

% --- Executes during object creation, after setting all properties.
function func_antanna_CreateFcn(hObject, eventdata, handles)
% hObject    handle to func_antanna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x_antanna_Callback(hObject, eventdata, handles)
% hObject    handle to x_antanna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_antanna as text
%        str2double(get(hObject,'String')) returns contents of x_antanna as a double


% --- Executes during object creation, after setting all properties.
function x_antanna_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_antanna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y_antanna_Callback(hObject, eventdata, handles)
% hObject    handle to y_antanna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_antanna as text
%        str2double(get(hObject,'String')) returns contents of y_antanna as a double


% --- Executes during object creation, after setting all properties.
function y_antanna_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_antanna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function z_antanna_Callback(hObject, eventdata, handles)
% hObject    handle to z_antanna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of z_antanna as text
%        str2double(get(hObject,'String')) returns contents of z_antanna as a double


% --- Executes during object creation, after setting all properties.
function z_antanna_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z_antanna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function AzAngle_antanna_Callback(hObject, eventdata, handles)
% hObject    handle to AzAngle_antanna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AzAngle_antanna as text
%        str2double(get(hObject,'String')) returns contents of AzAngle_antanna as a double


% --- Executes during object creation, after setting all properties.
function AzAngle_antanna_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AzAngle_antanna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ELAngle_antanna_Callback(hObject, eventdata, handles)
% hObject    handle to ELAngle_antanna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ELAngle_antanna as text
%        str2double(get(hObject,'String')) returns contents of ELAngle_antanna as a double


% --- Executes during object creation, after setting all properties.
function ELAngle_antanna_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ELAngle_antanna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.




% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 set(handles.x_antanna,'string',' ');
set(handles.y_antanna,'string',' ');
 set(handles.z_antanna,'string',' ');
 set(handles.Gt_antanna,'string',' ');
 set(handles.Gr_antanna,'string',' ');
set(handles.AzAngle_antanna,'string',' ');
set(handles.ELAngle_antanna,'string',' ');
set(handles.AdjacentAngle_antanna,'string',' ');
set(handles.beamwidth_antanna,'string',' ');

 h2=msgbox('重置参数成功');
 pause(1);
close (h2);


% --- Executes on selection change in anglemeasureMode_antanna.
function anglemeasureMode_antanna_Callback(hObject, eventdata, handles)
% hObject    handle to anglemeasureMode_antanna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns anglemeasureMode_antanna contents as cell array
%        contents{get(hObject,'Value')} returns selected item from anglemeasureMode_antanna


% --- Executes during object creation, after setting all properties.
function anglemeasureMode_antanna_CreateFcn(hObject, eventdata, handles)
% hObject    handle to anglemeasureMode_antanna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function beamwidth_antanna_Callback(hObject, eventdata, handles)
% hObject    handle to beamwidth_antanna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of beamwidth_antanna as text
%        str2double(get(hObject,'String')) returns contents of beamwidth_antanna as a double


% --- Executes during object creation, after setting all properties.
function beamwidth_antanna_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beamwidth_antanna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function AdjacentAngle_antanna_Callback(hObject, eventdata, handles)
% hObject    handle to AdjacentAngle_antanna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AdjacentAngle_antanna as text
%        str2double(get(hObject,'String')) returns contents of AdjacentAngle_antanna as a double


% --- Executes during object creation, after setting all properties.
function AdjacentAngle_antanna_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AdjacentAngle_antanna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% --- Executes on selection change in func_antenna.
function func_antenna_Callback(hObject, eventdata, handles)
% hObject    handle to func_antenna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns func_antenna contents as cell array
%        contents{get(hObject,'Value')} returns selected item from func_antenna


% --- Executes during object creation, after setting all properties.
function func_antenna_CreateFcn(hObject, eventdata, handles)
% hObject    handle to func_antenna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to anglemeasureMode_antanna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of anglemeasureMode_antanna as text
%        str2double(get(hObject,'String')) returns contents of anglemeasureMode_antanna as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to anglemeasureMode_antanna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
