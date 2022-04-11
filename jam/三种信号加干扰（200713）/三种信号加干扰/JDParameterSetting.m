function varargout = JDParameterSetting(varargin)
% JDPARAMETERSETTING MATLAB code for JDParameterSetting.fig
%      JDPARAMETERSETTING, by itself, creates a new JDPARAMETERSETTING or raises the existing
%      singleton*.
%
%      H = JDPARAMETERSETTING returns the handle to a new JDPARAMETERSETTING or the handle to
%      the existing singleton*.
%
%      JDPARAMETERSETTING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in JDPARAMETERSETTING.M with the given input arguments.
%
%      JDPARAMETERSETTING('Property','Value',...) creates a new JDPARAMETERSETTING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before JDParameterSetting_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to JDParameterSetting_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help JDParameterSetting

% Last Modified by GUIDE v2.5 07-Aug-2016 20:03:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @JDParameterSetting_OpeningFcn, ...
                   'gui_OutputFcn',  @JDParameterSetting_OutputFcn, ...
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


% --- Executes just before JDParameterSetting is made visible.
function JDParameterSetting_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to JDParameterSetting (see VARARGIN)

% Choose default command line output for JDParameterSetting
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes JDParameterSetting wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = JDParameterSetting_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function B_JD_Callback(hObject, eventdata, handles) %#ok<*DEFNU,*INUSD>
% hObject    handle to B_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of B_JD as text
%        str2double(get(hObject,'String')) returns contents of B_JD as a double


% --- Executes during object creation, after setting all properties.
function B_JD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Gt_JD_Callback(hObject, eventdata, handles)
% hObject    handle to Gt_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Gt_JD as text
%        str2double(get(hObject,'String')) returns contents of Gt_JD as a double


% --- Executes during object creation, after setting all properties.
function Gt_JD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gt_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Gr_JD_Callback(hObject, eventdata, handles)
% hObject    handle to Gr_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Gr_JD as text
%        str2double(get(hObject,'String')) returns contents of Gr_JD as a double


% --- Executes during object creation, after setting all properties.
function Gr_JD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gr_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fz_JD_Callback(hObject, eventdata, handles)
% hObject    handle to fz_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fz_JD as text
%        str2double(get(hObject,'String')) returns contents of fz_JD as a double


% --- Executes during object creation, after setting all properties.
function fz_JD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fz_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function F_JD_Callback(hObject, eventdata, handles)
% hObject    handle to F_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of F_JD as text
%        str2double(get(hObject,'String')) returns contents of F_JD as a double


% --- Executes during object creation, after setting all properties.
function F_JD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to F_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pt_JD_Callback(hObject, eventdata, handles)
% hObject    handle to Pt_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pt_JD as text
%        str2double(get(hObject,'String')) returns contents of Pt_JD as a double


% --- Executes during object creation, after setting all properties.
function Pt_JD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pt_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function L_JD_Callback(hObject, eventdata, handles)
% hObject    handle to L_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of L_JD as text
%        str2double(get(hObject,'String')) returns contents of L_JD as a double


% --- Executes during object creation, after setting all properties.
function L_JD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to L_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tau_JD_Callback(hObject, eventdata, handles)
% hObject    handle to tau_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tau_JD as text
%        str2double(get(hObject,'String')) returns contents of tau_JD as a double


% --- Executes during object creation, after setting all properties.
function tau_JD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tau_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tr_JD_Callback(hObject, eventdata, handles)
% hObject    handle to tr_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tr_JD as text
%        str2double(get(hObject,'String')) returns contents of tr_JD as a double


% --- Executes during object creation, after setting all properties.
function tr_JD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tr_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f0_JD_Callback(hObject, eventdata, handles)
% hObject    handle to f0_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f0_JD as text
%        str2double(get(hObject,'String')) returns contents of f0_JD as a double


% --- Executes during object creation, after setting all properties.
function f0_JD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f0_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fs_JD_Callback(hObject, eventdata, handles)
% hObject    handle to fs_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fs_JD as text
%        str2double(get(hObject,'String')) returns contents of fs_JD as a double


% --- Executes during object creation, after setting all properties.
function fs_JD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frame_JD_Callback(hObject, eventdata, handles)
% hObject    handle to frame_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_JD as text
%        str2double(get(hObject,'String')) returns contents of frame_JD as a double


% --- Executes during object creation, after setting all properties.
function frame_JD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to Gt_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Gt_JD as text
%        str2double(get(hObject,'String')) returns contents of Gt_JD as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gt_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to Gr_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Gr_JD as text
%        str2double(get(hObject,'String')) returns contents of Gr_JD as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gr_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to fz_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fz_JD as text
%        str2double(get(hObject,'String')) returns contents of fz_JD as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fz_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to F_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of F_JD as text
%        str2double(get(hObject,'String')) returns contents of F_JD as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to F_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to Pt_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pt_JD as text
%        str2double(get(hObject,'String')) returns contents of Pt_JD as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pt_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to L_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of L_JD as text
%        str2double(get(hObject,'String')) returns contents of L_JD as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to L_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to tau_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tau_JD as text
%        str2double(get(hObject,'String')) returns contents of tau_JD as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tau_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to tr_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tr_JD as text
%        str2double(get(hObject,'String')) returns contents of tr_JD as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tr_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to f0_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f0_JD as text
%        str2double(get(hObject,'String')) returns contents of f0_JD as a double


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f0_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to fs_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fs_JD as text
%        str2double(get(hObject,'String')) returns contents of fs_JD as a double


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to B_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of B_JD as text
%        str2double(get(hObject,'String')) returns contents of B_JD as a double


% --- Executes during object creation, after setting all properties.
function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Te_JD_Callback(hObject, eventdata, handles)
% hObject    handle to Te_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Te_JD as text
%        str2double(get(hObject,'String')) returns contents of Te_JD as a double


% --- Executes during object creation, after setting all properties.
function Te_JD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Te_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f1_JD_Callback(hObject, eventdata, handles)
% hObject    handle to f1_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f1_JD as text
%        str2double(get(hObject,'String')) returns contents of f1_JD as a double


% --- Executes during object creation, after setting all properties.
function f1_JD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f1_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_jilei_JD_Callback(hObject, eventdata, handles)
% hObject    handle to num_jilei_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_jilei_JD as text
%        str2double(get(hObject,'String')) returns contents of num_jilei_JD as a double


% --- Executes during object creation, after setting all properties.
function num_jilei_JD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_jilei_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_tongdao_JD_Callback(hObject, eventdata, handles)
% hObject    handle to num_tongdao_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_tongdao_JD as text
%        str2double(get(hObject,'String')) returns contents of num_tongdao_JD as a double


% --- Executes during object creation, after setting all properties.
function num_tongdao_JD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_tongdao_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_cankao_JD_Callback(hObject, eventdata, handles)
% hObject    handle to num_cankao_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_cankao_JD as text
%        str2double(get(hObject,'String')) returns contents of num_cankao_JD as a double


% --- Executes during object creation, after setting all properties.
function num_cankao_JD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_cankao_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_baohu_JD_Callback(hObject, eventdata, handles)
% hObject    handle to num_baohu_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_baohu_JD as text
%        str2double(get(hObject,'String')) returns contents of num_baohu_JD as a double


% --- Executes during object creation, after setting all properties.
function num_baohu_JD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_baohu_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pfa_JD_Callback(hObject, eventdata, handles)
% hObject    handle to Pfa_JD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pfa_JD as text
%        str2double(get(hObject,'String')) returns contents of Pfa_JD as a double


% --- Executes during object creation, after setting all properties.
function Pfa_JD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pfa_JD (see GCBO)
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
B_JD=str2double(get(handles.B_JD,'string'));  
Gt_JD=str2double(get(handles.Gt_JD,'string'));   
Gr_JD=str2double(get(handles.Gr_JD,'string'));   
fz_JD=str2double(get(handles.fz_JD,'string'));   
F_JD=str2double(get(handles.F_JD,'string'));   
Pt_JD=str2double(get(handles.Pt_JD,'string')); 
L_JD=str2double(get(handles.L_JD,'string')); 
Te_JD=str2double(get(handles.Te_JD,'string')); 
tau_JD=str2double(get(handles.tau_JD,'string')); 
tr_JD=str2double(get(handles.tr_JD,'string')); 
f0_JD=str2double(get(handles.f0_JD,'string')); 
fs_JD=str2double(get(handles.fs_JD,'string')); 
frame_JD=str2double(get(handles.frame_JD,'string')); 
f1_JD=str2double(get(handles.f1_JD,'string')); 
num_jilei_JD=str2double(get(handles.num_jilei_JD,'string'));
num_tongdao_JD=str2double(get(handles.num_tongdao_JD,'string'));
num_cankao_JD=str2double(get(handles.num_cankao_JD,'string'));
num_baohu_JD=str2double(get(handles.num_baohu_JD,'string'));
Pfa_JD=str2double(get(handles.Pfa_JD,'string'));
save data/data_JDParameter B_JD Gt_JD Gr_JD  fz_JD  F_JD Pt_JD L_JD Te_JD tau_JD  tr_JD f0_JD fs_JD  frame_JD  f1_JD num_jilei_JD num_tongdao_JD num_cankao_JD num_baohu_JD Pfa_JD
% R=str2double(get(handles.R,'string')); 
% v=str2double(get(handles.v,'string')); 
% sigma0=str2double(get(handles.sigma0,'string')); 
% rcsk=str2double(get(handles.rcsk,'string')); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%生成发射信号
ts=1/fs_JD;
tm=0:1/fs_JD:tr_JD-1/fs_JD;  
N=length(tm);
[y]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(单位：秒)'), ylabel('y(单位：伏)'),title('简单脉冲发射信号');
figure,plot((0:fs_JD/length(abs(fft(abs(fftshift(y))))):fs_JD-fs_JD/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('简单脉冲发射信号的频谱');




% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load data/data_JDParameter;
set(handles.B_JD,'string',B_JD);
set(handles.Gt_JD,'string',Gt_JD);
set(handles.Gr_JD,'string',Gr_JD);
set(handles.fz_JD,'string',fz_JD);
set(handles.F_JD,'string',F_JD);
set(handles.Pt_JD,'string',Pt_JD);
set(handles.L_JD,'string',L_JD);
set(handles.tau_JD,'string',tau_JD);
set(handles.tr_JD,'string',tr_JD);
set(handles.f0_JD,'string',f0_JD);
set(handles.fs_JD,'string',fs_JD);
set(handles.frame_JD,'string',frame_JD);
set(handles.Te_JD,'string',Te_JD);
set(handles.f1_JD,'string',f1_JD);
set(handles.num_jilei_JD,'string',num_jilei_JD);
set(handles.num_tongdao_JD,'string',num_tongdao_JD);
set(handles.num_cankao_JD,'string',num_cankao_JD);
set(handles.num_baohu_JD,'string',num_cankao_JD);
set(handles.Pfa_JD,'string',num_cankao_JD);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.B_JD,'string','');
set(handles.Gt_JD,'string','');
set(handles.Gr_JD,'string','');
set(handles.fz_JD,'string','');
set(handles.F_JD,'string','');
set(handles.Pt_JD,'string','');
set(handles.L_JD,'string','');
set(handles.tau_JD,'string','');
set(handles.tr_JD,'string','');
set(handles.f0_JD,'string','');
set(handles.fs_JD,'string','');
set(handles.frame_JD,'string','');
set(handles.Te_JD,'string','');
set(handles.f1_JD,'string','');
set(handles.num_jilei_JD,'string','');
set(handles.num_tongdao_JD,'string','');
set(handles.num_cankao_JD,'string','');
set(handles.num_baohu_JD,'string','');
set(handles.Pfa_JD,'string','');
