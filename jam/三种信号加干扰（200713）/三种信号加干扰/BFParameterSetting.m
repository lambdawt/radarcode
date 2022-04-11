function varargout = BFParameterSetting(varargin)
% BFPARAMETERSETTING MATLAB code for BFParameterSetting.fig
%      BFPARAMETERSETTING, by itself, creates a new BFPARAMETERSETTING or raises the existing
%      singleton*.
%
%      H = BFPARAMETERSETTING returns the handle to a new BFPARAMETERSETTING or the handle to
%      the existing singleton*.
%
%      BFPARAMETERSETTING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BFPARAMETERSETTING.M with the given input arguments.
%
%      BFPARAMETERSETTING('Property','Value',...) creates a new BFPARAMETERSETTING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BFParameterSetting_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BFParameterSetting_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BFParameterSetting

% Last Modified by GUIDE v2.5 15-Jun-2020 15:07:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BFParameterSetting_OpeningFcn, ...
                   'gui_OutputFcn',  @BFParameterSetting_OutputFcn, ...
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


% --- Executes just before BFParameterSetting is made visible.
function BFParameterSetting_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BFParameterSetting (see VARARGIN)

% Choose default command line output for BFParameterSetting
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BFParameterSetting wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BFParameterSetting_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in maxingxuanze.
function maxingxuanze_Callback(hObject, eventdata, handles)
% hObject    handle to maxingxuanze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global flag;
n=get(hObject,'value');
% str=get(hObject,'string');
if n==3
    flag = 3;
    set(handles.maxulie_B,'Enable','on');
    
%      set(handles.text5,'Enable','on'),flag=3;
else
    if n==2
   flag=2;
  
    else
     flag=1; 
     
    end
    set(handles.maxulie_B,'Enable','off');
%     set(handles.text5,'Enable','off')
     set(handles.maxulie_B,'string','');
end

% Hints: contents = cellstr(get(hObject,'String')) returns maxingxuanze contents as cell array
%        contents{get(hObject,'Value')} returns selected item from maxingxuanze


% --- Executes during object creation, after setting all properties.
function maxingxuanze_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxingxuanze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
global flag;
flag=1;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tau_B_Callback(hObject, eventdata, handles) %#ok<*INUSD,*DEFNU>
% hObject    handle to tau_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tau_B as text
%        str2double(get(hObject,'String')) returns contents of tau_B as a double


% --- Executes during object creation, after setting all properties.
function tau_B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tau_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fs_B_Callback(hObject, eventdata, handles)
% hObject    handle to fs_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fs_B as text
%        str2double(get(hObject,'String')) returns contents of fs_B as a double


% --- Executes during object creation, after setting all properties.
function fs_B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f0_B_Callback(hObject, eventdata, handles)
% hObject    handle to f0_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f0_B as text
%        str2double(get(hObject,'String')) returns contents of f0_B as a double


% --- Executes during object creation, after setting all properties.
function f0_B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f0_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxulie_B_Callback(hObject, eventdata, handles)
% hObject    handle to maxulie_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxulie_B as text
%        str2double(get(hObject,'String')) returns contents of maxulie_B as a double


% --- Executes during object creation, after setting all properties.
function maxulie_B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxulie_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% n=get(handles.maxingxuanze,'value');
% str=get(handles.maxingxuanze,'string');
% 
% if str{n}=='自定义'
%    set(hObject,'Visible','on');
% else
    set(hObject,'Enable','off');
%     
% end
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tr_B_Callback(hObject, eventdata, handles)
% hObject    handle to tr_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tr_B as text
%        str2double(get(hObject,'String')) returns contents of tr_B as a double


% --- Executes during object creation, after setting all properties.
function tr_B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tr_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pt_Callback(hObject, eventdata, handles)
% hObject    handle to Pt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pt as text
%        str2double(get(hObject,'String')) returns contents of Pt as a double


% --- Executes during object creation, after setting all properties.
function Pt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frame_B_Callback(hObject, eventdata, handles)
% hObject    handle to frame_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_B as text
%        str2double(get(hObject,'String')) returns contents of frame_B as a double


% --- Executes during object creation, after setting all properties.
function frame_B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_B (see GCBO)
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

global code;
tau_B=str2num(get(handles.tau_B,'string'));
fs_B=str2num(get(handles.fs_B,'string'));
f0_B=str2num(get(handles.f0_B,'string'));
tr_B=str2num(get(handles.tr_B,'string'));
frame_B=str2num(get(handles.frame_B,'string'));
Te_B=str2num(get(handles.Te_B,'string'));
f1_B=str2num(get(handles.f1_B,'string'));
num_jilei_B=str2num(get(handles.num_jilei_B,'string'));
num_tongdao_B=str2num(get(handles.num_tongdao_B,'string'));
num_cankao_B=str2num(get(handles.num_cankao_B,'string'));
num_baohu_B=str2num(get(handles.num_baohu_B,'string'));
Pfa_B=str2num(get(handles.Pfa_B,'string'));
L_B=str2num(get(handles.L_B,'string'));

% % % str_B=get(handles.maxulie_B,'string');
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%将输入框中的字符串转换为数字，存储于数组code中
% % % len=length(str_B);
% % % code1=zeros(1,len);  
% % % flag2=0;
% % % i=1;
% % % while i<=len
% % %     if str_B(i)=='-'
% % %          code1(i)=-1;
% % %          i=i+1;
% % %          flag2=1;
% % %     else
% % %        code1(i)=str2num(str_B(i));
% % %     end
% % %      i=i+1;
% % % end
% % % 
% % % if flag2==1
% % %     [i,j,code]=find(code1);
% % %     code=code,
% % % else if flag==1
% % %         code=[1,1,1,-1,-1,-1,1],%^7位巴克码
% % %     else
% % %         code=[1,1,1,1,1,-1,-1,1,1,-1,1,-1,1],%13位巴克码 
% % %     end
% % % end

% if flag==1
%     code=[1,1,1,-1,-1,-1,1];%^7位巴克码
% else if flag==2
%         code=[1,1,1,1,1,-1,-1,1,1,-1,1,-1,1];%13位巴克码 
%     else
%         
%         str_B=get(handles.maxulie_B,'string');
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%将输入框中的字符串转换为数字，存储于数组code中
%         len=length(str_B);
%         code1=zeros(1,len);  
%         flag2=0;
%         i=1;
%         while i<=len
%             if str_B(i)=='-'
%                  code1(i)=-1;
%                  i=i+1;
%                  flag2=1;
%             else
%                code1(i)=str2num(str_B(i));
%             end
%              i=i+1;
%         end
%         if flag2==1
%             [i,j,code]=find(code1);
%             code=code;
%         else
%           code=code1;  
%         end
%     end
% end

global flag;

 if flag==1
    code=[1,1,1,-1,-1,-1,1];%^7位巴克码
 else if flag==2
        code=[1,1,1,1,1,-1,-1,1,1,-1,1,-1,1];%13位巴克码 
     else  
        str_B=get(handles.maxulie_B,'string');
        code=str2num(str_B);
     end
end       








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c=3e8;
F_B=str2num(get(handles.F_B,'string'));
B_B=str2num(get(handles.B_B,'string'));
fz_B=str2num(get(handles.fz_B,'string'));
Gt_B=str2num(get(handles.Gt_B,'string'));
Gr_B=str2num(get(handles.Gr_B,'string'));
Pt_B=str2num(get(handles.Pt_B,'string'));
save  data/data_BFParameter tau_B fs_B f0_B tr_B frame_B Te_B f1_B num_jilei_B num_tongdao_B num_cankao_B num_baohu_B Pfa_B L_B F_B B_B fz_B Gt_B Gr_B Pt_B code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tm_B=0:1/fs_B:tr_B-1/fs_B;%一个脉冲重复周期采样序列
N=length(tm_B);%一个脉冲重复周期采样点数长度
number1=length(code);
ts=1/fs_B;
[y]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(单位：秒)'), ylabel('y(单位：伏)'),title('相位编码发射信号');
figure,plot((0:fs_B/length(abs(fft(abs(fftshift(y))))):fs_B-fs_B/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('频率f(单位：Hz)'), ylabel('y(单位：伏)'),title('相位编码发射信号的频谱');





function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.tau_B,'string','');
set(handles.fs_B,'string','');
set(handles.f0_B,'string','');
set(handles.tr_B,'string','');
set(handles.frame_B,'string','');
set(handles.maxulie_B,'string','');
set(handles.F_B,'string','');
set(handles.B_B,'string','');
set(handles.fz_B,'string','');
set(handles.Gt_B,'string','');
set(handles.Pt_B,'string','');
set(handles.Te_B,'string','');
set(handles.f1_B,'string','');
set(handles.L_B,'string','');
set(handles.num_jilei_B,'string','');
set(handles.num_tongdao_B,'string','');
set(handles.num_cankao_B,'string','');
set(handles.num_baohu_B,'string','');
set(handles.Pfa_B,'string','');
set(handles.num_jilei_B,'string','');
set(handles.Gr_B,'string','');
% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load data/data_BFParameter
set(handles.tau_B,'string',tau_B);
set(handles.fs_B,'string',fs_B);
set(handles.f0_B,'string',f0_B);
set(handles.tr_B,'string',tr_B);
set(handles.frame_B,'string',frame_B);
set(handles.maxulie_B,'Enable','on');
maxulie_B=num2str(code);
set(handles.maxulie_B,'string',maxulie_B);
set(handles.F_B,'string',F_B);
set(handles.B_B,'string',B_B);
set(handles.fz_B,'string',fz_B);
set(handles.Gt_B,'string',Gt_B);
set(handles.Pt_B,'string',Pt_B);
set(handles.Te_B,'string',Te_B);
set(handles.f1_B,'string',f1_B);
set(handles.L_B,'string',L_B);
set(handles.num_jilei_B,'string',num_jilei_B);
set(handles.num_tongdao_B,'string',num_tongdao_B);
set(handles.num_baohu_B,'string',num_baohu_B);
set(handles.num_cankao_B,'string',num_cankao_B);
set(handles.Pfa_B,'string',Pfa_B);
set(handles.Gr_B,'string',Gr_B);
% --- Executes during object creation, after setting all properties.
function text5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
%  set(hObject,'Enable','off');
   
    
    
    



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function F_B_Callback(hObject, eventdata, handles)
% hObject    handle to F_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of F_B as text
%        str2double(get(hObject,'String')) returns contents of F_B as a double


% --- Executes during object creation, after setting all properties.
function F_B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to F_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function B_B_Callback(hObject, eventdata, handles)
% hObject    handle to B_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of B_B as text
%        str2double(get(hObject,'String')) returns contents of B_B as a double


% --- Executes during object creation, after setting all properties.
function B_B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Rmax_B_Callback(hObject, eventdata, handles)
% hObject    handle to Rmax_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Rmax_B as text
%        str2double(get(hObject,'String')) returns contents of Rmax_B as a double


% --- Executes during object creation, after setting all properties.
function Rmax_B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rmax_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigma_B_Callback(hObject, eventdata, handles)
% hObject    handle to sigma_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigma_B as text
%        str2double(get(hObject,'String')) returns contents of sigma_B as a double


% --- Executes during object creation, after setting all properties.
function sigma_B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigma_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fz_B_Callback(hObject, eventdata, handles)
% hObject    handle to fz_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fz_B as text
%        str2double(get(hObject,'String')) returns contents of fz_B as a double


% --- Executes during object creation, after setting all properties.
function fz_B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fz_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Gt_B_Callback(hObject, eventdata, handles)
% hObject    handle to Gt_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Gt_B as text
%        str2double(get(hObject,'String')) returns contents of Gt_B as a double


% --- Executes during object creation, after setting all properties.
function Gt_B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gt_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pt_B_Callback(hObject, eventdata, handles)
% hObject    handle to Pt_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pt_B as text
%        str2double(get(hObject,'String')) returns contents of Pt_B as a double


% --- Executes during object creation, after setting all properties.
function Pt_B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pt_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function L_B_Callback(hObject, eventdata, handles)
% hObject    handle to L_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of L_B as text
%        str2double(get(hObject,'String')) returns contents of L_B as a double


% --- Executes during object creation, after setting all properties.
function L_B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to L_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Te_B_Callback(hObject, eventdata, handles)
% hObject    handle to Te_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Te_B as text
%        str2double(get(hObject,'String')) returns contents of Te_B as a double


% --- Executes during object creation, after setting all properties.
function Te_B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Te_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f1_B_Callback(hObject, eventdata, handles)
% hObject    handle to f1_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of f1_B as text
%        str2double(get(hObject,'String')) returns contents of f1_B as a double


% --- Executes during object creation, after setting all properties.
function f1_B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f1_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_jilei_B_Callback(hObject, eventdata, handles)
% hObject    handle to num_jilei_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_jilei_B as text
%        str2double(get(hObject,'String')) returns contents of num_jilei_B as a double


% --- Executes during object creation, after setting all properties.
function num_jilei_B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_jilei_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_cankao_B_Callback(hObject, eventdata, handles)
% hObject    handle to num_cankao_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_cankao_B as text
%        str2double(get(hObject,'String')) returns contents of num_cankao_B as a double


% --- Executes during object creation, after setting all properties.
function num_cankao_B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_cankao_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_baohu_B_Callback(hObject, eventdata, handles)
% hObject    handle to num_baohu_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_baohu_B as text
%        str2double(get(hObject,'String')) returns contents of num_baohu_B as a double


% --- Executes during object creation, after setting all properties.
function num_baohu_B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_baohu_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pfa_B_Callback(hObject, eventdata, handles)
% hObject    handle to Pfa_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pfa_B as text
%        str2double(get(hObject,'String')) returns contents of Pfa_B as a double


% --- Executes during object creation, after setting all properties.
function Pfa_B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pfa_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_tongdao_B_Callback(hObject, eventdata, handles)
% hObject    handle to num_tongdao_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_tongdao_B as text
%        str2double(get(hObject,'String')) returns contents of num_tongdao_B as a double


% --- Executes during object creation, after setting all properties.
function num_tongdao_B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_tongdao_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Gr_B_Callback(hObject, eventdata, handles)
% hObject    handle to Gr_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Gr_B as text
%        str2double(get(hObject,'String')) returns contents of Gr_B as a double


% --- Executes during object creation, after setting all properties.
function Gr_B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gr_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to L_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of L_B as text
%        str2double(get(hObject,'String')) returns contents of L_B as a double


% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to L_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on maxulie_B and none of its controls.
function maxulie_B_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to maxulie_B (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
