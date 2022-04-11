function varargout = other(varargin)
% OTHER MATLAB code for other.fig
%      OTHER, by itself, creates a new OTHER or raises the existing
%      singleton*.
%
%      H = OTHER returns the handle to a new OTHER or the handle to
%      the existing singleton*.
%
%      OTHER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OTHER.M with the given input arguments.
%
%      OTHER('Property','Value',...) creates a new OTHER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before other_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to other_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help other

% Last Modified by GUIDE v2.5 29-May-2020 15:17:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @other_OpeningFcn, ...
                   'gui_OutputFcn',  @other_OutputFcn, ...
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


% --- Executes just before other is made visible.
function other_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to other (see VARARGIN)

% Choose default command line output for other
handles.output = hObject;

global sigma0;
sigma0 = 2;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes other wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = other_OutputFcn(hObject, eventdata, handles) 
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

str=get(handles.JammingS,'string');
n=get(handles.JammingS,'Value');
switch str{n}

case  '多普勒闪烁干扰'
set(handles.fd_dopplerblink,'string',' ');
set(handles.Td_dopplerblink,'string',' ');
set(handles.Pj_dopplerblink,'string',' ');
set(handles.fs_dopplerblink,'string',' ');
set(handles.flagT_dopplerblink,'string',' ');

case '箔条干扰'
  
set(handles.tf_botiao,'string',' ');
set(handles.sf_botiao,'string',' ');
set(handles.vl_botiao,'string',' ');
set(handles.vf_botiao,'string',' ');
set(handles.ts_botiao,'string',' ');
set(handles.bt_botiao,'string','  ');
set(handles.al_botiao,'string',' ');
set(handles.sref_botiao,'string',' ');
set(handles.smax_botiao,'string',' ');

case 'AGC干扰'

  set(handles.fs_AGC,'string',' ');
set(handles.CurrentT_AGC,'string',' ');
set(handles.Pj_AGC,'string',' ');
set(handles.Period_AGC,'string',' ');
set(handles.D_AGC,'string',' ');
set(handles.radio_AGC,'string',' ');
end
% h3=msgbox('初始化目标参数成功');
% pause(1);
% close (h3);

% --- Executes on selection change in JammingS.
function JammingS_Callback(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to JammingS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns JammingS contents as cell array
%        contents{get(hObject,'Value')} returns selected item from JammingS


% --- Executes during object creation, after setting all properties.
function JammingS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to JammingS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

strOJS=get(handles.JammingS,'string');
nOJS=get(handles.JammingS,'Value');
switch strOJS{nOJS}
    
   case  '多普勒闪烁干扰'
                  
fd_dopplerblink=str2double(get(handles.fd_dopplerblink,'string'));
Td_dopplerblink=str2double(get(handles.Td_dopplerblink,'string'));
Pj_dopplerblink=str2double(get(handles.Pj_dopplerblink,'string'));
fs_dopplerblink=str2double(get(handles.fs_dopplerblink,'string'));
flagT_dopplerblink=str2double(get(handles.flagT_dopplerblink,'string'));
save data/data_dopplerblink fd_dopplerblink Td_dopplerblink Pj_dopplerblink fs_dopplerblink flagT_dopplerblink

case '箔条干扰'
  
% set(handles.tf_botiao,'string','0.99');
% set(handles.sf_botiao,'string','50');
% set(handles.vl_botiao,'string','50');
% set(handles.vf_botiao,'string','25');
% set(handles.ts_botiao,'string','0.5');
% set(handles.bt_botiao,'string','30*pi*180');
% set(handles.al_botiao,'string','30*pi*180');
% set(handles.sref_botiao,'string','1');
% set(handles.smax_botiao,'string','1');

tf_botiao=str2double(get(handles.tf_botiao,'string'));
sf_botiao=str2double(get(handles.sf_botiao,'string'));
vl_botiao=str2double(get(handles.vl_botiao,'string'));
vf_botiao=str2double(get(handles.vf_botiao,'string'));
ts_botiao=str2double(get(handles.ts_botiao,'string'));
bt_botiao=str2double(get(handles.bt_botiao,'string'));
al_botiao=str2double(get(handles.al_botiao,'string'));
sref_botiao=str2double(get(handles.sref_botiao,'string'));
smax_botiao=str2double(get(handles.smax_botiao,'string'));
save data/data_botiao tf_botiao sf_botiao vl_botiao vf_botiao ts_botiao bt_botiao al_botiao sref_botiao smax_botiao

case 'AGC干扰'
% set(handles.fs_AGC,'string',' ');
% set(handles.CurrentT_AGC,'string',' ');
% set(handles.Pj_AGC,'string','16');
% set(handles.Period_AGC,'string',' ');
% set(handles.D_AGC,'string',' ');
% set(handles.radio_AGC,'string',' ');
fs_AGC=str2double(get(handles.fs_AGC,'string'));
CurrentT_AGC=str2double(get(handles.CurrentT_AGC,'string'));
Pj_AGC=str2double(get(handles.Pj_AGC,'string'));
Period_AGC=str2double(get(handles.Period_AGC,'string'));
D_AGC=str2double(get(handles.D_AGC,'string'));
radio_AGC=str2double(get(handles.radio_AGC,'string'));
save data/data_AGC fs_AGC CurrentT_AGC Pj_AGC Period_AGC D_AGC radio_AGC
end

h1=msgbox('载入目标参数成功');
pause(1);
close(h1);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str=get(handles.JammingS,'string');
n=get(handles.JammingS,'Value');
switch str{n}
    
case  '多普勒闪烁干扰'
set(handles.fd_dopplerblink,'string','0.6e4');
set(handles.Td_dopplerblink,'string','100e-6');
set(handles.Pj_dopplerblink,'string','4');
set(handles.fs_dopplerblink,'string','2e6');
set(handles.flagT_dopplerblink,'string','800e-6');

case '箔条干扰'
  
set(handles.tf_botiao,'string','0.99');
set(handles.sf_botiao,'string','50');
set(handles.vl_botiao,'string','50');
set(handles.vf_botiao,'string','25');
set(handles.ts_botiao,'string','0.5');
set(handles.bt_botiao,'string','30');
set(handles.al_botiao,'string','30');
set(handles.sref_botiao,'string','1');
set(handles.smax_botiao,'string','1');
case 'AGC干扰'
set(handles.fs_AGC,'string','4e6');
set(handles.CurrentT_AGC,'string','900e-6');
set(handles.Pj_AGC,'string','16');
set(handles.Period_AGC,'string','600e-6');
set(handles.D_AGC,'string','0.5');
set(handles.radio_AGC,'string','0.25');
end
h3=msgbox('初始化目标参数成功');
pause(1);
close (h3);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sigma;
global sigma0;
global rcsk;
sigma = rcs(rcsk,sigma0);
global strO;
global nO;
strO=get(handles.JammingS,'string');
nO=get(handles.JammingS,'Value');
switch strO{nO}
 
 case '多普勒闪烁干扰'
%             load data/data_tf
            load data/data_JDParameter %#ok<*LOAD>
            load data/data_target0Parameter
            load data/data_dopplerblink
%             global sigma;
%             global rcsk;
%             sigma = rcs(rcsk,sigma0);
            c=3e8;
            lamta=c/fz_JD;
            Prs=(Pt_JD*(10^((Gt_JD*0.1)))*(10^((Gr_JD*0.1)))*(lamta^2)*sigma)/((4*pi)^3*R^4*10^(L_JD*0.1)); %#ok<*NODEF> %目标回波信号功率
            A=sqrt(Prs);
%             ts=1/fs_JD;
            tm=0:1/fs_JD:tr_JD-1/fs_JD;  
            N=length(tm);
            f_doppler=2*v/lamta;
            R0=rand(100);
            [s_echo_2]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD);
            [ sig_jam,t_jam ] = jam_dopplerblink( fd_dopplerblink,Td_dopplerblink,R0,s_echo_2,fs_dopplerblink,Pj_dopplerblink,flagT_dopplerblink ); 
            view_jam_dopplerblink( s_echo_2,sig_jam,t_jam,fs_dopplerblink );
            
    case '箔条干扰'
        load data/data_botiao
        fs=2e6;f0=0.5e6;
        R=2e3;c=3e8;A=5;Tr=600e-6;tau=100e-6;frame=1;tm=0:1/fs:Tr-1/fs;N=length(tm);f_doppler=1e4;
        [s_echo_2]=JDhuiboxinhao(R,c,A,N,frame,fs,f_doppler,tm,f0,tau);

        %目标状态参数
        px=1e3;py=1e3;pz=1e3;%目标位置
        vx=10;vy=10;vz=0;%目标速度
        ax=0;ay=0;az=0;%目标加速度
        phi=pi/180;
        [ TargetStatus ] = paraset_targetstatus( px,py,pz,vx,vy,vz,ax,ay,az,phi );

        [ PassivePara ] = paraset_passivejaming( tf_botiao,sf_botiao,vl_botiao,vf_botiao,ts_botiao,bt_botiao,al_botiao,sref_botiao,smax_botiao );

          vx=1;vy=1;
        [ WindV ] = paraset_windvelocity( vx,vy );


        echo=s_echo_2;fc=f0;CurrentT=2.2;
        [ sig_jam,t_jam ] = jam_passive( echo,fs,fc,CurrentT,TargetStatus,WindV,PassivePara );

        figure;
        subplot(211)
        plot(t_jam,real(sig_jam));grid;xlabel('时间(s)');ylabel('幅度(V)');title('箔条干扰信号时域波形'); 
        subplot(212)
        Fsig=fftshift(fft(sig_jam));
        K=length(t_jam);
        k=floor(-K/2+0.5:K/2-0.5);
        dfs=fs/K;
        plot(k*dfs,abs(Fsig));grid;xlabel('频率(Hz)');ylabel('幅度(V/Hz)');title('箔条干扰信号频域波形');

    case 'AGC干扰'
        fs=2e6;f0=0.5e6;
        R=2e3;c=3e8;A=5;Tr=600e-6;tau=100e-6;frame=1;tm=0:1/fs:Tr-1/fs;N=length(tm);f_doppler=1e4;
        [s_echo_2]=JDhuiboxinhao(R,c,A,N,frame,fs,f_doppler,tm,f0,tau);
        load data/data_AGC
        % Period=Tr;CurrentT=Tr+0.6*Tr;Pj=16;D=0.5;radio=0.25;echo=s_echo_2;
        echo=s_echo_2;
        [ sig_jam,t_jam ] = jam_AGC( CurrentT_AGC,Pj_AGC,Period_AGC,D_AGC,radio_AGC,echo,fs_AGC );

        figure;
        subplot(211)
        plot(t_jam,real(sig_jam));grid;xlabel('时间(s)');ylabel('幅度(V)');title('AGC干扰信号时域波形'); 
        subplot(212)
        Fsig=fftshift(fft(sig_jam));
        K=length(t_jam);
        k=floor(-K/2+0.5:K/2-0.5);
        dfs=fs/K;
        plot(k*dfs,abs(Fsig));grid;xlabel('频率(Hz)');ylabel('幅度(V/Hz)');title('AGC干扰信号频域波形');

end


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
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



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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



function flagT_dopplerblink_Callback(hObject, eventdata, handles)
% hObject    handle to flagT_dopplerblink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of flagT_dopplerblink as text
%        str2double(get(hObject,'String')) returns contents of flagT_dopplerblink as a double


% --- Executes during object creation, after setting all properties.
function flagT_dopplerblink_CreateFcn(hObject, eventdata, handles)
% hObject    handle to flagT_dopplerblink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fd_dopplerblink_Callback(hObject, eventdata, handles)
% hObject    handle to fd_dopplerblink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fd_dopplerblink as text
%        str2double(get(hObject,'String')) returns contents of fd_dopplerblink as a double


% --- Executes during object creation, after setting all properties.
function fd_dopplerblink_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fd_dopplerblink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Td_dopplerblink_Callback(hObject, eventdata, handles)
% hObject    handle to Td_dopplerblink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Td_dopplerblink as text
%        str2double(get(hObject,'String')) returns contents of Td_dopplerblink as a double


% --- Executes during object creation, after setting all properties.
function Td_dopplerblink_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Td_dopplerblink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pj_dopplerblink_Callback(hObject, eventdata, handles)
% hObject    handle to Pj_dopplerblink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pj_dopplerblink as text
%        str2double(get(hObject,'String')) returns contents of Pj_dopplerblink as a double


% --- Executes during object creation, after setting all properties.
function Pj_dopplerblink_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pj_dopplerblink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fs_dopplerblink_Callback(hObject, eventdata, handles)
% hObject    handle to fs_dopplerblink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fs_dopplerblink as text
%        str2double(get(hObject,'String')) returns contents of fs_dopplerblink as a double


% --- Executes during object creation, after setting all properties.
function fs_dopplerblink_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs_dopplerblink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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



function fs_sz_Callback(hObject, eventdata, handles)
% hObject    handle to fs_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fs_sz as text
%        str2double(get(hObject,'String')) returns contents of fs_sz as a double


% --- Executes during object creation, after setting all properties.
function fs_sz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Bj_sz_Callback(hObject, eventdata, handles)
% hObject    handle to Bj_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Bj_sz as text
%        str2double(get(hObject,'String')) returns contents of Bj_sz as a double


% --- Executes during object creation, after setting all properties.
function Bj_sz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bj_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fj_sz_Callback(hObject, eventdata, handles)
% hObject    handle to fj_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fj_sz as text
%        str2double(get(hObject,'String')) returns contents of fj_sz as a double


% --- Executes during object creation, after setting all properties.
function fj_sz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fj_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frame_sz_Callback(hObject, eventdata, handles)
% hObject    handle to frame_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_sz as text
%        str2double(get(hObject,'String')) returns contents of frame_sz as a double


% --- Executes during object creation, after setting all properties.
function frame_sz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Prj_sz_Callback(hObject, eventdata, handles)
% hObject    handle to Prj_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Prj_sz as text
%        str2double(get(hObject,'String')) returns contents of Prj_sz as a double


% --- Executes during object creation, after setting all properties.
function Prj_sz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Prj_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tr_sz_Callback(hObject, eventdata, handles)
% hObject    handle to Tr_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tr_sz as text
%        str2double(get(hObject,'String')) returns contents of Tr_sz as a double


% --- Executes during object creation, after setting all properties.
function Tr_sz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tr_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ns_sz_Callback(hObject, eventdata, handles)
% hObject    handle to Ns_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ns_sz as text
%        str2double(get(hObject,'String')) returns contents of Ns_sz as a double


% --- Executes during object creation, after setting all properties.
function Ns_sz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ns_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T_fr_saopin_Callback(hObject, eventdata, handles)
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



function fs_AGC_Callback(hObject, eventdata, handles)
% hObject    handle to fs_AGC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fs_AGC as text
%        str2double(get(hObject,'String')) returns contents of fs_AGC as a double


% --- Executes during object creation, after setting all properties.
function fs_AGC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs_AGC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CurrentT_AGC_Callback(hObject, eventdata, handles)
% hObject    handle to CurrentT_AGC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CurrentT_AGC as text
%        str2double(get(hObject,'String')) returns contents of CurrentT_AGC as a double


% --- Executes during object creation, after setting all properties.
function CurrentT_AGC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurrentT_AGC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pj_AGC_Callback(hObject, eventdata, handles)
% hObject    handle to Pj_AGC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pj_AGC as text
%        str2double(get(hObject,'String')) returns contents of Pj_AGC as a double


% --- Executes during object creation, after setting all properties.
function Pj_AGC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pj_AGC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Period_AGC_Callback(hObject, eventdata, handles)
% hObject    handle to Period_AGC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Period_AGC as text
%        str2double(get(hObject,'String')) returns contents of Period_AGC as a double


% --- Executes during object creation, after setting all properties.
function Period_AGC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Period_AGC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function D_AGC_Callback(hObject, eventdata, handles)
% hObject    handle to D_AGC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of D_AGC as text
%        str2double(get(hObject,'String')) returns contents of D_AGC as a double


% --- Executes during object creation, after setting all properties.
function D_AGC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to D_AGC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function radio_AGC_Callback(hObject, eventdata, handles)
% hObject    handle to radio_AGC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of radio_AGC as text
%        str2double(get(hObject,'String')) returns contents of radio_AGC as a double


% --- Executes during object creation, after setting all properties.
function radio_AGC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radio_AGC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit46_Callback(hObject, eventdata, handles)
% hObject    handle to edit46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit46 as text
%        str2double(get(hObject,'String')) returns contents of edit46 as a double


% --- Executes during object creation, after setting all properties.
function edit46_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tf_botiao_Callback(hObject, eventdata, handles)
% hObject    handle to tf_botiao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tf_botiao as text
%        str2double(get(hObject,'String')) returns contents of tf_botiao as a double


% --- Executes during object creation, after setting all properties.
function tf_botiao_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tf_botiao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sf_botiao_Callback(hObject, eventdata, handles)
% hObject    handle to sf_botiao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sf_botiao as text
%        str2double(get(hObject,'String')) returns contents of sf_botiao as a double


% --- Executes during object creation, after setting all properties.
function sf_botiao_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sf_botiao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vl_botiao_Callback(hObject, eventdata, handles)
% hObject    handle to vl_botiao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vl_botiao as text
%        str2double(get(hObject,'String')) returns contents of vl_botiao as a double


% --- Executes during object creation, after setting all properties.
function vl_botiao_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vl_botiao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vf_botiao_Callback(hObject, eventdata, handles)
% hObject    handle to vf_botiao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vf_botiao as text
%        str2double(get(hObject,'String')) returns contents of vf_botiao as a double


% --- Executes during object creation, after setting all properties.
function vf_botiao_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vf_botiao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ts_botiao_Callback(hObject, eventdata, handles)
% hObject    handle to ts_botiao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ts_botiao as text
%        str2double(get(hObject,'String')) returns contents of ts_botiao as a double


% --- Executes during object creation, after setting all properties.
function ts_botiao_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ts_botiao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bt_botiao_Callback(hObject, eventdata, handles)
% hObject    handle to bt_botiao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bt_botiao as text
%        str2double(get(hObject,'String')) returns contents of bt_botiao as a double


% --- Executes during object creation, after setting all properties.
function bt_botiao_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bt_botiao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function al_botiao_Callback(hObject, eventdata, handles)
% hObject    handle to al_botiao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of al_botiao as text
%        str2double(get(hObject,'String')) returns contents of al_botiao as a double


% --- Executes during object creation, after setting all properties.
function al_botiao_CreateFcn(hObject, eventdata, handles)
% hObject    handle to al_botiao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function smax_botiao_Callback(hObject, eventdata, handles)
% hObject    handle to smax_botiao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of smax_botiao as text
%        str2double(get(hObject,'String')) returns contents of smax_botiao as a double


% --- Executes during object creation, after setting all properties.
function smax_botiao_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smax_botiao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sref_botiao_Callback(hObject, eventdata, handles)
% hObject    handle to sref_botiao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sref_botiao as text
%        str2double(get(hObject,'String')) returns contents of sref_botiao as a double


% --- Executes during object creation, after setting all properties.
function sref_botiao_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sref_botiao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit49_Callback(hObject, eventdata, handles)
% hObject    handle to edit49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit49 as text
%        str2double(get(hObject,'String')) returns contents of edit49 as a double


% --- Executes during object creation, after setting all properties.
function edit49_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit50_Callback(hObject, eventdata, handles)
% hObject    handle to edit50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit50 as text
%        str2double(get(hObject,'String')) returns contents of edit50 as a double


% --- Executes during object creation, after setting all properties.
function edit50_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit51_Callback(hObject, eventdata, handles)
% hObject    handle to edit51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit51 as text
%        str2double(get(hObject,'String')) returns contents of edit51 as a double


% --- Executes during object creation, after setting all properties.
function edit51_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit52_Callback(hObject, eventdata, handles)
% hObject    handle to edit52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit52 as text
%        str2double(get(hObject,'String')) returns contents of edit52 as a double


% --- Executes during object creation, after setting all properties.
function edit52_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit53_Callback(hObject, eventdata, handles)
% hObject    handle to edit53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit53 as text
%        str2double(get(hObject,'String')) returns contents of edit53 as a double


% --- Executes during object creation, after setting all properties.
function edit53_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit54_Callback(hObject, eventdata, handles)
% hObject    handle to edit54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit54 as text
%        str2double(get(hObject,'String')) returns contents of edit54 as a double


% --- Executes during object creation, after setting all properties.
function edit54_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
