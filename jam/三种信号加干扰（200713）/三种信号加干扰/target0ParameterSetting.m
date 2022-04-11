function varargout = target0ParameterSetting(varargin)
% TARGET0PARAMETERSETTING MATLAB code for target0ParameterSetting.fig
%      TARGET0PARAMETERSETTING, by itself, creates a new TARGET0PARAMETERSETTING or raises the existing
%      singleton*.
%
%      H = TARGET0PARAMETERSETTING returns the handle to a new TARGET0PARAMETERSETTING or the handle to
%      the existing singleton*.
%
%      TARGET0PARAMETERSETTING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TARGET0PARAMETERSETTING.M with the given input arguments.
%
%      TARGET0PARAMETERSETTING('Property','Value',...) creates a new TARGET0PARAMETERSETTING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before target0ParameterSetting_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to target0ParameterSetting_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help target0ParameterSetting

% Last Modified by GUIDE v2.5 27-Oct-2016 11:37:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @target0ParameterSetting_OpeningFcn, ...
                   'gui_OutputFcn',  @target0ParameterSetting_OutputFcn, ...
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


% --- Executes just before target0ParameterSetting is made visible.
function target0ParameterSetting_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to target0ParameterSetting (see VARARGIN)

% Choose default command line output for target0ParameterSetting
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes target0ParameterSetting wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = target0ParameterSetting_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function R_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R as text
%        str2double(get(hObject,'String')) returns contents of R as a double


% --- Executes during object creation, after setting all properties.
function R_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function v_Callback(hObject, eventdata, handles)
% hObject    handle to v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of v as text
%        str2double(get(hObject,'String')) returns contents of v as a double


% --- Executes during object creation, after setting all properties.
function v_CreateFcn(hObject, eventdata, handles)
% hObject    handle to v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigma_Callback(hObject, eventdata, handles)
% hObject    handle to sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigma as text
%        str2double(get(hObject,'String')) returns contents of sigma as a double


% --- Executes during object creation, after setting all properties.
function sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigma0_Callback(hObject, eventdata, handles)
% hObject    handle to sigma0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigma0 as text
%        str2double(get(hObject,'String')) returns contents of sigma0 as a double


% --- Executes during object creation, after setting all properties.
function sigma0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigma0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�״��ź���ʽ��־λ

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��ȡ����Ŀ�����
global RadarS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Ŀ�����

global code;
c=3e8;
% load    �״����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LFM�źŲ���
load data/data_LFMParameter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BK�źŲ���
load data/data_BFParameter;
%global flag;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%JD�źŲ���
load data/data_JDParameter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���߲���
load data/antannaParameter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���߲���
load data/data_target0Parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����Ŀ��ز�
global sigma;
global rcsk;
global sigma0;
sigma0 = 2;
sigma = rcs(rcsk,sigma0);
if  RadarS==1
        %%%%%%%%%%%%LFM�źŻز�����    
        lamta=c/fz_LFM;%����
        ts=1/fs_LFM;
        k=B1_LFM/tau_LFM;                                 %���Ե�Ƶ�źŵ���ϵ��
        tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;  
        N=length(tm);
        f_doppler=2*v/lamta;%��Ŀ�������Ƶ��
        Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %Ŀ��ز��źŹ���
        A=sqrt(Prs);%�ز��źŷ���
        [s_echo_2,echo]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k);
        
        [~,RTAz,RTEl] = xyz2radar(radar.x,radar.y,radar.z,...
                    target.x,target.y,target.z);                        %����Ŀ��Ƕ�
%         [FTR,FTAz,FTEl] = xyz2radar(radar.x1,radar.y1,radar.z1,...
%                     jammer.x3,jammer.y3,jammer.z3);                        %������ŽǶ�
%         Gt = 10^(radar.Gt/10);                                                           %�����״﷢������
%         Gj = jammer.Gain;                                                          %������Ż���������
        Gr = gain(radar,RTAz,RTEl);                                                %Ŀ���Ӧ�Ľ�������
%         Grj = gain(radar,FTAz,FTEl)/Gj;/Gt                                               %���Ż���Ӧ�Ľ�������
        echo1.sum=s_echo_2*Gr.Gainsum;
        echo1.Az=s_echo_2*Gr.GainAz;
        echo1.Ez=s_echo_2*Gr.GainEl;
   
%      figure,plot(0:ts:(N-1)*ts,real(echo1.sum(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('LFM�ز��ź�');
%      figure,plot((0:fs_LFM/length(abs(fft(abs(fftshift(echo(1,:)))))):fs_LFM-fs_LFM/length(abs(fft(abs(fftshift(echo(1,:))))))),abs((fft(echo(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('LFM�ز��źŵ�Ƶ��');     
%        
%         
        figure,plot(0:ts:(N-1)*ts,real(s_echo_2(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('LFM�ز��ź�');
        figure,plot((0:fs_LFM/length(abs(fft(abs(fftshift(echo(1,:)))))):fs_LFM-fs_LFM/length(abs(fft(abs(fftshift(echo(1,:))))))),abs((fft(echo(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('LFM�ز��źŵ�Ƶ��');     
    
else if RadarS==2
        %%%%%%%%%%%��λ�����źŻز�����
        tm_B=0:1/fs_B:tr_B-1/fs_B;%һ�������ظ����ڲ�������
        N=length(tm_B);%һ�������ظ����ڲ�����������
        number1=length(code);
        ts=1/fs_B;
        lamta=c/fz_B;
        Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %Ŀ��ز��źŹ���
        A=sqrt(Prs);%�ز��źŷ���
        f_doppler=2*v/lamta;%��Ŀ�������Ƶ��
        [~,y1]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
        [s_echo_2,echo]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
        [~,RTAz,RTEl] = xyz2radar(radar.x,radar.y,radar.z,...
                    target.x,target.y,target.z);                        %����Ŀ��Ƕ�
%         [FTR,FTAz,FTEl] = xyz2radar(radar.x1,radar.y1,radar.z1,...
%                     jammer.x3,jammer.y3,jammer.z3);                        %������ŽǶ�
%         Gt = 10^(radar.Gt/10);                                                           %�����״﷢������
%         Gj = jammer.Gain;                                                          %������Ż���������
        Gr = gain(radar,RTAz,RTEl);                                                %Ŀ���Ӧ�Ľ�������
%         Grj = gain(radar,FTAz,FTEl)/Gj;/Gt                                               %���Ż���Ӧ�Ľ�������
        echo1.sum=s_echo_2*Gr.Gainsum;
        echo1.Az=s_echo_2*Gr.GainAz;
        echo1.Ez=s_echo_2*Gr.GainEl;
        figure,plot(0:ts:(N-1)*ts,real(echo1.sum(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('��λ����ز��ź�');
        figure,plot((0:fs_B/length(abs(fft(abs(fftshift(echo(1,:)))))):fs_B-fs_B/length(abs(fft(abs(fftshift(echo(1,:))))))),abs((fft(echo(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('��λ����ز��źŵ�Ƶ��');
     
%         figure,plot(0:ts:(N-1)*ts,real(s_echo_2(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('��λ����ز��ź�');
%         figure,plot((0:fs_B/length(abs(fft(abs(fftshift(echo(1,:)))))):fs_B-fs_B/length(abs(fft(abs(fftshift(echo(1,:))))))),abs((fft(echo(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('��λ����ز��źŵ�Ƶ��');
%      
       
    else
        %%%%%%%%%%%�������źŻز�����
        ts=1/fs_JD;
        lamta=c/fz_JD;
        f_doppler=2*v/lamta;
        Prs=(Pt_JD*(10^((Gt_JD*0.1)))*(10^((Gr_JD*0.1)))*(lamta^2)*sigma)/((4*pi)^3*R^4*10^(L_JD*0.1)); %Ŀ��ز��źŹ���
        A=sqrt(Prs);%�ز��źŷ���
        tm=0:1/fs_JD:tr_JD-1/fs_JD;%һ�������ظ����ڲ�������
        N=length(tm);%һ�������ظ����ڲ�����������
%         [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
        [s_echo_2,echo]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
         [~,RTAz,RTEl] = xyz2radar(radar.x,radar.y,radar.z,...
                    target.x,target.y,target.z);                        %����Ŀ��Ƕ�
%         [FTR,FTAz,FTEl] = xyz2radar(radar.x1,radar.y1,radar.z1,...
%                     jammer.x3,jammer.y3,jammer.z3);                        %������ŽǶ�
%         Gt = 10^(radar.Gt/10);                                                           %�����״﷢������
%         Gj = jammer.Gain;                                                          %������Ż���������
        Gr = gain(radar,RTAz,RTEl);                                                %Ŀ���Ӧ�Ľ�������
%         Grj = gain(radar,FTAz,FTEl)/Gj;/Gt                                               %���Ż���Ӧ�Ľ�������
        echo1.sum=s_echo_2*Gr.Gainsum;
        echo1.Az=s_echo_2*Gr.GainAz;
        echo1.Ez=s_echo_2*Gr.GainEl;
%         figure,plot(0:ts:(N-1)*ts,real(echo1.sum(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('������ز��ź�');
%         figure,plot((0:fs_JD/length(abs(fft(abs(fftshift(echo(1,:)))))):fs_JD-fs_JD/length(abs(fft(abs(fftshift(echo(1,:))))))),abs((fft(echo(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('������ز��źŵ�Ƶ��');
%       

        figure,plot(0:ts:(N-1)*ts,real(s_echo_2(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('������ز��ź�');
        figure,plot((0:fs_JD/length(abs(fft(abs(fftshift(echo(1,:)))))):fs_JD-fs_JD/length(abs(fft(abs(fftshift(echo(1,:))))))),abs((fft(echo(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('������ز��źŵ�Ƶ��');
      
    end
    

end

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% load data/data_target0Parameter

 R=str2num(get(handles.R,'string'));
 v=str2num(get(handles.v,'string'));
sigma0=str2num(get(handles.sigma0,'string'));
target.x=str2num(get(handles.x_target,'string'));
target.y=str2num(get(handles.y_target,'string'));
target.z=str2num(get(handles.z_target,'string'));
save data/data_target0Parameter R v sigma0 target;
h1=msgbox('����Ŀ������ɹ�');
pause(1);
close(h1);


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.R,'string','');
set(handles.v,'string','');
set(handles.sigma0,'string','');
set(handles.x_target,'string','');
set(handles.y_target,'string','');
set(handles.z_target,'string','');

% set(handles.L,'string','');

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load data/data_target0Parameter
set(handles.R,'string',R);
set(handles.v,'string',v);
set(handles.sigma0,'string',sigma0);
set(handles.x_target,'string',target.x);
set(handles.y_target,'string',target.y);
set(handles.z_target,'string',target.z);


function L_Callback(hObject, eventdata, handles)
% hObject    handle to L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of L as text
%        str2double(get(hObject,'String')) returns contents of L as a double


% --- Executes during object creation, after setting all properties.
function L_CreateFcn(hObject, eventdata, handles)
% hObject    handle to L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x_target_Callback(hObject, eventdata, handles)
% hObject    handle to x_target (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x=str2num(get(handles.x_target,'string'));
y=str2num(get(handles.y_target,'string'));
z=str2num(get(handles.z_target,'string'));
r=sqrt(x^2+y^2+z^2);
set(handles.R,'string',r);
% Hints: get(hObject,'String') returns contents of x_target as text
%        str2double(get(hObject,'String')) returns contents of x_target as a double


% --- Executes during object creation, after setting all properties.
function x_target_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_target (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y_target_Callback(hObject, eventdata, handles)
% hObject    handle to y_target (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x=str2num(get(handles.x_target,'string'));
y=str2num(get(handles.y_target,'string'));
z=str2num(get(handles.z_target,'string'));
r=sqrt(x^2+y^2+z^2);
set(handles.R,'string',r);
% Hints: get(hObject,'String') returns contents of y_target as text
%        str2double(get(hObject,'String')) returns contents of y_target as a double


% --- Executes during object creation, after setting all properties.
function y_target_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_target (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function z_target_Callback(hObject, eventdata, handles)
% hObject    handle to z_target (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x=str2num(get(handles.x_target,'string'));
y=str2num(get(handles.y_target,'string'));
z=str2num(get(handles.z_target,'string'));
r=sqrt(x^2+y^2+z^2);
set(handles.R,'string',r);
% Hints: get(hObject,'String') returns contents of z_target as text
%        str2double(get(hObject,'String')) returns contents of z_target as a double


% --- Executes during object creation, after setting all properties.
function z_target_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z_target (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
