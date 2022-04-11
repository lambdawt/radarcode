function varargout = rl(varargin)
% RL MATLAB code for rl.fig
%      RL, by itself, creates a new RL or raises the existing
%      singleton*.
%
%      H = RL returns the handle to a new RL or the handle to
%      the existing singleton*.
%
%      RL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RL.M with the given input arguments.
%
%      RL('Property','Value',...) creates a new RL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before rl_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to rl_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help rl

% Last Modified by GUIDE v2.5 24-Apr-2019 13:09:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @rl_OpeningFcn, ...
                   'gui_OutputFcn',  @rl_OutputFcn, ...
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


% --- Executes just before rl is made visible.
function rl_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to rl (see VARARGIN)

% Choose default command line output for rl
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes rl wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = rl_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all;
run('Radarclutter');


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

azi_num=2000;   %取2000个点
fr=1000;        %雷达重复频率

lamda0=0.05;   %杂波波长
sigmav=str2double(get(handles.edit1,'string')); %杂波方差   
sigmaf=2*sigmav/lamda0;  %杂波功率谱的标准偏差


rand('state',sum(100*clock)); %产生服从U(0,1)分布的随机序列
d1=rand(1,azi_num);            
rand('state',7*sum(100*clock)+3);
d2=rand(1,azi_num);
xi=2*sqrt(-2*log(d1)).*cos(2*pi*d2);  %正交且独立的高斯序列N(0，1)
xq=2*sqrt(-2*log(d1)).*sin(2*pi*d2);
%形成滤波器频率响应
coe_num=12;           %求滤波器系数，用傅里叶级数展开法
for n=0:coe_num
    coeff(n+1)=2*sigmaf*sqrt(pi)*exp(-4*sigmaf^2*pi^2*n^2/fr^2)/fr;  
end
for n=1:2*coe_num+1
    if n<=coe_num+1
        b(n)=1/2*coeff(coe_num+2-n);
    else
        b(n)=1/2*coeff(n-coe_num);
    end
end
%生成高斯谱杂波
xxi=conv(b,xi);   
xxq=conv(b,xq);   
xxi=xxi(coe_num*2+1:azi_num+coe_num*2);%目的是去掉暂态响应
xxq=xxq(coe_num*2+1:azi_num+coe_num*2);
xisigmac=std(xxi);     
ximuc=mean(xxi);       
yyi=(xxi-ximuc)/xisigmac;    
xqsigmac=std(xxq);     
xqmuc=mean(xxq);       
yyq=(xxq-xqmuc)/xqsigmac;    %归一化
sigmac=1.2 ;           %杂波的速度标准差
yyi=sigmac*yyi;        %使瑞利分布杂波具有指定的标准差
yyq=sigmac*yyq;        %使瑞利分布虚部杂波
ydata=yyi+j*yyq;       %瑞利分布杂波形成

plot(handles.axes1,real(ydata)); %瑞利分布杂波实部
hold(handles.axes1,'off');
title(handles.axes1,'瑞利杂波时域波形，实部');  
plot(handles.axes4,imag(ydata)); %瑞利分布杂波虚部
title(handles.axes4,'瑞利杂波时域波形，虚部');
hold(handles.axes4,'off');

num=100;                   %求概率密度函数的参数
maxdat=max(abs(ydata));
mindat=min(abs(ydata));
NN=hist(abs(ydata),num);   
xpdf1=num*NN/((sum(NN))*(maxdat-mindat));        %用直方图估计的概率密度函数
xaxisl=mindat:(maxdat-mindat)/num:maxdat-(maxdat-mindat)/num;  
th_val=(xaxisl./sigmac.^2).*exp(-xaxisl.^2./(2*sigmac.^2));   %概率密度函数理论值

plot(handles.axes2,xaxisl,xpdf1);               %做出仿真结果的概率密度函数曲线
hold(handles.axes2,'on');plot(handles.axes2,xaxisl,th_val,'r:'); %做出理论概率密度函数曲线
title(handles.axes2,'杂波幅度分布');
xlabel(handles.axes2,'幅度');ylabel(handles.axes2,'概率密度');
hold(handles.axes2,'off');

signal=ydata;
signal=signal-mean(signal);      %求功率谱密度，先去掉直流分量
M=256;                     %用burg法估计功率谱密度
psd_dat=pburg(real(signal),32,M,fr);
psd_dat=psd_dat/(max(psd_dat));   %归一化处理
freqx=0:0.5*M;
freqx=freqx*fr/M;
plot(handles.axes3,freqx,psd_dat);title(handles.axes3,'杂波频谱');
xlabel(handles.axes3,'频率/HZ');ylabel(handles.axes3,'功率谱密度');

%做出理想高斯谱曲线
powerf=exp(-freqx.^2/(2*sigmaf.^2));
hold(handles.axes3,'on');plot(handles.axes3,freqx,powerf,'r:');
hold(handles.axes3,'off');
