function varargout = ds(varargin)
% DS MATLAB code for ds.fig
%      DS, by itself, creates a new DS or raises the existing
%      singleton*.
%
%      H = DS returns the handle to a new DS or the handle to
%      the existing singleton*.
%
%      DS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DS.M with the given input arguments.
%
%      DS('Property','Value',...) creates a new DS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ds_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ds_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ds

% Last Modified by GUIDE v2.5 25-Apr-2019 19:35:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ds_OpeningFcn, ...
                   'gui_OutputFcn',  @ds_OutputFcn, ...
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


% --- Executes just before ds is made visible.
function ds_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ds (see VARARGIN)

% Choose default command line output for ds
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ds wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ds_OutputFcn(hObject, eventdata, handles) 
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

azi_num=2000;
fr=1000;

lamda0=0.05;
sigmav=str2double(get(handles.edit1,'string'));
% sigmav=1.0;
sigmaf=2*sigmav/lamda0;

rand('state',sum(100*clock));
d1=rand(1,azi_num);
rand('state',7*sum(100*clock)+3);
d2=rand(1,azi_num);
xi=1*(sqrt(-2*log(d1)).*cos(2*pi*d2));
xq=2*sqrt(-2*log(d1)).*sin(2*pi*d2);

coe_num=12;
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
%Gaussion clutter generation
xxi=conv(b,xi);
xxi=xxi(coe_num*2+1:azi_num+coe_num*2);
xisigmac=std(xxi);
ximuc=mean(xxi);
yyi=(xxi-ximuc)/xisigmac;
muc=str2double(get(handles.edit3,'string'));
sigmac=str2double(get(handles.edit2,'string'));
yyi=sigmac*yyi+log(muc);
xdata=exp(yyi);

plot(handles.axes2,xdata);
xlabel(handles.axes2,'时间');
ylabel(handles.axes2,'幅度');
title(handles.axes2,'对数正态分布杂波时域波形');
hold(handles.axes2,'off');

num=100;
maxdat=max(abs(xdata));
mindat=min(abs(xdata));
NN=hist(abs(xdata),num);
xpdf1=num*NN/((sum(NN)*(maxdat-mindat)));
xaxis1=mindat:(maxdat-mindat)/num:maxdat-(maxdat-mindat)/num;
th_val=lognpdf(xaxis1,log(muc),sigmac);

plot(handles.axes3,xaxis1,xpdf1);
hold(handles.axes3,'on'),plot(handles.axes3,xaxis1,th_val,':r');
title(handles.axes3,'杂波幅度分布');xlabel(handles.axes3,'幅度');ylabel(handles.axes3,'概率密度');
hold(handles.axes3,'off');

signal=xdata;
signal=signal-mean(signal);

M=128;
psd_dat=pburg(real(signal),16,M,fr);
psd_dat=psd_dat/(max(psd_dat));
freqx=0:0.5*M;
freqx=freqx*fr/M;
plot(handles.axes4,freqx,psd_dat);
title(handles.axes4,'杂波频谱');
xlabel(handles.axes4,'频率(Hz)');ylabel(handles.axes4,'功率谱密度');

powerf=exp(-freqx.^2/(2*sigmaf.^2));
hold(handles.axes4,'on');plot(handles.axes4,freqx,powerf,':r');
hold(handles.axes4,'off');



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
