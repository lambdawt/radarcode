function varargout = DeceptionJammingParameter(varargin)
% DECEPTIONJAMMINGPARAMETER MATLAB code for DeceptionJammingParameter.fig
%      DECEPTIONJAMMINGPARAMETER, by itself, creates a new DECEPTIONJAMMINGPARAMETER or raises the existing
%      singleton*.
%
%      H = DECEPTIONJAMMINGPARAMETER returns the handle to a new DECEPTIONJAMMINGPARAMETER or the handle to
%      the existing singleton*.
%
%      DECEPTIONJAMMINGPARAMETER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DECEPTIONJAMMINGPARAMETER.M with the given input arguments.
%
%      DECEPTIONJAMMINGPARAMETER('Property','Value',...) creates a new DECEPTIONJAMMINGPARAMETER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DeceptionJammingParameter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DeceptionJammingParameter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DeceptionJammingParameter

% Last Modified by GUIDE v2.5 02-Nov-2016 22:19:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DeceptionJammingParameter_OpeningFcn, ...
                   'gui_OutputFcn',  @DeceptionJammingParameter_OutputFcn, ...
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


% --- Executes just before DeceptionJammingParameter is made visible.
function DeceptionJammingParameter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DeceptionJammingParameter (see VARARGIN)

% Choose default command line output for DeceptionJammingParameter
handles.output = hObject;

global sigma0;
sigma0 = 2;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DeceptionJammingParameter wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DeceptionJammingParameter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global RadarS;
global sigma;
global rcsk;
rcsk = 2;
global sigma0;
sigma = rcs(rcsk,sigma0);
global string1;
global string2;
global temp1;
global temp2;
c=3e8;
Gj=str2num(get(handles.Gj,'string'));
Gjr=str2num(get(handles.Gjr,'string'));
Pj=str2num(get(handles.Pj,'string'));
tf=str2num(get(handles.tf,'string'));
vf=str2num(get(handles.vf,'string'));
L=3;%�ۺ����
load data/data_target0Parameter
if RadarS == 1
string1=get(handles.n1,'string');
string2=get(handles.n2,'string');
temp1=get(handles.n1,'Value');
temp2=get(handles.n2,'Value');
load data/data_LFMParameter;
switch string1{temp1}
    case '�ٶ�'
        switch string2{temp2}
            case {'��' ,'��','�ܼ�' }
               R1=str2num(get(handles.R1,'string'));
               v1=str2num(get(handles.v1,'string'));
               congmubiao=str2num(get(handles.congmubiao,'string'));
               fr=1/tr_LFM;%�����ظ�Ƶ
                lamta=c/fz_LFM;%����
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%һ�������ظ����ڲ�������
                N=length(tm);%һ�������ظ����ڲ�����������
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%����ǿ��
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %#ok<*NODEF> %Ŀ��ز��źŹ���
                A=sqrt(Prs);%�ز��źŷ���
                f_doppler=2*v/lamta;%��Ŀ�������Ƶ��%���Ե�Ƶ�źŵ���ϵ��
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
                Aj=sqrt(Prj);
                %�ٶ�
                f_doppler1=2*v1/lamta;
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [~,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y,tr_LFM,ts,N);
                [s_echo_2,echo]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                [s_ft,echo3]=LFMDeceptionJamming(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao,temp1);
                s_echo_1=s_echo_2+s_noise+s_ft;%%%%�����ź�
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('�״�����ź�');
%                 figure,plot(linspace(1.2e-5,1.3e-5,N),real(s_echo_1(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('�״�����ź�');
                figure,plot((0:fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_LFM-fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�״�����źŵ�Ƶ��');
                save data/data_DeceptionJammingParameter R1 v1 congmubiao Gj Gjr Pj   L
            case '����'
%                  load data/data_DeceptionJammingParameter
                tf=10;
                c=3e8;
                lamta=c/fz_LFM;
                Rj=2e3;
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L_LFM/10)));
                Aj=sqrt(Prj);
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%һ�������ظ����ڲ�������
                N=length(tm);%һ�������ظ����ڲ�����������
                f_doppler=2*v/lamta;
                f_doppler1=f_doppler;
                k=B1_LFM/tau_LFM;%���Ե�Ƶ�źŵ���ϵ��
                ts=1/fs_LFM;%�������
                %���ɷ����ź�
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L/10))); %Ŀ��ز��źŹ���
                A=sqrt(Prs);
                Te=1143;%�¶�
                F=4;%����ϵ��
                An=10*log10((1.382e-23)*Te*B_LFM*10^(F/10));%����ǿ��
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k); 
                [~,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y,tr_LFM,ts,N);
                %���ɻز��ź�
                [s_echo_2,echo]=LFMtuoyinhuiboxinhao(y,R,tr_LFM,ts,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                %figure,plot(0:ts:(N-1)*ts,real(s_echo_2(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('�ز��ź�');
                %figure,plot((0:fs/length(abs(fft(abs(fftshift(echo(1,:)))))):fs-fs/length(abs(fft(abs(fftshift(echo(1,:))))))),abs((fft(echo(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�ز��źŵ�Ƶ��');

                %���ɸ����ź�
                [s_ft,echo3]=LFMtuoyinganrao(R,frame_LFM,tf,Aj,tm,tau_LFM,f0_LFM,B1_LFM,N,f_doppler1,f_doppler,vf,lamta,k,ts,c,fs_LFM,temp1);

                %��������
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                %Ŀ��ز��źš���Ŀ���źš�����������һ��������ջ�
                s_echo_1=s_echo_2+s_noise+s_ft;
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('�״�����ź�');
                figure,plot((0:fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_LFM-fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�״�����źŵ�Ƶ��');
                save data/data_DeceptionJammingParameter  Gj Gjr Pj tf vf L
        end
    case '����'
         switch string2{temp2}
            case {'��' ,'��','�ܼ�' }
                R1=str2num(get(handles.R1,'string'));
                v1=str2num(get(handles.v1,'string'));
                congmubiao=str2num(get(handles.congmubiao,'string'));
                fr=1/tr_LFM;%�����ظ�Ƶ
                lamta=c/fz_LFM;%����
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%һ�������ظ����ڲ�������
                N=length(tm);%һ�������ظ����ڲ�����������
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%����ǿ��
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %Ŀ��ز��źŹ���
                A=sqrt(Prs);%�ز��źŷ���
                f_doppler=2*v/lamta;%��Ŀ�������Ƶ��%���Ե�Ƶ�źŵ���ϵ��
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
                Aj=sqrt(Prj);
                %�ٶ�
                f_doppler1=2*v1/lamta;
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [~,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y,tr_LFM,ts,N);
                [s_echo_2,echo]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                [s_ft,echo3]=LFMDeceptionJamming(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao,temp1);
                s_echo_1=s_echo_2+s_noise+s_ft;%%%%�����ź�
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('�״�����ź�');
                figure,plot((0:fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_LFM-fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�״�����źŵ�Ƶ��');
                save data/data_DeceptionJammingParameter R1 v1 congmubiao Gj Gjr Pj tf vf L
           
            case '����'
                load data/data_target0Parameter
                c=3e8;
                lamta=c/fz_LFM;
                Rj=2e3;
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L_LFM/10)));
                Aj=sqrt(Prj);
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%һ�������ظ����ڲ�������
                N=length(tm);%һ�������ظ����ڲ�����������
                f_doppler=2*v/lamta;
                f_doppler1=f_doppler;
                k=B1_LFM/tau_LFM;%���Ե�Ƶ�źŵ���ϵ��
                ts=1/fs_LFM;%�������
                %���ɷ����ź�
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L/10))); %Ŀ��ز��źŹ���
                A=sqrt(Prs);
                Te=1143;%�¶�
                F=4;%����ϵ��
                An=10*log10((1.382e-23)*Te*B_LFM*10^(F/10));%����ǿ��
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k); 
                [~,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y,tr_LFM,ts,N);
                %���ɻز��ź�
                [s_echo_2,echo]=LFMtuoyinhuiboxinhao(y,R,tr_LFM,ts,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                %figure,plot(0:ts:(N-1)*ts,real(s_echo_2(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('�ز��ź�');
                %figure,plot((0:fs/length(abs(fft(abs(fftshift(echo(1,:)))))):fs-fs/length(abs(fft(abs(fftshift(echo(1,:))))))),abs((fft(echo(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�ز��źŵ�Ƶ��');

                %���ɸ����ź�
                [s_ft,echo3]=LFMtuoyinganrao(R,frame_LFM,tf,Aj,tm,tau_LFM,f0_LFM,B1_LFM,N,f_doppler1,f_doppler,vf,lamta,k,ts,c,fs_LFM,temp1);

                %��������
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                %Ŀ��ز��źš���Ŀ���źš�����������һ��������ջ�
                s_echo_1=s_echo_2+s_noise+s_ft;
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('�״�����ź�');
                figure,plot((0:fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_LFM-fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�״�����źŵ�Ƶ��');
                save data/data_DeceptionJammingParameter  Gj Gjr Pj tf vf L
         end
    case '����'
         switch string2{temp2}
            case {'��' ,'��','�ܼ�' }
                 R1=str2num(get(handles.R1,'string'));
                v1=str2num(get(handles.v1,'string'));
                congmubiao=str2num(get(handles.congmubiao,'string'));
                fr=1/tr_LFM;%�����ظ�Ƶ
                lamta=c/fz_LFM;%����
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%һ�������ظ����ڲ�������
                N=length(tm);%һ�������ظ����ڲ�����������
                An=10*log10((1.382e-23)*Te_LFM*B_LFM*10^(F_LFM/10));%����ǿ��
                ts=1/fs_LFM;
                k=B1_LFM/tau_LFM;   
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %Ŀ��ز��źŹ���
                A=sqrt(Prs);%�ز��źŷ���
                f_doppler=2*v/lamta;%��Ŀ�������Ƶ��%���Ե�Ƶ�źŵ���ϵ��
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
                Aj=sqrt(Prj);
                %�ٶ�
                f_doppler1=2*v1/lamta;
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k);
                [~,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y,tr_LFM,ts,N);
                [s_echo_2,echo]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                [s_ft,echo3]=LFMDeceptionJamming(D,y,R1,tr_LFM,ts,c,Aj,N,frame_LFM,fs_LFM,f_doppler1,tm,f0_LFM,B1_LFM,tau_LFM,k,congmubiao,temp1);
                s_echo_1=s_echo_2+s_noise+s_ft;%%%%�����ź�
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('�״�����ź�');
                figure,plot((0:fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_LFM-fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�״�����źŵ�Ƶ��');
                save data/data_DeceptionJammingParameter R1 v1 congmubiao Gj Gjr Pj tf vf L
            
            case '����'
                load data/data_target0Parameter
                c=3e8;
                lamta=c/fz_LFM;
                Rj=2e3;
                Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L_LFM/10)));
                Aj=sqrt(Prj);
                tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;%һ�������ظ����ڲ�������
                N=length(tm);%һ�������ظ����ڲ�����������
                f_doppler=2*v/lamta;
                f_doppler1=f_doppler;
                k=B1_LFM/tau_LFM;%���Ե�Ƶ�źŵ���ϵ��
                ts=1/fs_LFM;%�������
                %���ɷ����ź�
                Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L/10))); %Ŀ��ز��źŹ���
                A=sqrt(Prs);
                Te=1143;%�¶�
                F=4;%����ϵ��
                An=10*log10((1.382e-23)*Te*B_LFM*10^(F/10));%����ǿ��
                [y,D]=shengchengLFMxinhao(B1_LFM,Pt_LFM,tau_LFM,f0_LFM,tm,k); 
                [~,match_filter_fft]=maiyaxishu(f0_LFM,fs_LFM,y,tr_LFM,ts,N);
                %���ɻز��ź�
                [s_echo_2,echo]=LFMtuoyinhuiboxinhao(y,R,tr_LFM,ts,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
                %figure,plot(0:ts:(N-1)*ts,real(s_echo_2(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('�ز��ź�');
                %figure,plot((0:fs/length(abs(fft(abs(fftshift(echo(1,:)))))):fs-fs/length(abs(fft(abs(fftshift(echo(1,:))))))),abs((fft(echo(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�ز��źŵ�Ƶ��');

                %���ɸ����ź�
                [s_ft,echo3]=LFMtuoyinganrao(R,frame_LFM,tf,Aj,tm,tau_LFM,f0_LFM,B1_LFM,N,f_doppler1,f_doppler,vf,lamta,k,ts,c,fs_LFM,temp1);

                %��������
                [s_noise]=zaosheng(frame_LFM,N,An,B_LFM,fs_LFM);
                %Ŀ��ز��źš���Ŀ���źš�����������һ��������ջ�
                s_echo_1=s_echo_2+s_noise+s_ft;
                figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('�״�����ź�');
                figure,plot((0:fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_LFM-fs_LFM/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�״�����źŵ�Ƶ��');
                save data/data_DeceptionJammingParameter  Gj Gjr Pj tf vf L
         end
         
end
elseif RadarS==2
                        
      string1=get(handles.n1,'string');
      string2=get(handles.n2,'string');
      temp1=get(handles.n1,'Value');
      temp2=get(handles.n2,'Value');
      load data/data_BFParameter;
      global code; %#ok<*TLEV>
     switch string1{temp1}
    case '�ٶ�'
        switch string2{temp2}
        case {'��' ,'��','�ܼ�' }
        R1=str2num(get(handles.R1,'string'));
        v1=str2num(get(handles.v1,'string'));
        congmubiao=str2num(get(handles.congmubiao,'string'));
        tm_B=0:1/fs_B:tr_B-1/fs_B;%һ�������ظ����ڲ�������
        N=length(tm_B);%һ�������ظ����ڲ�����������
        number1=length(code);
        ts=1/fs_B;
        lamta=c/fz_B;
        fr=1/tr_B;
        An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%����ǿ��
        Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %Ŀ��ز��źŹ���
        A=sqrt(Prs);%�ز��źŷ���
        f_doppler=2*v/lamta;%��Ŀ�������Ƶ��
        Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
        Aj=sqrt(Prj);
        %�ٶ�
        f_doppler1=2*v1/lamta; 
        [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
        [s_echo_2,echo]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
        [~,match_filter_fft]=maiyaxishu(f0_B,fs_B,y,tr_B,ts,N);
        [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
        [s_ft,echo3]=BKDeceptionJamming(D,y,R1,tr_B,ts,c,Aj,N,frame_B,fs_B,f_doppler1,tm_B,f0_B,tau_B,congmubiao,y1,temp1);
        s_echo_1=s_echo_2+s_noise+s_ft;%%%%�����ź�
        figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('�״�����ź�');
        figure,plot((0:fs_B/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_B-fs_B/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�״�����źŵ�Ƶ��');
        save data/data_DeceptionJammingParameter R1 v1 congmubiao Gj Gjr Pj tf vf L
        case '����'
            c=3e8;
            Rj=2e3;
            number1=length(code);
            ts=1/fs_B;
            lamta=c/fz_B;
            tm=0:1/fs_B:tr_B-1/fs_B;%һ�������ظ����ڲ�������
            N=length(tm);%һ�������ظ����ڲ�����������
            f_doppler=2*v/lamta;%��Ŀ�������Ƶ��
            f_doppler1= f_doppler;%��Ŀ�������Ƶ��
            Te=1143;%�¶�
            An=10*log10((1.382e-23)*Te*B_B*10^(F_B/10));%����ǿ��
            Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %Ŀ��ز��źŹ���
            A=sqrt(Prs);%�ز��źŷ���
            Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L/10)));
            Aj=sqrt(Prj);
            [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
%             figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(��λ����)'), ylabel('y(��λ����)'),title('�����ź�');
%             figure,plot((0:fs_B/length(abs(fft(abs(fftshift(y))))):fs_B-fs_B/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�����źŵ�Ƶ��');
%             %%%%%%%%%%%1.1.2������������ѹ��ϵ��
            [~,match_filter_fft]=maiyaxishu(f0_B,fs_B,y,tr_B,ts,N);
             %%%%%%%%%%%%%1.3���ɻز��ź�%%%%%%%%%%%%
            [s_echo_2,echo]=BKtuoyinhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tm,tau_B); 
            %��������
            [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
            %�����ź�
            [s_ft,echo3]=BKtuoyinganrao(R,frame_B,tf,Aj,tau_B,N,f_doppler1,f_doppler,vf,lamta,ts,c,fs_B,y1,tr_B,temp1);
            %Ŀ��ز��źš���Ŀ���źš�����������һ��������ջ�
            s_echo_1=s_echo_2+s_noise+s_ft;
            figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('�״�����ź�');
            figure,plot((0:fs_B/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_B-fs_B/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�״�����źŵ�Ƶ��');
            save data/data_DeceptionJammingParameter  Gj Gjr Pj tf vf L
        end
  case '����'
        switch string2{temp2}
        case {'��' ,'��','�ܼ�' }
        R1=str2num(get(handles.R1,'string'));
        v1=str2num(get(handles.v1,'string'));
        congmubiao=str2num(get(handles.congmubiao,'string'));
        load data/data_BFParameter;
        tm_B=0:1/fs_B:tr_B-1/fs_B;%һ�������ظ����ڲ�������
        N=length(tm_B);%һ�������ظ����ڲ�����������
        number1=length(code);
        ts=1/fs_B;
        lamta=c/fz_B;
        fr=1/tr_B;
        An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%����ǿ��
        Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %Ŀ��ز��źŹ���
        A=sqrt(Prs);%�ز��źŷ���
        f_doppler=2*v/lamta;%��Ŀ�������Ƶ��
        Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
        Aj=sqrt(Prj);
        %�ٶ�
        f_doppler1=2*v1/lamta; 
        [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
        [s_echo_2,echo]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
        [~,match_filter_fft]=maiyaxishu(f0_B,fs_B,y,tr_B,ts,N);
        [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
        [s_ft,echo3]=BKDeceptionJamming(D,y,R1,tr_B,ts,c,Aj,N,frame_B,fs_B,f_doppler1,tm_B,f0_B,tau_B,congmubiao,y1,temp1);
        s_echo_1=s_echo_2+s_noise+s_ft;%%%%�����ź�
        figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('�״�����ź�');
        figure,plot((0:fs_B/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_B-fs_B/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�״�����źŵ�Ƶ��');
        save data/data_DeceptionJammingParameter R1 v1 congmubiao Gj Gjr Pj tf vf L
        case '����'
           c=3e8;
            Rj=2e3;
            number1=length(code);
            ts=1/fs_B;
            lamta=c/fz_B;
            tm=0:1/fs_B:tr_B-1/fs_B;%һ�������ظ����ڲ�������
            N=length(tm);%һ�������ظ����ڲ�����������
            f_doppler=2*v/lamta;%��Ŀ�������Ƶ��
            f_doppler1= f_doppler;%��Ŀ�������Ƶ��
            Te=1143;%�¶�
            An=10*log10((1.382e-23)*Te*B_B*10^(F_B/10));%����ǿ��
            Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %Ŀ��ز��źŹ���
            A=sqrt(Prs);%�ز��źŷ���
            Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L/10)));
            Aj=sqrt(Prj);
            [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
%             figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(��λ����)'), ylabel('y(��λ����)'),title('�����ź�');
%             figure,plot((0:fs_B/length(abs(fft(abs(fftshift(y))))):fs_B-fs_B/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�����źŵ�Ƶ��');
%             %%%%%%%%%%%1.1.2������������ѹ��ϵ��
            [~,match_filter_fft]=maiyaxishu(f0_B,fs_B,y,tr_B,ts,N);
             %%%%%%%%%%%%%1.3���ɻز��ź�%%%%%%%%%%%%
            [s_echo_2,echo]=BKtuoyinhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tm,tau_B); 
            %��������
            [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
            %�����ź�
            [s_ft,echo3]=BKtuoyinganrao(R,frame_B,tf,Aj,tau_B,N,f_doppler1,f_doppler,vf,lamta,ts,c,fs_B,y1,tr_B,temp1);
            %Ŀ��ز��źš���Ŀ���źš�����������һ��������ջ�
            s_echo_1=s_echo_2+s_noise+s_ft;
            figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('�״�����ź�');
            figure,plot((0:fs_B/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_B-fs_B/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�״�����źŵ�Ƶ��');
            save data/data_DeceptionJammingParameter  Gj Gjr Pj tf vf L
        end
       case '����'
        switch string2{temp2}
        case {'��' ,'��','�ܼ�' }
        R1=str2num(get(handles.R1,'string'));
        v1=str2num(get(handles.v1,'string'));
        congmubiao=str2num(get(handles.congmubiao,'string'));
        load data/data_BFParameter;
        tm_B=0:1/fs_B:tr_B-1/fs_B;%һ�������ظ����ڲ�������
        N=length(tm_B);%һ�������ظ����ڲ�����������
        number1=length(code);
        ts=1/fs_B;
        lamta=c/fz_B;
        fr=1/tr_B;
        An=10*log10((1.382e-23)*Te_B*B_B*10^(F_B/10));%����ǿ��
        Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %Ŀ��ز��źŹ���
        A=sqrt(Prs);%�ز��źŷ���
        f_doppler=2*v/lamta;%��Ŀ�������Ƶ��
        Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
        Aj=sqrt(Prj);
        %�ٶ�
        f_doppler1=2*v1/lamta; 
        [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
        [s_echo_2,echo]=BKhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tau_B); 
        [~,match_filter_fft]=maiyaxishu(f0_B,fs_B,y,tr_B,ts,N);
        [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
        [s_ft,echo3]=BKDeceptionJamming(D,y,R1,tr_B,ts,c,Aj,N,frame_B,fs_B,f_doppler1,tm_B,f0_B,tau_B,congmubiao,y1,temp1);
        s_echo_1=s_echo_2+s_noise+s_ft;%%%%�����ź�
        figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('�״�����ź�');
        figure,plot((0:fs_B/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_B-fs_B/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�״�����źŵ�Ƶ��');
        save data/data_DeceptionJammingParameter R1 v1 congmubiao Gj Gjr Pj tf vf L
      case '����'
           c=3e8;
            Rj=2e3;
            number1=length(code);
            ts=1/fs_B;
            lamta=c/fz_B;
            tm=0:1/fs_B:tr_B-1/fs_B;%һ�������ظ����ڲ�������
            N=length(tm);%һ�������ظ����ڲ�����������
            f_doppler=2*v/lamta;%��Ŀ�������Ƶ��
            f_doppler1= f_doppler;%��Ŀ�������Ƶ��
            Te=1143;%�¶�
            An=10*log10((1.382e-23)*Te*B_B*10^(F_B/10));%����ǿ��
            Prs=((Pt_B*(10^((Gt_B/10)))*(10^((Gr_B/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_B/10))); %Ŀ��ز��źŹ���
            A=sqrt(Prs);%�ز��źŷ���
            Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L/10)));
            Aj=sqrt(Prj);
            [y,y1,D]=shengchengBKxinhao(tau_B,fs_B,f0_B,flag,number1,code,Pt_B,tr_B,ts);
%             figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(��λ����)'), ylabel('y(��λ����)'),title('�����ź�');
%             figure,plot((0:fs_B/length(abs(fft(abs(fftshift(y))))):fs_B-fs_B/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�����źŵ�Ƶ��');
%             %%%%%%%%%%%1.1.2������������ѹ��ϵ��
            [~,match_filter_fft]=maiyaxishu(f0_B,fs_B,y,tr_B,ts,N);
             %%%%%%%%%%%%%1.3���ɻز��ź�%%%%%%%%%%%%
            [s_echo_2,echo]=BKtuoyinhuiboxinhao(y1,R,tr_B,ts,A,N,frame_B,fs_B,f_doppler,tm,tau_B); 
            %��������
            [s_noise]=zaosheng(frame_B,N,An,B_B,fs_B);
            %�����ź�
            [s_ft,echo3]=BKtuoyinganrao(R,frame_B,tf,Aj,tau_B,N,f_doppler1,f_doppler,vf,lamta,ts,c,fs_B,y1,tr_B,temp1);
            %Ŀ��ز��źš���Ŀ���źš�����������һ��������ջ�
            s_echo_1=s_echo_2+s_noise+s_ft;
            figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('�״�����ź�');
            figure,plot((0:fs_B/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_B-fs_B/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�״�����źŵ�Ƶ��');
            save data/data_DeceptionJammingParameter  Gj Gjr Pj tf vf L
        end
     end
elseif RadarS==3
     string1=get(handles.n1,'string');
      string2=get(handles.n2,'string');
      temp1=get(handles.n1,'Value');
      temp2=get(handles.n2,'Value');
      load data/data_JDParameter;
     switch string1{temp1}
     case '�ٶ�'
        switch string2{temp2}
        case {'��' ,'��','�ܼ�' }
        R1=str2num(get(handles.R1,'string'));
        v1=str2num(get(handles.v1,'string'));
        congmubiao=str2num(get(handles.congmubiao,'string'));
        fr=1/tr_JD;
        lamta=c/fz_JD;%����
        tm=0:1/fs_JD:tr_JD-1/fs_JD;%һ�������ظ����ڲ�������
        N=length(tm);%һ�������ظ����ڲ�����������
        An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%����ǿ��
        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        ts=1/fs_JD; 
        Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %Ŀ��ز��źŹ���
        A=sqrt(Prs);%�ز��źŷ���
        f_doppler=2*v/lamta;%��Ŀ�������Ƶ��
        Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
        Aj=sqrt(Prj);
        %�ٶ�
         f_doppler1=2*v1/lamta; 
        [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
        [~,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y,tr_JD,ts,N);
        [s_echo_2,echo]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        [s_ft,echo3]=JDDeceptionJamming(D,R1,tr_JD,c,Aj,N,frame_JD,fs_JD,f_doppler1,tm,f0_JD,tau_JD,congmubiao,temp1);
        s_echo_1=s_echo_2+s_noise+s_ft;%%%%�����ź�
        figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('�״�����ź�');
        figure,plot((0:fs_JD/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_JD-fs_JD/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�״�����źŵ�Ƶ��');
        save data/data_DeceptionJammingParameter R1 v1 congmubiao Gj Gjr Pj tf vf L
        
        case '����'
             %���ɷ����ź�
        Rj=2e3;
        c=3e8;
        lamta=c/fz_JD;%����
        tm=0:1/fs_JD:tr_JD-1/fs_JD;%һ�������ظ����ڲ�������
        N=length(tm);%һ�������ظ����ڲ�����������
        An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%����ǿ��
        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        ts=1/fs_JD; 
        Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %Ŀ��ز��źŹ���
        A=sqrt(Prs);%�ز��źŷ���
        f_doppler=2*v/lamta;%��Ŀ�������Ƶ��
   
        %�ٶ�
        f_doppler=2*v/lamta;
        f_doppler1=f_doppler;
        lamta=c/fz_JD;
        f_doppler=2*v/lamta;
        [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
        ts=1/fs_JD;%�������
        tm=0:1/fs_JD:tr_JD-1/fs_JD;%һ�������ظ����ڲ�������
        N=length(tm);%һ�������ظ����ڲ�����������
      
        Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L/10)));
        Aj=sqrt(Prj);
%         figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(��λ����)'), ylabel('y(��λ����)'),title('�����ź�');
%         figure,plot((0:fs/length(abs(fft(abs(fftshift(y))))):fs-fs/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�����źŵ�Ƶ��');
        %������������ѹ��ϵ��
        [~,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y,tr_JD,ts,N);

        %���ɻز��ź�
        [s_echo_2,echo]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
             [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        %���ɸ����ź�
        [s_ft,echo3]=JDtuoyinganrao(R,frame_JD,tf,Aj,tm,tau_JD,f0_JD,N,f_doppler1,f_doppler,vf,lamta,ts,c,fs_JD,temp1);

        %Ŀ��ز��źš���Ŀ���źš�����������һ��������ջ�
        s_echo_1=s_echo_2+s_noise+s_ft;
        figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('�״�����ź�');
        figure,plot((0:fs_JD/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_JD-fs_JD/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�״�����źŵ�Ƶ��');
        save data/data_DeceptionJammingParameter  Gj Gjr Pj tf vf L
        
        end
    case '����'
        switch string2{temp2}
        case {'��' ,'��','�ܼ�' }
        R1=str2num(get(handles.R1,'string'));
        v1=str2num(get(handles.v1,'string'));
        congmubiao=str2num(get(handles.congmubiao,'string'));
        load data/data_JDParameter;
        fr=1/tr_JD;
        lamta=c/fz_JD;%����
        tm=0:1/fs_JD:tr_JD-1/fs_JD;%һ�������ظ����ڲ�������
        N=length(tm);%һ�������ظ����ڲ�����������
        An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%����ǿ��
        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        ts=1/fs_JD; 
        Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %Ŀ��ز��źŹ���
        A=sqrt(Prs);%�ز��źŷ���
        f_doppler=2*v/lamta;%��Ŀ�������Ƶ��
        Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
        Aj=sqrt(Prj);
        %�ٶ�
         f_doppler1=2*v1/lamta; 
        [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
        [~,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y,tr_JD,ts,N);
        [s_echo_2,echo]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        [s_ft,echo3]=JDDeceptionJamming(D,R1,tr_JD,c,Aj,N,frame_JD,fs_JD,f_doppler1,tm,f0_JD,tau_JD,congmubiao,temp1);
        s_echo_1=s_echo_2+s_noise+s_ft;%%%%�����ź�
        figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('�״�����ź�');
        figure,plot((0:fs_JD/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_JD-fs_JD/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�״�����źŵ�Ƶ��');
        save data/data_DeceptionJammingParameter R1 v1 congmubiao Gj Gjr Pj tf vf L
        
        case '����'
             %���ɷ����ź�
        Rj=2e3;
        c=3e8;
        lamta=c/fz_JD;%����
        tm=0:1/fs_JD:tr_JD-1/fs_JD;%һ�������ظ����ڲ�������
        N=length(tm);%һ�������ظ����ڲ�����������
        An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%����ǿ��
        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        ts=1/fs_JD; 
        Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %Ŀ��ز��źŹ���
        A=sqrt(Prs);%�ز��źŷ���
        f_doppler=2*v/lamta;%��Ŀ�������Ƶ��
   
        %�ٶ�
        f_doppler=2*v/lamta;
        f_doppler1=f_doppler;
        lamta=c/fz_JD;
        f_doppler=2*v/lamta;
        [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
        ts=1/fs_JD;%�������
        tm=0:1/fs_JD:tr_JD-1/fs_JD;%һ�������ظ����ڲ�������
        N=length(tm);%һ�������ظ����ڲ�����������
      
        Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L/10)));
        Aj=sqrt(Prj);
%         figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(��λ����)'), ylabel('y(��λ����)'),title('�����ź�');
%         figure,plot((0:fs/length(abs(fft(abs(fftshift(y))))):fs-fs/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�����źŵ�Ƶ��');
        %������������ѹ��ϵ��
        [~,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y,tr_JD,ts,N);

        %���ɻز��ź�
        [s_echo_2,echo]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
             [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        %���ɸ����ź�
        [s_ft,echo3]=JDtuoyinganrao(R,frame_JD,tf,Aj,tm,tau_JD,f0_JD,N,f_doppler1,f_doppler,vf,lamta,ts,c,fs_JD,temp1);

        %Ŀ��ز��źš���Ŀ���źš�����������һ��������ջ�
        s_echo_1=s_echo_2+s_noise+s_ft;
        figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('�״�����ź�');
        figure,plot((0:fs_JD/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_JD-fs_JD/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�״�����źŵ�Ƶ��');
        save data/data_DeceptionJammingParameter  Gj Gjr Pj tf vf L
        
        end
        case '����'
        switch string2{temp2}
        case {'��' ,'��','�ܼ�' }
        R1=str2num(get(handles.R1,'string'));
        v1=str2num(get(handles.v1,'string'));
        congmubiao=str2num(get(handles.congmubiao,'string'));
        load data/data_JDParameter;
        fr=1/tr_JD;
        lamta=c/fz_JD;%����
        tm=0:1/fs_JD:tr_JD-1/fs_JD;%һ�������ظ����ڲ�������
        N=length(tm);%һ�������ظ����ڲ�����������
        An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%����ǿ��
        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        ts=1/fs_JD; 
        Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %Ŀ��ز��źŹ���
        A=sqrt(Prs);%�ز��źŷ���
        f_doppler=2*v/lamta;%��Ŀ�������Ƶ��
        Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(R*R)*10^(L/10)));
        Aj=sqrt(Prj);
        %�ٶ�
         f_doppler1=2*v1/lamta; 
        [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
        [~,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y,tr_JD,ts,N);
        [s_echo_2,echo]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        [s_ft,echo3]=JDDeceptionJamming(D,R1,tr_JD,c,Aj,N,frame_JD,fs_JD,f_doppler1,tm,f0_JD,tau_JD,congmubiao,temp1);
        s_echo_1=s_echo_2+s_noise+s_ft;%%%%�����ź�
        figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('�״�����ź�');
        figure,plot((0:fs_JD/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_JD-fs_JD/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�״�����źŵ�Ƶ��');
        save data/data_DeceptionJammingParameter R1 v1 congmubiao Gj Gjr Pj tf vf L
        
        case '����'
             %���ɷ����ź�
        Rj=2e3;
        c=3e8;
        lamta=c/fz_JD;%����
        tm=0:1/fs_JD:tr_JD-1/fs_JD;%һ�������ظ����ڲ�������
        N=length(tm);%һ�������ظ����ڲ�����������
        An=10*log10((1.382e-23)*Te_JD*B_JD*10^(F_JD/10));%����ǿ��
        [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        ts=1/fs_JD; 
        Prs=((Pt_JD*(10^((Gt_JD/10)))*(10^((Gr_JD/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_JD/10))); %Ŀ��ز��źŹ���
        A=sqrt(Prs);%�ز��źŷ���
        f_doppler=2*v/lamta;%��Ŀ�������Ƶ��
   
        %�ٶ�
        f_doppler=2*v/lamta;
        f_doppler1=f_doppler;
        lamta=c/fz_JD;
        f_doppler=2*v/lamta;
        [y,D]=shengchengJDxinhao(Pt_JD,tau_JD,f0_JD,tm);
        ts=1/fs_JD;%�������
        tm=0:1/fs_JD:tr_JD-1/fs_JD;%һ�������ظ����ڲ�������
        N=length(tm);%һ�������ظ����ڲ�����������
      
        Prj=((Pj*(10^((Gj/10)))*(10^((Gjr/10)))*lamta*lamta*0.5)/((4*pi*4*pi)*(Rj*Rj)*10^(L/10)));
        Aj=sqrt(Prj);
%         figure,plot(0:ts:(N-1)*ts,real(y)),xlabel('t(��λ����)'), ylabel('y(��λ����)'),title('�����ź�');
%         figure,plot((0:fs/length(abs(fft(abs(fftshift(y))))):fs-fs/length(abs(fft(abs(fftshift(y)))))),abs(fft(y))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�����źŵ�Ƶ��');
        %������������ѹ��ϵ��
        [~,match_filter_fft]=maiyaxishu(f0_JD,fs_JD,y,tr_JD,ts,N);

        %���ɻز��ź�
        [s_echo_2,echo]=JDhuiboxinhao(R,c,A,N,frame_JD,fs_JD,f_doppler,tm,f0_JD,tau_JD); 
             [s_noise]=zaosheng(frame_JD,N,An,B_JD,fs_JD);
        %���ɸ����ź�
        [s_ft,echo3]=JDtuoyinganrao(R,frame_JD,tf,Aj,tm,tau_JD,f0_JD,N,f_doppler1,f_doppler,vf,lamta,ts,c,fs_JD,temp1);

        %Ŀ��ز��źš���Ŀ���źš�����������һ��������ջ�
        s_echo_1=s_echo_2+s_noise+s_ft;
        figure,plot(0:ts:(N-1)*ts,real(s_echo_1(1:N))),xlabel('t(��λ��s)'), ylabel('y(��λ����)'),title('�״�����ź�');
        figure,plot((0:fs_JD/length(abs(fft(abs(fftshift(echo3(1,:)))))):fs_JD-fs_JD/length(abs(fft(abs(fftshift(echo3(1,:))))))),abs((fft(echo3(1,:))))),xlabel('Ƶ��f(��λ��Hz)'), ylabel('y(��λ����)'),title('�״�����źŵ�Ƶ��');
        save data/data_DeceptionJammingParameter  Gj Gjr Pj tf vf L
        
        end
     end
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
string1=get(handles.n1,'string');
string2=get(handles.n2,'string');
temp1=get(handles.n1,'Value');
temp2=get(handles.n2,'Value');
Gj=str2num(get(handles.Gj,'string'));
Gjr=str2num(get(handles.Gjr,'string')); %#ok<*NASGU,*ST2NM>
Pj=str2num(get(handles.Pj,'string'));
tf=str2num(get(handles.tf,'string'));
vf=str2num(get(handles.vf,'string'));
L=3;%�ۺ����
switch string1{temp1}
    case '�ٶ�'
        switch string2{temp2}
            case '��'
                R1=2100;
                v1=270;
                congmubiao=1;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                 
            case '��'
                R1=2100;
                v1=270;
                congmubiao=2;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                 
            case '�ܼ�'
                R1=2100;
                v1=270;
                congmubiao=10;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                 
            case '����'
                R1=2100;
                v1=270;
                tf=10;
                vf=10;
                set(handles.R1,'string','');
                set(handles.v1,'string','');
                set(handles.congmubiao,'string','');
                set(handles.R1,'Enable','off');
                set(handles.v1,'Enable','off');
                set(handles.congmubiao,'Enable','off');
                set(handles.vf,'Enable','on');
                set(handles.tf,'Enable','on');
                set(handles.tf,'string',tf);
                set(handles.vf,'string',vf);
                
        end
        
    case '����'
         switch string2{temp2}
            case '��'
                R1=2500;
                v1=300;
                congmubiao=1;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                
            case '��'
                R1=2500;
                v1=300;
                congmubiao=2;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                
            case '�ܼ�'
                R1=2500;
                v1=300;
                congmubiao=10;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                
            case '����'
                R1=2500;
                v1=300;
                tf=10;
                vf=10;
                set(handles.R1,'string','');
                set(handles.v1,'string','');
                set(handles.congmubiao,'string','');
                set(handles.R1,'Enable','off');
                set(handles.v1,'Enable','off');
                set(handles.congmubiao,'Enable','off');
                set(handles.vf,'Enable','on');
                set(handles.tf,'Enable','on');
                set(handles.tf,'string',tf);
                set(handles.vf,'string',vf);
                
         end
    case '����'
        switch string2{temp2}
            case '��'
                R1=2500;
                v1=270;
                congmubiao=1;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                
            case '��'
                R1=2500;
                v1=270;
                congmubiao=2;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                
            case '�ܼ�'
                R1=2500;
                v1=270;
                congmubiao=10;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                
            case '����'
                set(handles.R1,'string','');
                set(handles.v1,'string','');
                set(handles.congmubiao,'string','');
                set(handles.R1,'Enable','off');
                set(handles.v1,'Enable','off');
                set(handles.congmubiao,'Enable','off');
                set(handles.vf,'Enable','on');
                set(handles.tf,'Enable','on');
                tf=10;
                vf=10;
                R1=2500;
                v1=270;
                set(handles.tf,'string',tf);
                set(handles.vf,'string',vf);
                
        end
end


function Gj_Callback(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to Gj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Gj as text
%        str2double(get(hObject,'String')) returns contents of Gj as a double


% --- Executes during object creation, after setting all properties.
function Gj_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Gjr_Callback(hObject, eventdata, handles)
% hObject    handle to Gjr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Gjr as text
%        str2double(get(hObject,'String')) returns contents of Gjr as a double


% --- Executes during object creation, after setting all properties.
function Gjr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gjr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pj_Callback(hObject, eventdata, handles)
% hObject    handle to Pj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pj as text
%        str2double(get(hObject,'String')) returns contents of Pj as a double


% --- Executes during object creation, after setting all properties.
function Pj_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vf_Callback(hObject, eventdata, handles)
% hObject    handle to vf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vf as text
%        str2double(get(hObject,'String')) returns contents of vf as a double


% --- Executes during object creation, after setting all properties.
function vf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function v1_Callback(hObject, eventdata, handles)
% hObject    handle to v1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of v1 as text
%        str2double(get(hObject,'String')) returns contents of v1 as a double


% --- Executes during object creation, after setting all properties.
function v1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to v1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function R1_Callback(hObject, eventdata, handles)
% hObject    handle to R1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R1 as text
%        str2double(get(hObject,'String')) returns contents of R1 as a double


% --- Executes during object creation, after setting all properties.
function R1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function congmubiao_Callback(hObject, eventdata, handles)
% hObject    handle to congmubiao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of congmubiao as text
%        str2double(get(hObject,'String')) returns contents of congmubiao as a double


% --- Executes during object creation, after setting all properties.
function congmubiao_CreateFcn(hObject, eventdata, handles)
% hObject    handle to congmubiao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in n1.
function n1_Callback(hObject, eventdata, handles)
% hObject    handle to n1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns n1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from n1


% --- Executes during object creation, after setting all properties.
function n1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in n2.
function n2_Callback(hObject, eventdata, handles)
% hObject    handle to n2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns n2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from n2
string1=get(handles.n1,'string');
string2=get(handles.n2,'string');
temp1=get(handles.n1,'Value');
temp2=get(handles.n2,'Value');
Gj=str2num(get(handles.Gj,'string'));
Gjr=str2num(get(handles.Gjr,'string')); %#ok<*NASGU,*ST2NM>
Pj=str2num(get(handles.Pj,'string'));
tf=str2num(get(handles.tf,'string'));
vf=str2num(get(handles.vf,'string'));
L=3;%�ۺ����
switch string1{temp1}
    case '�ٶ�'
        switch string2{temp2}
            case '��'
                R1=2100;
                v1=270;
                congmubiao=1;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                 
            case '��'
                R1=2100;
                v1=270;
                congmubiao=2;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                 
            case '�ܼ�'
                R1=2100;
                v1=270;
                congmubiao=10;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                 
            case '����'
                R1=2100;
                v1=270;
                tf=10;
                vf=10;
                set(handles.R1,'string','');
                set(handles.v1,'string','');
                set(handles.congmubiao,'string','');
                set(handles.R1,'Enable','off');
                set(handles.v1,'Enable','off');
                set(handles.congmubiao,'Enable','off');
                set(handles.vf,'Enable','on');
                set(handles.tf,'Enable','on');
                set(handles.tf,'string',tf);
                set(handles.vf,'string',vf);
                
        end
        
    case '����'
         switch string2{temp2}
            case '��'
                R1=2500;
                v1=300;
                congmubiao=1;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                
            case '��'
                R1=2500;
                v1=300;
                congmubiao=2;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                
            case '�ܼ�'
                R1=2500;
                v1=300;
                congmubiao=10;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                
            case '����'
                R1=2500;
                v1=300;
                tf=10;
                vf=10;
                set(handles.R1,'string','');
                set(handles.v1,'string','');
                set(handles.congmubiao,'string','');
                set(handles.R1,'Enable','off');
                set(handles.v1,'Enable','off');
                set(handles.congmubiao,'Enable','off');
                set(handles.vf,'Enable','on');
                set(handles.tf,'Enable','on');
                set(handles.tf,'string',tf);
                set(handles.vf,'string',vf);
                
         end
    case '����'
        switch string2{temp2}
            case '��'
                R1=2500;
                v1=270;
                congmubiao=1;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                
            case '��'
                R1=2500;
                v1=270;
                congmubiao=2;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                
            case '�ܼ�'
                R1=2500;
                v1=270;
                congmubiao=10;
                set(handles.R1,'Enable','on');
                set(handles.v1,'Enable','on');
                set(handles.congmubiao,'Enable','on');
                set(handles.R1,'string',R1);
                set(handles.v1,'string',v1);
                set(handles.congmubiao,'string',congmubiao);
                set(handles.vf,'string','');
                set(handles.tf,'string','');
                set(handles.vf,'Enable','off');
                set(handles.tf,'Enable','off');
                
            case '����'
                set(handles.R1,'string','');
                set(handles.v1,'string','');
                set(handles.congmubiao,'string','');
                set(handles.R1,'Enable','off');
                set(handles.v1,'Enable','off');
                set(handles.congmubiao,'Enable','off');
                set(handles.vf,'Enable','on');
                set(handles.tf,'Enable','on');
                tf=10;
                vf=10;
                R1=2500;
                v1=270;
                set(handles.tf,'string',tf);
                set(handles.vf,'string',vf);
                
        end
end

% --- Executes during object creation, after setting all properties.
function n2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tf_Callback(hObject, eventdata, handles)
% hObject    handle to tf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tf as text
%        str2double(get(hObject,'String')) returns contents of tf as a double


% --- Executes during object creation, after setting all properties.
function tf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
