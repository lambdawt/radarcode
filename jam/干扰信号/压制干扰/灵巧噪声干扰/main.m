load data/data_smart
load data/data_LFMParameter;
        load data/data_target0Parameter
Kfm=4e6;tau=1e-6; c=3e8;
rcsk=0;%起伏类型
sigma = rcs(rcsk,sigma0);
lamta=c/fz_LFM;%波长 
        k=B1_LFM/tau_LFM;                                 %线性调频信号调制系数
        tm=0:1/fs_LFM:tr_LFM-1/fs_LFM;  
        N=length(tm);
        f_doppler=2*v/lamta;%真目标多普勒频率
        Prs=((Pt_LFM*(10^((Gt_LFM/10)))*(10^((Gr_LFM/10)))*lamta*lamta*sigma)/((4*pi*4*pi*4*pi)*(R*R*R*R)*10^(L_LFM/10))); %目标回波信号功率
        A=sqrt(Prs);%回波信号幅度
        [vRadarSig]=LFMhuiboxinhao(R,c,A,N,frame_LFM,fs_LFM,f_doppler,tm,f0_LFM,B1_LFM,tau_LFM,k); 
        Pn=(B1_LFM/(2.5*Kfm))^2;Bn=B1_LFM/2;
        [vSmartNoiseSig]=jam_smartnoise( vRadarSig,Pn,Prj_smart,Bn,Kfm,fs_smart );
        view_SmartNoise( vSmartNoiseSig,fs_smart,tau,R); 
        