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