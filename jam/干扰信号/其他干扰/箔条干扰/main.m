load data/data_botiao
 rcsk=0;
            sigma = rcs(rcsk,sigma0);
        fs=2e6;f0=0.5e6;
        R=2e3;c=3e8;A=5;Tr=600e-6;tau=100e-6;frame=1;tm=0:1/fs:Tr-1/fs;N=length(tm);f_doppler=1e4;
        [s_echo_2]=JDhuiboxinhao(R,c,A,N,frame,fs,f_doppler,tm,f0,tau);

        %Ŀ��״̬����
        px=1e3;py=1e3;pz=1e3;%Ŀ��λ��
        vx=10;vy=10;vz=0;%Ŀ���ٶ�
        ax=0;ay=0;az=0;%Ŀ����ٶ�
        phi=pi/180;
        [ TargetStatus ] = paraset_targetstatus( px,py,pz,vx,vy,vz,ax,ay,az,phi );

        [ PassivePara ] = paraset_passivejaming( tf_botiao,sf_botiao,vl_botiao,vf_botiao,ts_botiao,bt_botiao,al_botiao,sref_botiao,smax_botiao );

          vx=1;vy=1;
        [ WindV ] = paraset_windvelocity( vx,vy );


        echo=s_echo_2;fc=f0;CurrentT=2.2;
        [ sig_jam,t_jam ] = jam_passive( echo,fs,fc,CurrentT,TargetStatus,WindV,PassivePara );

        figure;
        subplot(211)
        plot(t_jam,real(sig_jam));grid;xlabel('ʱ��(s)');ylabel('����(V)');title('���������ź�ʱ����'); 
        subplot(212)
        Fsig=fftshift(fft(sig_jam));
        K=length(t_jam);
        k=floor(-K/2+0.5:K/2-0.5);
        dfs=fs/K;
        plot(k*dfs,abs(Fsig));grid;xlabel('Ƶ��(Hz)');ylabel('����(V/Hz)');title('���������ź�Ƶ����');