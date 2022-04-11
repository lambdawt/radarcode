 load data/data_JDParameter %#ok<*LOAD>
            load data/data_target0Parameter
            load data/data_dopplerblink
            rcsk=0;
            sigma = rcs(rcsk,sigma0);
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