load data/data_saopin
 %T_fr_saopin=2*Tr_saopin
 [ sig_noise,t_noise ] = jam_sweepfrequency( fs_saopin,Bj_saopin,fj_saopin,frame_saopin,Prj_saopin,Tr_saopin,T_fr_saopin,Time_begin_saopin,K_sweep_saopin );
view_jam_sweepfrequency( sig_noise,t_noise,fs_saopin );