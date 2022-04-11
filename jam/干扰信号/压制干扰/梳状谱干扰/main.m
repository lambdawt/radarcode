load data/data_shuzhuangpu
fj=[1e7,4e7,7e7];
[sig_noise,t_noise] = jam_combspectrum(fs_shuzhuangpu,Bj_shuzhuangpu,Ns_shuzhuangpu,fj,frame_shuzhuangpu,Prj_shuzhuangpu,Tr_shuzhuangpu);
view_jam_combspectrum( sig_noise,t_noise,fs_shuzhuangpu );