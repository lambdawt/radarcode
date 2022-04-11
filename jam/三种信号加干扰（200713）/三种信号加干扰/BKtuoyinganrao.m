function [s_ft,echo3]=BKtuoyinganrao(R,frame,tf,Aj,tau,N,f_doppler1,f_doppler,vf,lamta,ts,c,fs,y1,tr,temp1)
CPI=0;
P=0;
if temp1 == 3%联合
    for i=0:frame-1
        if mod(i,16)==0
           P=tf*(i+1)*ts;%(c/f0/100)
        else
        end
        %echo3(i+1,:)=Aj*rectpuls(tm-delt1-P-tau/2,tau).*exp(j*2*pi*(((f0-B1/2)+k*(tm-delt1-P)/2).*(tm-delt1-P)));
       echo3(i+1,:)= Aj.*[zeros(1,floor((2*R/c+P)/ts)),y1,zeros(1,ceil((tr-(2*R/c+P)-tau)/ts))];
    end
    for i=1:frame
        s_ft1_1(1,(i-1)*N+1:i*N)=echo3(i,:);%将frame行回波信号重复脉冲接成一行
    end

    s_ft=zeros(1,frame*N);%生成1行，frame*N列个0矩阵，便于产生多个假目标用
    s_doppler_1(1,1:frame*N)=exp(1j*(1:frame*N)*2*pi*f_doppler(1,:)/fs);
    s_ft1_1(1,CPI*16*N+1:(CPI+1)*16*N)=s_ft1_1(1,CPI*16*N+1:(CPI+1)*16*N).*s_doppler_1(1,CPI*16*N+1:(CPI+1)*16*N);
    for i=0:frame-1
        if   mod(i,16)==0 && CPI<frame/16-2
             CPI=CPI+1;
             if CPI==1
                 f_doppler1(1,:)=f_doppler1-2*(vf/lamta);

             else
                 f_doppler1(CPI,1)=f_doppler1(CPI-1,1)-2*(vf/lamta);%*((i-s_d_t_end)*N*ts);%(20*deltf)*i*N*ts;
             end
    %          (20*deltf)*i*N*ts%*(c/f0/100)/2;
             s_doppler_1(CPI,1:frame*N)=exp(1j*(1:frame*N)*2*pi*f_doppler1(CPI,:)/fs);
             s_ft1_1(1,(CPI-1)*16*N+1:CPI*16*N)=s_ft1_1(1,(CPI-1)*16*N+1:CPI*16*N).*s_doppler_1(CPI,(CPI-1)*16*N+1:CPI*16*N);


        else 
        end
    end
    s_ft=s_ft+s_ft1_1;
   
elseif temp1 == 1%速度
    for i=0:frame-1
        if mod(i,16)==0
           P=0;%(c/f0/100)
        else
        end
        %echo3(i+1,:)=Aj*rectpuls(tm-delt1-P-tau/2,tau).*exp(j*2*pi*(((f0-B1/2)+k*(tm-delt1-P)/2).*(tm-delt1-P)));
       echo3(i+1,:)= Aj.*[zeros(1,floor((2*R/c+P)/ts)),y1,zeros(1,ceil((tr-(2*R/c+P)-tau)/ts))];
    end
    for i=1:frame
        s_ft1_1(1,(i-1)*N+1:i*N)=echo3(i,:);%将frame行回波信号重复脉冲接成一行
    end

    s_ft=zeros(1,frame*N);%生成1行，frame*N列个0矩阵，便于产生多个假目标用
    s_doppler_1(1,1:frame*N)=exp(1j*(1:frame*N)*2*pi*f_doppler(1,:)/fs);
    s_ft1_1(1,CPI*16*N+1:(CPI+1)*16*N)=s_ft1_1(1,CPI*16*N+1:(CPI+1)*16*N).*s_doppler_1(1,CPI*16*N+1:(CPI+1)*16*N);
    for i=0:frame-1
        if   mod(i,16)==0 && CPI<frame/16-2
             CPI=CPI+1;
             if CPI==1
                 f_doppler1(1,:)=f_doppler1-2*(vf/lamta);

             else
                 f_doppler1(CPI,1)=f_doppler1(CPI-1,1)-2*(vf/lamta);%*((i-s_d_t_end)*N*ts);%(20*deltf)*i*N*ts;
             end
    %          (20*deltf)*i*N*ts%*(c/f0/100)/2;
             s_doppler_1(CPI,1:frame*N)=exp(1j*(1:frame*N)*2*pi*f_doppler1(CPI,:)/fs);
             s_ft1_1(1,(CPI-1)*16*N+1:CPI*16*N)=s_ft1_1(1,(CPI-1)*16*N+1:CPI*16*N).*s_doppler_1(CPI,(CPI-1)*16*N+1:CPI*16*N);


        else 
        end
    end
    s_ft=s_ft+s_ft1_1;
    
elseif temp1 == 2%距离
    for i=0:frame-1
        if mod(i,16)==0
           P=tf*(i+1)*ts;%(c/f0/100)
        else
        end
        %echo3(i+1,:)=Aj*rectpuls(tm-delt1-P-tau/2,tau).*exp(j*2*pi*(((f0-B1/2)+k*(tm-delt1-P)/2).*(tm-delt1-P)));
       echo3(i+1,:)= Aj.*[zeros(1,floor((2*R/c+P)/ts)),y1,zeros(1,ceil((tr-(2*R/c+P)-tau)/ts))];
    end
    for i=1:frame
        s_ft1_1(1,(i-1)*N+1:i*N)=echo3(i,:);%将frame行回波信号重复脉冲接成一行
    end

    s_ft=zeros(1,frame*N);%生成1行，frame*N列个0矩阵，便于产生多个假目标用
    s_doppler_1(1,1:frame*N)=exp(1j*(1:frame*N)*2*pi*f_doppler(1,:)/fs);
    s_ft1_1(1,CPI*16*N+1:(CPI+1)*16*N)=s_ft1_1(1,CPI*16*N+1:(CPI+1)*16*N).*s_doppler_1(1,CPI*16*N+1:(CPI+1)*16*N);
    for i=0:frame-1
        if   mod(i,16)==0 && CPI<frame/16-2
             CPI=CPI+1;
             if CPI==1
                 f_doppler1(1,:)=f_doppler1;

             else
                 f_doppler1(CPI,1)=f_doppler1(CPI-1,1);%*((i-s_d_t_end)*N*ts);%(20*deltf)*i*N*ts;
             end
    %          (20*deltf)*i*N*ts%*(c/f0/100)/2;
             s_doppler_1(CPI,1:frame*N)=exp(1j*(1:frame*N)*2*pi*f_doppler1(CPI,:)/fs);
             s_ft1_1(1,(CPI-1)*16*N+1:CPI*16*N)=s_ft1_1(1,(CPI-1)*16*N+1:CPI*16*N).*s_doppler_1(CPI,(CPI-1)*16*N+1:CPI*16*N);


        else 
        end
    end
    s_ft=s_ft+s_ft1_1;    
    
end