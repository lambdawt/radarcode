function s=nonline_eq_sirp(z,nu)
%
for m=1:length(z)
    max_s=10;
    min_s=0;
    isexistmax=0;
    F_z=0.5+0.5*erf(z(m)/sqrt(2));
    if (F_z>0.9999999)
        s(m)=inf;
    else
        for k=1:10000
            if (gammainc(nu*max_s*max_s,nu)<F_z)
                min_s=max_s;
                max_s=max_s*2;
            else
                isexistmax=1;break;
            end
        end
        if (isexistmax==1)
            s(m)=0.5*(max_s-min_s);
            s0=0;
            for k=1:10000
                if (gammainc(nu*s(m)*s(m),nu)>F_z)
                    max_s=s(m);
                else
                    min_s=s(m);
                end
                s0=s(m);
                s(m)=0.5*(max_s+min_s);
                if (abs(s(m)-s0)<0.0001)
                    break;
                end
            end
        end
    end
end
figure;plot(abs(s));