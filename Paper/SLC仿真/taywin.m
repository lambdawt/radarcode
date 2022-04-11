function tay=taywin(n,SLL,N);

% 1 n is the number of coefficients in the Taylor windows
% 2 N is the number of nearly constant-level sidelobes adjacent to mainlobe;
% 2 SLL is the peak sidelobe level (in dB) relative to the mainlobe peak;
% tay=taywin(256,-35,5)

B=10.^(-SLL/20);
A=log(B+sqrt(B*B-1))/pi;
cgma2=N.^2/(A*A+(N-1/2)*(N-1/2));

for m=1:N-1
   NUM=1;
   DEN=1;
   for i=1:N-1
      NUM=NUM*(1-m*m/cgma2/(A*A+(i-1/2)*(i-1/2)));
      if (i~=m)
         DEN=DEN*(1-m*m/(i*i));
      end
   end
   F(m)=(-1).^(m+1)*NUM/(2*DEN)
end
for k=1:n
  tay(k)=1;
  for m=1:N-1
     tay(k)=tay(k)+2*F(m)*cos(2*pi*m*(k-n/2-1+1/2)/n);%-1changed by wzp,add   to (k-n/2+1/2)
  end
end
tay=tay.';