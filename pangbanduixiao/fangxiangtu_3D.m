function[pattern]=rect_array(Ny,Nz,d_lamda,theta0,fai0,theta)
%The calculation of three plane antenna pattern
d_lamda=0.5;
fai0=0;%azimuth
theta0=0;%Pitching Angle
theta=-90:1:90;
Ny=8;
Nz=8;
ay0=exp(j*2*pi*d_lamda*(0:Ny-1)'*sin(fai0*pi/180)*cos(theta0*pi/180));%y airspace direction vector
az0=exp(j*2*pi*d_lamda*(0:Nz-1)'*sin(theta0*pi/180));%z airspace direction vector
aa0=kron(ay0,az0);
for m=1:length(theta)
    for n=1:length(theta)
        ay=exp(j*2*pi*d_lamda*(0:Ny-1)'*sin(theta(m)*pi/180)*cos(theta(n)*pi/180));%y airspace direction vector
        az=exp(j*2*pi*d_lamda*(0:Nz-1)'*sin(theta(n)*pi/180));%z airspace direction vector
        aa=kron(ay,az);
        pattern(m,n)=abs(aa0'*aa);
    end
end
figure;mesh(theta,theta,pattern);
xlabel('azimuth\circ');ylabel('Pitching Angle\circ');zlabel('The normalized pattern')
figure;
[c,h]=contour(theta,theta,pattern,[3.5,8.5,14:5:64]);%Draw a rectangular planar contours
xlabel('azimuthcirc');ylabel('Pitching Angle\circ');