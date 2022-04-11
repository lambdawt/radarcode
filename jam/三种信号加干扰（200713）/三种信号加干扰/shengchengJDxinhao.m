function [y,D]=shengchengJDxinhao(Pt,tau,f0,tm)

D=5;
y=sqrt(Pt)*rectpuls(tm-tau/2,tau).*exp(1j*2*pi*(f0.*tm));