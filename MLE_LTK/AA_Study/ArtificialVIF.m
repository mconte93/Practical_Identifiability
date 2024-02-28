%Code to generate insilico VIF
clear all
close all

t=linspace(0,300,50);
A1=0.309*60;
A2=0.330*60;
sigma1=0.0563*60;
sigma2=0.132*60;
T1=0.17046*60;
T2=0.365*60;
alfa=1.050;
beta=0.1685/60;
s=38.078/60;
tau=0.483*60;

sum1=A1./(sigma1.*sqrt(2.*pi)).*exp(-(t-T1).^2./(2.*sigma1.^2));
sum2=A2./(sigma2.*sqrt(2.*pi)).*exp(-(t-T2).^2./(2.*sigma2.^2));

VIF=(sum1+sum2+alfa.*exp(-beta*t)./(1+exp(-s.*(t-tau))))./5;

%Saving data option 
%save('ArtificialVIF.mat','VIF','t')