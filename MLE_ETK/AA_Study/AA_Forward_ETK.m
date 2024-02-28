%Code for the forward ETK model to generate insilico data with insilico VIF

clear all
close all

ktr=0.015;
ve=0.5;
vp=0.2;

%Parameter vector
par=[ktr,ve,vp];

%Initial conditions
C_0=0;  
Init=C_0;

%Load VIF produced with the code 'ArtificialVIF.m'
load('ArtificialVIF.mat')
time=t;

%Forward ETK model
[t,y] = ode45(@(t,y) ETK_ODE(t,y,par,VIF,time),time,Init,[]);
C_for=y+vp*VIF';
t_for=t;


%Saving data option 
%save ('AA_ETK_Persistent.mat','t_for','C_for','par','VIF')
