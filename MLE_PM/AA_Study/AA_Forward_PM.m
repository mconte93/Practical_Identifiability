%Code for the forward PM model to generate insilico data with insilico VIF

clear all
close all

ktr=0.002;
vp=0.45;

%Parameter vector
par=[ktr,vp];

%Initial condition
C_0=0;  
Init=C_0;

%Load VIF produced with the code 'ArtificialVIF.m'
load('ArtificialVIF.mat')
time=t;

%Forward TPM model
[t,y] = ode45(@(t,y) PM_ODE(t,y,par,VIF,time),time,Init,[]);
C_for=y+vp*VIF';
t_for=t;

%Saving data option 
%save ('AA_PM_Persistent.mat','t_for','C_for','par','VIF')
