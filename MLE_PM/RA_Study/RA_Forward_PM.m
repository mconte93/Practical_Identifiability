%Code for the forward PM model to generate insilico data with real VIF

clear all
close all

ktr=0.0015;
vp=0.01;

%Parameter vector
par=[ktr,vp];

%Initial conditions
C_0=0;  
Init=C_0;

%Load real VIF estimated from the data
load('VIF.mat')
time=time_VIF;
VIF=VIF';

%Forward PM model
[t,y] = ode45(@(t,y) PM_ODE(t,y,par,VIF,time),time,Init,[]);
C_for=y+vp*VIF;
t_for=t;

%Saving data option
%save ('RA_PM_Persistent.mat','t_for','C_for','par','VIF')

